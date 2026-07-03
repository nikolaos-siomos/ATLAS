#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 16:20:32 2026

@author: nikos
"""

import re
import os
import glob
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Union, Dict, Any, Optional
from utils.error_classes import CustomWarning

WYOMING_RE = re.compile(
    # Supports both legacy Wyoming files and downloaded BUFR text files, e.g.
    #   20240428_1200_wyoming_16716.dat
    #   20240428_1200_wyoming_bufr_16716.txt
    r"^(?P<date>\d{8})_(?P<time>\d{4})_wyoming(?:_bufr)?_(?P<identifier>\d{5})\.(?:dat|txt)$"
)

ECMWF_RE = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{4})_ecmwf_(?P<identifier>[A-Za-z0-9][A-Za-z0-9_-]*)\.nc$"
)

SCC_RE = re.compile(
    r"^rs_(?P<date>\d{8})(?P<identifier>[A-Za-z]+?)(?P<time>\d{2})\.nc$"
)

ASCII_RE = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{4})_custom_ascii_(?P<identifier>[A-Za-z0-9][A-Za-z0-9_-]*)\.(?:dat|txt|csv|asc)$"
)

RADIOSONDE_PATTERNS = {
    "wyoming": (WYOMING_RE, "%Y%m%d%H%M"),
    "ecmwf": (ECMWF_RE, "%Y%m%d%H%M"),
    "scc": (SCC_RE, "%Y%m%d%H"),
    "custom_ascii": (ASCII_RE, "%Y%m%d%H%M"),
}


def _normalise_identifier(identifier: Any) -> Optional[str]:
    """
    Convert identifier values from caller_info / metadata to comparable strings.

    Handles plain strings, numbers, numpy scalars, and scalar xarray DataArrays
    without importing xarray here.
    """

    if identifier is None:
        return None

    if hasattr(identifier, "values"):
        identifier = identifier.values

    if hasattr(identifier, "item"):
        try:
            identifier = identifier.item()
        except ValueError:
            return None

    identifier = str(identifier).strip()

    if identifier == "" or identifier.lower() == "nan":
        return None

    # Match downloaded/generated filenames robustly.
    # Cloudnet downloads use lowercase station identifiers, e.g. Barcelona -> barcelona.
    identifier = identifier.lower()

    return identifier


def _identifier_matches(
    match: re.Match,
    file_format: str,
    identifiers: Optional[Dict[str, Any]],
) -> bool:
    """
    Check whether the identifier parsed from the filename matches the expected
    identifier for this radiosonde format.

    If identifiers is None, no identifier filtering is applied. This keeps
    manual mode backward-compatible.
    """

    if identifiers is None:
        return True

    filename_identifier = _normalise_identifier(match.group("identifier"))
    expected_identifier = _normalise_identifier(identifiers.get(file_format))

    if expected_identifier is None:
        return False

    return filename_identifier == expected_identifier


def parse_radiosonde_filename(
    filename: Union[str, Path],
    identifiers: Optional[Dict[str, Any]] = None,
):
    """
    Parse a supported radiosonde filename.

    Returns
    -------
    parsed, status:
        parsed is a dictionary with keys:
          - radiosonde_format
          - radiosonde_time
          - radiosonde_identifier

        status is True if the filename matched a supported ATLAS radiosonde
        naming pattern and, when identifiers are provided, the filename
        identifier matches the expected identifier for that format. Otherwise
        status is False.
    """

    name = os.path.basename(str(filename))

    for file_format, (pattern, strip_format) in RADIOSONDE_PATTERNS.items():
        match = pattern.fullmatch(name)

        if match is None:
            continue

        if not _identifier_matches(match, file_format, identifiers):
            continue

        date = match.group("date")
        time = match.group("time")
        dt = datetime.strptime(date + time, strip_format)

        parsed = {
            "radiosonde_format": file_format,
            "radiosonde_time": np.datetime64(dt, "ns"),
            "radiosonde_identifier": _normalise_identifier(match.group("identifier")),
        }

        return parsed, True

    return {}, False


def infer_manual_radiosonde_format(filename: Union[str, Path]):
    """
    Infer the loader format from the standard ATLAS filename patterns only.

    This is kept as a lightweight helper/backward-compatible fallback. For
    manually provided files, prefer `detect_manual_radiosonde_format`, which
    tries the actual readers and therefore does not depend on the filename.
    """

    parsed, status = parse_radiosonde_filename(filename)

    if status:
        return parsed["radiosonde_format"]

    return "unknown"


def _as_float_or_default(value: Any, default: float = 0.0) -> float:
    """Convert scalar metadata values to float, falling back to default."""

    if value is None:
        return default

    if hasattr(value, "values"):
        value = value.values

    if hasattr(value, "item"):
        try:
            value = value.item()
        except ValueError:
            return default

    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _make_radiosonde_info_da(
    radiosonde_file: Union[str, Path],
    radiosonde_format: str,
    radiosonde_time: object,
    radiosonde_source: str = "manual",
):
    """Create the minimal radiosonde_info DataArray expected by readers."""

    import xarray as xr

    radiosonde_info = {
        "radiosonde_file": os.path.normpath(str(radiosonde_file)),
        "radiosonde_format": radiosonde_format,
        "radiosonde_source": radiosonde_source,
        "radiosonde_time": radiosonde_time,
        "measurement_time": radiosonde_time,
        "radiosonde_status": 0,
    }

    return xr.DataArray(
        data=list(radiosonde_info.values()),
        dims=["parameters"],
        coords={"parameters": list(radiosonde_info.keys())},
    )


def _manual_reader_order(filename: Union[str, Path]):
    """Return a practical reader order for probing a manual radiosonde file."""

    parsed, status = parse_radiosonde_filename(filename)
    preferred = [parsed["radiosonde_format"]] if status else []

    suffix = Path(filename).suffix.lower()

    if suffix == ".nc":
        candidates = ["ecmwf", "scc", "wyoming", "custom_ascii"]
    else:
        candidates = ["custom_ascii", "wyoming", "ecmwf", "scc"]

    return preferred + [fmt for fmt in candidates if fmt not in preferred]


def _reader_probe_succeeded(meteo: Any) -> bool:
    """Basic sanity check for the object returned by a radiosonde reader."""

    if meteo is None:
        return False

    if not hasattr(meteo, "dims"):
        return False

    if "height_asl" not in meteo.dims:
        return False

    if "atmo_parameters" not in meteo.dims:
        return False

    if meteo.sizes.get("height_asl", 0) == 0:
        return False

    return True


def detect_manual_radiosonde_format(
    target: object,
    radiosonde_file: Union[str, Path],
    caller_info: Optional[Dict[str, Any]] = None,
    station_altitude: Any = None,
):
    """
    Detect the format of a manually provided radiosonde by trying readers.

    This avoids relying on the manual filename. Each candidate format is placed
    into a temporary radiosonde_info object and the corresponding reader is
    called. The first reader that returns a valid meteo DataArray determines the
    manual `radiosonde_format`. The actual profile is intentionally discarded;
    downstream code still reads the selected file through `load_radiosonde`.
    """

    from utils.read_radiosondes import (
        read_radiosonde_ecmwf,
        read_radiosonde_wyoming,
        read_radiosonde_scc,
        read_radiosonde_ascii,
    )

    station_altitude = _as_float_or_default(station_altitude, default=0.0)
    last_error = None

    readers = {
        "ecmwf": lambda info: read_radiosonde_ecmwf(station_altitude, info),
        "wyoming": read_radiosonde_wyoming,
        "scc": read_radiosonde_scc,
        "custom_ascii": lambda info: read_radiosonde_ascii(caller_info, info),
    }

    for radiosonde_format in _manual_reader_order(radiosonde_file):
        if radiosonde_format == "custom_ascii" and caller_info is None:
            continue

        radiosonde_info = _make_radiosonde_info_da(
            radiosonde_file=radiosonde_file,
            radiosonde_format=radiosonde_format,
            radiosonde_time=target,
            radiosonde_source="manual",
        )

        try:
            meteo = readers[radiosonde_format](radiosonde_info)
        except Exception as exc:
            last_error = exc
            continue

        if _reader_probe_succeeded(meteo):
            return radiosonde_format, True, None

    return "unknown", False, last_error


def select_manual_radiosonde_file(
    target: object,
    radiosonde_file: Union[str, Path],
    caller_info: Optional[Dict[str, Any]] = None,
    station_altitude: Any = None,
):
    """
    Use a manually provided radiosonde file without folder search or downloads.

    Manual mode does not reject files based on a time window. The file format is
    detected by trying the available radiosonde readers. If the timestamp can be
    parsed from the filename, it is reported. Otherwise the measurement midpoint
    is used as a metadata placeholder.
    """

    radiosonde_file = os.path.normpath(str(radiosonde_file))

    if not os.path.isfile(radiosonde_file):
        CustomWarning(
            f"The manually provided radiosonde_file does not point to an existing file: {radiosonde_file}"
        )
        print()
        return {}, 2

    parsed, name_status = parse_radiosonde_filename(radiosonde_file)
    radiosonde_time = parsed["radiosonde_time"] if name_status else target

    radiosonde_format, read_status, last_error = detect_manual_radiosonde_format(
        target=radiosonde_time,
        radiosonde_file=radiosonde_file,
        caller_info=caller_info,
        station_altitude=station_altitude,
    )

    if not read_status:
        msg = (
            "The manually provided radiosonde file could not be parsed by any "
            "available radiosonde reader."
        )

        if last_error is not None:
            msg += f" Last reader error: {last_error}"

        CustomWarning(msg)
        print()
        return {}, 2

    if not name_status:
        CustomWarning(
            "The manually provided radiosonde filename does not match the standard "
            "ATLAS timestamp patterns. The measurement midpoint will be stored as "
            "radiosonde_time."
        )
        print()

    output = {
        "radiosonde_file": radiosonde_file,
        "radiosonde_format": radiosonde_format,
        "radiosonde_source": "manual",
        "radiosonde_time": radiosonde_time,
    }

    return output, 0

def select_radiosonde_filename(
    target: object,
    folder: Union[str, Path],
    identifiers: Dict[str, Any],
    time_limit: int = 18,
    priority_time_limit: int = 3,
) -> str:
    """
    Find the best radiosonde file in `folder` using filename timestamps.

    Only files matching these filename formats are considered:
      - Wyoming:      <yyyymmdd>_<hhmm>_wyoming_<wmo_id>.dat
      - Wyoming BUFR: <yyyymmdd>_<hhmm>_wyoming_bufr_<wmo_id>.txt
      - ECMWF:        <yyyymmdd>_<hhmm>_ecmwf_<site>.nc
      - SCC:          rs_<yyyymmdd><site><hh>.nc
      - Custom ASCII: <yyyymmdd>_<hhmm>_custom_ascii_<site>.<txt|dat|csv|asc>

    In automatic mode, the identifier parsed from the filename must match the
    identifier provided in the `identifiers` dictionary for the corresponding
    radiosonde format. The expected keys are:
      - "wyoming"
      - "ecmwf"
      - "scc"
      - "custom_ascii"

    Selection rule:
      1. Prefer the closest Wyoming file within +/- `priority_time_limit` hours.
      2. Otherwise, select the closest supported file within +/- `time_limit` hours.
      3. Ignore unrelated files, folders, and files with mismatching identifiers.

    Parameters
    ----------
    target:
        Target time, usually numpy.datetime64.
    folder:
        Folder containing radiosonde files.
    identifiers:
        Unique identifiers related to the lidar system to link with radiosondes
        which correspond to the same station.
    time_limit:
        Fallback time window in hours. Default: 18.
    priority_time_limit:
        Wyoming priority window in hours. Default: 3.

    Returns
    -------
    output, status:
        output is a dictionary with keys:
          - "radiosonde_file": selected file path
          - "radiosonde_time": datetime parsed from filename
          - "radiosonde_format": "wyoming", "ecmwf", "scc", or "custom_ascii"
          - "radiosonde_source": "auto"

        status codes:
          - 0: file found
          - 1: supported files exist, but none within `time_limit`
          - 2: no supported files found

        If status is 1 or 2, output is an empty dictionary.
    """

    filenames = glob.glob(os.path.join(folder, "*"))

    valid_name_ind = []
    time_stamps = []
    delta_hours = []
    file_formats = []

    output = {}

    for i, filename in enumerate(filenames):
        if not os.path.isfile(filename):
            continue

        parsed, status = parse_radiosonde_filename(filename, identifiers)

        if not status:
            continue

        valid_name_ind.append(i)

        radiosonde_time = parsed["radiosonde_time"]
        radiosonde_format = parsed["radiosonde_format"]

        time_stamps.append(radiosonde_time)
        delta_hours.append((radiosonde_time - target) / np.timedelta64(1, "h"))
        file_formats.append(radiosonde_format)

    if len(valid_name_ind) > 0:
        file_formats = np.array(file_formats)
        delta_hours = np.array(delta_hours)
        valid_name_ind = np.array(valid_name_ind)

        mask_wyoming = (file_formats == "wyoming") & \
            (np.abs(delta_hours) <= priority_time_limit)

        mask_time = (np.abs(delta_hours) <= time_limit)

        if mask_wyoming.any():
            delta_hours_wyoming = delta_hours[mask_wyoming]
            candidate_indices = np.where(mask_wyoming)[0]
            min_time_ind = candidate_indices[np.argmin(np.abs(delta_hours_wyoming))]

            target_filename = filenames[valid_name_ind[min_time_ind]]
            target_timestamp = time_stamps[min_time_ind]
            target_format = file_formats[min_time_ind]

            output = {
                "radiosonde_file": target_filename,
                "radiosonde_format": target_format,
                "radiosonde_source": "auto",
                "radiosonde_time": target_timestamp,
            }

            status = 0

        elif mask_time.any():
            delta_hours_time = delta_hours[mask_time]
            candidate_indices = np.where(mask_time)[0]
            min_time_ind = candidate_indices[np.argmin(np.abs(delta_hours_time))]

            target_filename = filenames[valid_name_ind[min_time_ind]]
            target_timestamp = time_stamps[min_time_ind]
            target_format = file_formats[min_time_ind]

            output = {
                "radiosonde_file": target_filename,
                "radiosonde_format": target_format,
                "radiosonde_source": "auto",
                "radiosonde_time": target_timestamp,
            }

            status = 0

        else:
            status = 1

            CustomWarning(f"No valid radiosonde file detected within {time_limit} hours from the middle of the measurement inside: {folder}")
            print()

    else:
        status = 2

        CustomWarning(f"No valid radiosonde file detected inside: {folder}")
        print()

    return output, status
