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
from typing import Union
from pathlib import Path
from datetime import datetime
from utils.error_classes import CustomWarning

WYOMING_RE = re.compile(
    # Supports both legacy Wyoming files and downloaded BUFR text files, e.g.
    #   20240428_1200_wyoming_16716.dat
    #   20240428_1200_wyoming_bufr_16716.txt
    r"^(?P<date>\d{8})_(?P<time>\d{4})_wyoming(?:_bufr)?_(?P<wmo_id>\d{5})\.(?:dat|txt)$"
)

ECMWF_RE = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{4})_ecmwf_(?P<site>[A-Za-z0-9][A-Za-z0-9_-]*)\.nc$"
)

SCC_RE = re.compile(
    r"^rs_(?P<date>\d{8})(?P<site>[A-Za-z]+?)(?P<time>\d{2})\.nc$"
)

ASCII_RE = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{4}).*\.txt$"
)


def parse_radiosonde_filename(filename: Union[str, Path]):
    """
    Parse a supported radiosonde filename.

    Returns
    -------
    parsed, status:
        parsed is a dictionary with keys:
          - radiosonde_format
          - radiosonde_time

        status is True if the filename matched a supported ATLAS radiosonde
        naming pattern, otherwise False.
    """

    name = os.path.basename(str(filename))

    wyoming_match = WYOMING_RE.fullmatch(name)
    ecmwf_match = ECMWF_RE.fullmatch(name)
    scc_match = SCC_RE.fullmatch(name)
    ascii_match = ASCII_RE.fullmatch(name)

    if wyoming_match:
        match = wyoming_match
        strip_format = "%Y%m%d%H%M"
        file_format = "wyoming"

    elif ecmwf_match:
        match = ecmwf_match
        strip_format = "%Y%m%d%H%M"
        file_format = "ecmwf"

    elif scc_match:
        match = scc_match
        strip_format = "%Y%m%d%H"
        file_format = "scc"

    elif ascii_match:
        match = ascii_match
        strip_format = "%Y%m%d%H%M"
        file_format = "ascii"

    else:
        return {}, False

    date = match.group("date")
    time = match.group("time")
    dt = datetime.strptime(date + time, strip_format)

    parsed = {
        "radiosonde_format": file_format,
        "radiosonde_time": np.datetime64(dt, "ns"),
    }

    return parsed, True


def infer_manual_radiosonde_format(filename: Union[str, Path]):
    """
    Infer the loader format for a manually provided radiosonde file.

    First uses the standard ATLAS filename patterns. If the filename does not
    contain a parseable timestamp, fall back to the extension for common manual
    ASCII files.
    """

    parsed, status = parse_radiosonde_filename(filename)

    if status:
        return parsed["radiosonde_format"]

    suffix = Path(filename).suffix.lower()

    if suffix in [".txt", ".dat", ".csv", ".asc"]:
        return "ascii"

    if suffix == ".nc":
        # Unknown NetCDF radiosonde. Keep it explicit rather than pretending it
        # is SCC or ECMWF. The downstream loader can decide whether to support it.
        return "unknown"

    return "unknown"


def select_manual_radiosonde_file(
    target: object,
    radiosonde_file: Union[str, Path],
):
    """
    Use a manually provided radiosonde file without folder search or downloads.

    Manual mode does not reject files based on a time window. If the timestamp
    can be parsed from the filename, it is reported. Otherwise the measurement
    midpoint is used as a metadata placeholder.
    """

    radiosonde_file = os.path.normpath(str(radiosonde_file))

    if not os.path.isfile(radiosonde_file):
        CustomWarning(
            f"The manually provided radiosonde_file does not point to an existing file: {radiosonde_file}"
        )
        print()
        return {}, 2

    parsed, status = parse_radiosonde_filename(radiosonde_file)

    if status:
        radiosonde_format = parsed["radiosonde_format"]
        radiosonde_time = parsed["radiosonde_time"]
    else:
        radiosonde_format = infer_manual_radiosonde_format(radiosonde_file)
        radiosonde_time = target

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
    time_limit: int = 18,
    priority_time_limit: int = 3,
) -> str:
    """
    Find the best radiosonde file in `folder` using filename timestamps.

    Only files matching these filename formats are considered:
      - Wyoming: <yyyymmdd>_<hhmm>_wyoming_<wmo_id>.dat
      - Wyoming BUFR text: <yyyymmdd>_<hhmm>_wyoming_bufr_<wmo_id>.txt
      - ECMWF:   <yyyymmdd>_<hhmm>_ecmwf_<site>.nc
      - SCC:     rs_<yyyymmdd><site><hh>.nc
      - ASCII:   <yyyymmdd>_<hhmm>*.txt

    Selection rule:
      1. Prefer the closest Wyoming file within +/- `priority_time_limit` hours.
      2. Otherwise, select the closest supported file within +/- `time_limit` hours.
      3. Ignore unrelated files and folders.

    Parameters
    ----------
    target:
        Target time, usually numpy.datetime64.
    folder:
        Folder containing radiosonde files.
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
          - "radiosonde_format": "wyoming", "ecmwf", "scc", or "ascii"
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

        parsed, status = parse_radiosonde_filename(filename)

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
