#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 16:16:16 2026

@author: nikos
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

RadiosondeFormat = Literal["ecmwf", "scc", "wyoming"]


class RadiosondeFormatError(Exception):
    """Base exception for radiosonde format detection."""


class UnknownRadiosondeFormatError(RadiosondeFormatError):
    """Raised when a file is not one of the supported radiosonde formats."""


class NetcdfInspectionError(RadiosondeFormatError):
    """Raised when a NetCDF-looking file cannot be inspected."""


@dataclass(frozen=True)
class NetcdfSignature:
    """Small content signature extracted from a NetCDF file."""

    root_attrs: dict[str, str]
    variables: set[str]


def infer_radiosonde_format(path: str | Path) -> RadiosondeFormat:
    """
    Infer the radiosonde format from file contents only.

    Returns one of:
        "ecmwf", "scc", "wyoming"

    This function deliberately does not use the filename.
    """
    file_path = Path(path)

    if not file_path.is_file():
        raise FileNotFoundError(f"File does not exist: {file_path}")

    if looks_like_wyoming_ascii(file_path):
        return "wyoming"

    if is_netcdf_file(file_path):
        signature = read_netcdf_signature(file_path)
        return classify_netcdf_signature(signature, path=file_path)

    raise UnknownRadiosondeFormatError(
        f"Could not infer radiosonde format for {file_path}. "
        "Expected one of: ecmwf, scc, wyoming."
    )


def infer_radiosonde_format_or_none(path: str | Path) -> Optional[RadiosondeFormat]:
    """Return the inferred format, or None if the file is not recognized."""
    try:
        return infer_radiosonde_format(path)
    except (OSError, RadiosondeFormatError):
        return None


def looks_like_wyoming_ascii(path: str | Path) -> bool:
    """
    Detect University of Wyoming ASCII/CSV output without using the filename.
    """
    file_path = Path(path)

    try:
        head = file_path.read_bytes()[:16384]
    except OSError:
        return False

    if _is_netcdf_magic(head):
        return False

    try:
        text = head.decode("utf-8-sig")
    except UnicodeDecodeError:
        return False

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return False

    header = lines[0].lower()

    required_header_tokens = ("pressure", "height", "temperature")
    if not all(token in header for token in required_header_tokens):
        return False

    if header.count(",") < 4:
        return False

    numeric_rows = 0
    for line in lines[1:]:
        if line.count(",") >= 4 and _line_has_number(line):
            numeric_rows += 1

        if numeric_rows >= 3:
            return True

    return False


def is_netcdf_file(path: str | Path) -> bool:
    """
    Return True for NetCDF classic or NetCDF4/HDF5 files by magic bytes.
    """
    try:
        head = Path(path).read_bytes()[:8]
    except OSError:
        return False

    return _is_netcdf_magic(head)


def read_netcdf_signature(path: str | Path) -> NetcdfSignature:
    """
    Extract root attributes and variable names from a NetCDF file.

    Tries h5py first, then netCDF4, then scipy.io.netcdf_file.
    """
    file_path = Path(path)

    try:
        return _read_signature_with_h5py(file_path)
    except ImportError:
        pass
    except OSError:
        pass

    try:
        return _read_signature_with_netcdf4(file_path)
    except ImportError:
        pass
    except Exception:
        pass

    try:
        return _read_signature_with_scipy(file_path)
    except ImportError as exc:
        raise NetcdfInspectionError(
            "NetCDF inspection requires at least one of h5py, netCDF4, or scipy."
        ) from exc
    except Exception as exc:
        raise NetcdfInspectionError(
            f"Could not inspect NetCDF file {file_path}: {exc}"
        ) from exc


def classify_netcdf_signature(
    signature: NetcdfSignature,
    *,
    path: str | Path | None = None,
) -> Literal["ecmwf", "scc"]:
    """Classify a NetCDF signature as either ECMWF or SCC."""
    scc_score = _scc_score(signature)
    ecmwf_score = _ecmwf_score(signature)

    # The SCC sample can contain the string "ECMWF" in WMO_Station_Number,
    # so explicit SCC radiosonde metadata must take priority.
    if scc_score >= 6 and scc_score > ecmwf_score:
        return "scc"

    if ecmwf_score >= 6 and ecmwf_score >= scc_score:
        return "ecmwf"

    suffix = f" for {Path(path)}" if path is not None else ""
    raise UnknownRadiosondeFormatError(
        f"NetCDF file{suffix} is valid NetCDF, but it does not match the known "
        "ECMWF or SCC radiosonde signatures."
    )


def _is_netcdf_magic(head: bytes) -> bool:
    return head.startswith(b"CDF") or head.startswith(b"\x89HDF\r\n\x1a\n")


def _line_has_number(line: str) -> bool:
    for token in line.replace(",", " ").split():
        try:
            float(token)
            return True
        except ValueError:
            continue

    return False


def _read_signature_with_h5py(path: Path) -> NetcdfSignature:
    import h5py  # type: ignore

    attrs: dict[str, str] = {}
    variables: set[str] = set()

    with h5py.File(path, "r") as ds:
        attrs = {
            str(key): _attr_to_string(value)
            for key, value in ds.attrs.items()
        }

        def visitor(name: str, obj: object) -> None:
            if isinstance(obj, h5py.Dataset):
                variables.add(name.split("/")[-1])

        ds.visititems(visitor)

    return NetcdfSignature(root_attrs=attrs, variables=variables)


def _read_signature_with_netcdf4(path: Path) -> NetcdfSignature:
    import netCDF4  # type: ignore

    with netCDF4.Dataset(path) as ds:
        attrs = {
            name: _attr_to_string(getattr(ds, name))
            for name in ds.ncattrs()
        }
        variables = set(ds.variables.keys())

    return NetcdfSignature(root_attrs=attrs, variables=variables)


def _read_signature_with_scipy(path: Path) -> NetcdfSignature:
    from scipy.io import netcdf_file  # type: ignore

    with netcdf_file(path, mode="r", mmap=False) as ds:
        attrs = {
            str(name): _attr_to_string(value)
            for name, value in getattr(ds, "_attributes", {}).items()
        }
        variables = set(ds.variables.keys())

    return NetcdfSignature(root_attrs=attrs, variables=variables)


def _attr_to_string(value: object) -> str:
    """Convert common NetCDF/HDF5 attribute values to plain text."""
    if isinstance(value, (bytes, bytearray)):
        return bytes(value).decode("utf-8", errors="replace")

    if hasattr(value, "tolist"):
        try:
            value = value.tolist()
        except Exception:
            pass

    if isinstance(value, (list, tuple)):
        return " ".join(_attr_to_string(item) for item in value)

    if isinstance(value, (bytes, bytearray)):
        return bytes(value).decode("utf-8", errors="replace")

    return str(value)


def _normalised_attrs(signature: NetcdfSignature) -> dict[str, str]:
    return {
        key.lower(): value.lower()
        for key, value in signature.root_attrs.items()
    }


def _normalised_vars(signature: NetcdfSignature) -> set[str]:
    return {name.lower() for name in signature.variables}


def _scc_score(signature: NetcdfSignature) -> int:
    attrs = _normalised_attrs(signature)
    variables = _normalised_vars(signature)

    score = 0

    scc_metadata_keys = {
        "measurement_type",
        "sounding_start_date",
        "sounding_start_time_ut",
        "sounding_station_name",
        "wmo_station_number",
    }

    score += len(scc_metadata_keys.intersection(attrs.keys()))

    if attrs.get("measurement_type") in {"rs", "radiosonde"}:
        score += 3

    if {"altitude", "pressure", "temperature"}.issubset(variables):
        score += 4

    return score


def _ecmwf_score(signature: NetcdfSignature) -> int:
    attrs = _normalised_attrs(signature)
    variables = _normalised_vars(signature)

    score = 0

    source = attrs.get("source", "")
    institution = attrs.get("institution", "")
    title = attrs.get("title", "")
    cloudnet_file_type = attrs.get("cloudnet_file_type", "")

    if "ecmwf" in source or "integrated forecast system" in source or "ifs" in source:
        score += 4

    if "european centre for medium-range weather forecasting" in institution:
        score += 4

    if "ecmwf" in title or "ifs" in title:
        score += 2

    if cloudnet_file_type == "model":
        score += 3

    model_markers = {
        "forecast_time",
        "level",
        "flux_level",
        "uwind",
        "vwind",
        "specific_humidity",
        "cloud_fraction",
    }

    score += len(model_markers.intersection(variables))

    if {"time", "height", "pressure", "temperature"}.issubset(variables):
        score += 2

    return score