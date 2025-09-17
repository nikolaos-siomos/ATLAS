#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 14:42:46 2025

@author: nikos
"""

def is_licel_header(path, delimiter=b"\n", max_lines=100):
    """
    Check whether the header of a file (up to the first delimiter or max_lines)
    contains only ASCII characters (0–127).

    Parameters
    ----------
    path : str
        Path to the file.
    delimiter : bytes, optional
        Delimiter that marks the end of the ASCII header (default: empty line b"\\n").
    max_lines : int, optional
        Maximum number of lines to scan before giving up.

    Returns
    -------
    bool
        True if a delimiter is found within max_lines and all preceding
        lines are ASCII-only; False otherwise.
    """
    with open(path, "rb") as f:
        for i in range(max_lines):
            line = f.readline()
            if not line:  # EOF reached before finding delimiter
                return False
            if max(line) >= 128:  # found non-ASCII byte
                return False
            if line.strip() == delimiter.strip():  # delimiter found (e.g. empty line)
                return True
    return False  # delimiter not found within max_lines


def detect_netcdf(path):
    """
    Return one of:
      - 'NETCDF3_CLASSIC'
      - 'NETCDF3_64BIT_OFFSET'
      - 'NETCDF3_64BIT_DATA'   (aka CDF-5)
      - 'NETCDF4'              (HDF5-based, netCDF structure)
      - None                   (not a NetCDF file)
    """
    # Check magic bytes
    with open(path, "rb") as f:
        head8 = f.read(8)

    # NetCDF-3 signatures: CDF\001, CDF\002, CDF\005
    if head8.startswith(b"CDF"):
        ver = head8[3:4]
        return {
            b"\x01": "NETCDF3_CLASSIC",
            b"\x02": "NETCDF3_64BIT_OFFSET",
            b"\x05": "NETCDF3_64BIT_DATA",
        }.get(ver, "NETCDF3_CLASSIC")

    # HDF5 signature → could be NetCDF-4 or generic HDF5
    if head8 == b"\x89HDF\r\n\x1a\n":
        try:
            # Only import if needed
            from netCDF4 import Dataset
            with Dataset(path, "r") as _:
                return "NETCDF4"
        except Exception:
            # It's HDF5 but not a netCDF-4 file
            return None

    return None