"""
Read ECMWF profiles from Cloudnet.

This module provides functions to download or reuse ECMWF model files from the
Cloudnet API for a given station and datetime, then extract vertical profiles
of altitude, temperature, and pressure at the closest model time.

Inputs
------
station : str
    Station name or code, case insensitive.
date : "dd.mm.yyyy"
    Requested date string.
time : "hh:mm:ss"
    Requested time string (UTC).
save_dir : str
    Directory where the NetCDF file will be cached.
src_path : str, optional
    Alternate source .nc to copy into the cache if not present.
nc_path : str, optional
    Alternate source .nc to copy into the cache if not present.

Outputs
-------
altitude_m : list[float]
temperature_C : list[float]
pressure_hPa : list[float]
model_time_iso : str
    Model timestamp in ISO format "YYYY-MM-DDTHH:MM:SS".
nc_path_used : str
    Absolute path to the NetCDF file used.

Version
-------
11 (20.10.2025)

Author
------
Livio Belegante, INOE
Contact: livio@inoe.ro
"""

from __future__ import annotations

import os
import json
import ssl
import shutil
import tempfile
import urllib.request
import urllib.error
from datetime import datetime, date as _date, time as _time
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import xarray as xr

import argparse
import csv


__all__ = ["read_ecmwf_profile", "lv_get_profile"]

# ---------------- Utilities ----------------
def _closest_time_index(time_coord, target_dt: datetime) -> tuple[int, pd.Timestamp]:
    """Return the index of the closest time and the matched pandas Timestamp."""
    if hasattr(time_coord, "values"):
        tvals = pd.to_datetime(time_coord.values)
    else:
        tvals = pd.to_datetime(time_coord)
    diffs = np.abs(tvals - pd.Timestamp(target_dt))
    idx = int(np.argmin(diffs))
    return idx, pd.Timestamp(tvals[idx])


def _ensure_1d(arr, name: str) -> np.ndarray:
    """Squeeze and validate that an array is 1D."""
    a = np.asarray(arr).squeeze()
    if a.ndim != 1:
        raise ValueError(f"{name} is not 1-D after squeeze, shape={a.shape}")
    return a


def _expected_name(date_obj: _date, time_obj: _time, station: str) -> str:
    """Build the cache filename YYYYMMDDHHMM_station_ecmwf.nc."""
    hhmm = f"{time_obj.hour:02d}{time_obj.minute:02d}"
    return f"{date_obj.strftime('%Y%m%d')}{hhmm}_{station.lower()}_ecmwf.nc"


def _ensure_dir(path: str) -> str:
    """Create directory if missing and return absolute path."""
    os.makedirs(path, exist_ok=True)
    return os.path.abspath(path)


# ---------------- Cloudnet downloader ----------------
def _download_cloudnet_ecmwf(station: str, date_obj: _date, out_path: str) -> str:
    """
    Download the daily ECMWF model file from the Cloudnet API.

    The API endpoint is:
    GET https://cloudnet.fmi.fi/api/model-files?site=<site>&date=<YYYY-MM-DD>&model=ecmwf

    Parameters
    ----------
    station : str
        Cloudnet site, case insensitive.
    date_obj : datetime.date
        Target date.
    out_path : str
        Destination NetCDF filepath.

    Returns
    -------
    str
        The `out_path` where the file was saved.

    Raises
    ------
    RuntimeError
        If the API returns an error or the response is malformed.
    FileNotFoundError
        If the API returns no files.
    """
    site = station.strip().lower()
    date_iso = date_obj.strftime("%Y-%m-%d")
    api_url = f"https://cloudnet.fmi.fi/api/model-files?site={site}&date={date_iso}&model=ecmwf"

    ctx = ssl.create_default_context()
    try:
        with urllib.request.urlopen(api_url, context=ctx, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Cloudnet API HTTP error {e.code} for {api_url}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Cloudnet API connection error for {api_url}: {e.reason}") from e

    if not isinstance(data, list) or not data:
        raise FileNotFoundError(f"No ECMWF model files for site={site} date={date_iso}")

    chosen = None
    for item in data:
        dl = item.get("downloadUrl")
        if dl:
            chosen = dl
            break
    if not chosen:
        raise RuntimeError("downloadUrl missing in Cloudnet API response")

    tmpfd, tmpname = tempfile.mkstemp(prefix="cn_ecmwf_", suffix=".nc")
    os.close(tmpfd)
    try:
        with urllib.request.urlopen(chosen, context=ctx, timeout=300) as r, open(tmpname, "wb") as f:
            while True:
                chunk = r.read(8192)
                if not chunk:
                    break
                f.write(chunk)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        shutil.move(tmpname, out_path)
    except Exception:
        try:
            os.remove(tmpname)
        except Exception:
            pass
        raise

    return out_path


# ---------------- Storage and caching ----------------
def _ensure_nc_saved(
    date_obj: _date,
    time_obj: _time,
    station: str,
    save_dir: str,
    src_path: Optional[str] = None,
    nc_path: Optional[str] = None,
    try_download: bool = True,
) -> str:
    """
    Ensure the expected NetCDF is present in the cache, else seed or download it.

    Strategy
    1) Reuse the cached file if present.
    2) Else, copy from `nc_path` or `src_path` if provided.
    3) Else, download the daily Cloudnet file and store it under the expected name.

    Returns
    -------
    str
        Absolute path to the cached NetCDF.
    """
    save_dir = _ensure_dir(save_dir)
    expected_name = _expected_name(date_obj, time_obj, station)
    expected = os.path.join(save_dir, expected_name)

    if os.path.isfile(expected):
        return expected

    for label, candidate in (("nc_path", nc_path), ("src_path", src_path)):
        if candidate:
            src = os.path.abspath(candidate)
            if not os.path.isfile(src):
                raise FileNotFoundError(f"{label} not found: {src}")
            shutil.copy(src, expected)
            return expected

    if try_download:
        return _download_cloudnet_ecmwf(station, date_obj, expected)

    raise FileNotFoundError(
        f"Model file missing in cache and download disabled.\nExpected: {expected}"
    )


# ---------------- Profile extraction ----------------
def _extract_profile(nc_path: str, target_dt: datetime) -> Tuple[List[float], List[float], List[float], str]:
    """
    Open a Cloudnet NetCDF, select the closest model time, and extract the profile.

    Returns altitude [m], temperature [C], pressure [hPa], and the matched model time ISO.
    """
    with xr.open_dataset(nc_path) as ds:
        needed = ["time", "temperature", "pressure", "height"]
        present = set(list(ds.data_vars) + list(ds.coords))
        missing = [v for v in needed if v not in present]
        if missing:
            raise KeyError(f"Missing variables in {nc_path}: {missing}. Present: {sorted(present)}")

        tidx, used_dt = _closest_time_index(ds["time"], target_dt)

        z_m = _ensure_1d(ds["height"].isel(time=tidx).values, "height")
        T_K = _ensure_1d(ds["temperature"].isel(time=tidx).values, "temperature")
        P_Pa = _ensure_1d(ds["pressure"].isel(time=tidx).values, "pressure")

        if not (z_m.shape == T_K.shape == P_Pa.shape):
            raise ValueError(f"Shape mismatch. z={z_m.shape}, T={T_K.shape}, P={P_Pa.shape}")

    T_C = (T_K - 273.15).astype(float).tolist()
    P_hPa = (P_Pa / 100.0).astype(float).tolist()
    Z_m = z_m.astype(float).tolist()

    return Z_m, T_C, P_hPa, used_dt.isoformat(timespec="seconds")


# ---------------- Public library API ----------------
def read_ecmwf_profile(
    station: str,
    date_ddmmyyyy: str,
    time_hhmmss: str,
    save_dir: str,
    src_path: Optional[str] = None,
    nc_path: Optional[str] = None,
) -> Tuple[List[float], List[float], List[float], str, str]:
    """
    Read ECMWF profile for a Cloudnet site at the closest model time.

    Parameters
    ----------
    station : str
        Station name or code.
    date_ddmmyyyy : str
        Date in format "dd.mm.yyyy".
    time_hhmmss : str
        Time in format "hh:mm:ss" UTC.
    save_dir : str
        Cache directory for the NetCDF files.
    src_path : str, optional
        Seed NetCDF copied into cache if missing.
    nc_path : str, optional
        Seed NetCDF copied into cache if missing.

    Returns
    -------
    altitude_m : list[float]
    temperature_C : list[float]
    pressure_hPa : list[float]
    model_time_iso : str
        ISO timestamp of the matched model time.
    nc_path_used : str
        Absolute path to the NetCDF used.

    Raises
    ------
    ValueError
        If inputs are malformed.
    FileNotFoundError
        If no model file is available and download is disabled.
    RuntimeError
        If Cloudnet API fails.
    """
    station = str(station).strip()
    if not station:
        raise ValueError("station is empty")

    try:
        req_date = datetime.strptime(date_ddmmyyyy, "%d.%m.%Y").date()
    except ValueError:
        raise ValueError("date must be dd.mm.yyyy")

    try:
        req_time = datetime.strptime(time_hhmmss, "%H:%M:%S").time()
    except ValueError:
        raise ValueError("time must be hh:mm:ss (24h)")

    nc_used = _ensure_nc_saved(
        req_date, req_time, station, save_dir,
        src_path=src_path,
        nc_path=nc_path,
        try_download=True,
    )

    target_dt = datetime.combine(req_date, req_time)
    alt_m, temp_C, press_hPa, model_time_iso = _extract_profile(nc_used, target_dt)

    return alt_m, temp_C, press_hPa, model_time_iso, os.path.abspath(nc_used)


# ---------------- Backward compatibility for LabVIEW ----------------
def lv_get_profile(
    station: str,
    date_ddmmyyyy: str,
    time_hhmmss: str,
    save_dir: str,
    src_path: str = "",
    nc_path: str = "",
):
    """
    LabVIEW entry point wrapper around `read_ecmwf_profile`.

    Same signature used previously in LabVIEW integration.
    """
    return read_ecmwf_profile(
        station=station,
        date_ddmmyyyy=date_ddmmyyyy,
        time_hhmmss=time_hhmmss,
        save_dir=save_dir,
        src_path=src_path or None,
        nc_path=nc_path or None,
    )

def export_to_csv(outcsv, model_time_iso, station, alt_m, temp_C, press_hPa):
    # Exported csv filename
    dt = datetime.strptime(model_time_iso, "%Y-%m-%dT%H:%M:%S")
    filename_part = dt.strftime("%Y%m%d_%H%M")
    default_csv = f"{filename_part}_{station}.txt"

    csv_path = os.path.abspath(os.path.join(outcsv, default_csv))

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["altitude_m", "temperature_C", "pressure_hPa"])
        for z, T, P in zip(alt_m, temp_C, press_hPa):
            w.writerow([z, T, P])
    print("Saved CSV:",csv_path)

# ---------------- CLI for quick tests ----------------
def _cli():

    ap = argparse.ArgumentParser(description="Extract ECMWF profile for a Cloudnet site. Downloads if missing.")
    ap.add_argument("station", help='e.g. "Bucharest"')
    ap.add_argument("date", help='date "dd.mm.yyyy"')
    ap.add_argument("time", help='time "hh:mm:ss" UTC')
    ap.add_argument("--save-dir", required = True, help="cache directory for .nc files")
    ap.add_argument("--src", default=None, help="optional seed .nc to copy if cache empty")
    ap.add_argument("--nc", default=None, help="optional seed .nc to copy if cache empty")
    ap.add_argument("--outcsv-dir", default=None, help="CSV export folder path (defaults to the save-dir folder)")
    args = ap.parse_args()

    # Create directories inside the script parent folder if save-dir or outcsv are not provided
        
    if args.outcsv_dir is None:
        default_dir = args.save_dir        
        args.outcsv_dir = os.path.join(default_dir)

    os.makedirs(args.outcsv_dir, exist_ok=True)

    alt_m, temp_C, press_hPa, model_time_iso, nc_used = read_ecmwf_profile(
        args.station, args.date, args.time, args.save_dir, src_path=args.src, nc_path=args.nc
    )
    
    print()

    print("Model time:", model_time_iso)
    print(f"Levels: {len(alt_m)}")
    print()
    print("Saved netcdf:", nc_used)
    
    export_to_csv(
        outcsv = args.outcsv_dir, 
        model_time_iso = model_time_iso, 
        station = args.station,
        alt_m = alt_m, 
        temp_C = temp_C,
        press_hPa = press_hPa
        )
    print()



if __name__ == "__main__":
    _cli()
