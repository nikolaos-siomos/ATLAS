Cloudnet ECMWF Profile Reader

Overview
--------
This library provides tools to read ECMWF model profiles (temperature, pressure, altitude)
from the Cloudnet API (https://cloudnet.fmi.fi/api/model-files) for a given station and datetime.
It can automatically download, cache, and reuse the model files locally, returning
the temperature and pressure profiles at the closest model time.

You can use it:
- Programmatically inside Python (as a library)
- From the command line (CLI) for quick extraction and CSV export
- From LabVIEW via the compatibility wrapper

Installation
------------

Option 1 – Local script
Simply keep the file (for example):
    get_T_P_profiles_from_cloudnet_2025-10-20_v11.py
in your working directory and import it as:
    import get_T_P_profiles_from_cloudnet_2025-10-20_v11 as cloudnet

Option 2 – Package (optional)
If you plan to reuse it in multiple projects, create a directory structure like:
    cloudnet_profiles/
      __init__.py
      profiles.py
and move the code to profiles.py.
Then you can install it locally:
    pip install -e .

Python Usage
------------

Example:
    from get_T_P_profiles_from_cloudnet_2025-10-20_v11 import read_ecmwf_profile

    alt_m, temp_C, press_hPa, model_time_iso, nc_used = read_ecmwf_profile(
        station="Bucharest",
        date_ddmmyyyy="15.10.2025",
        time_hhmmss="12:00:00",
        save_dir="./cache"
    )

    print("Model time:", model_time_iso)
    print("File used:", nc_used)
    print("Temperature[0:5] =", temp_C[:5])

If the file already exists, it is read directly.
If not, the function downloads it automatically from Cloudnet.

Command Line Usage
------------------

You can call the script directly without writing Python code.

Usage:
    python get_T_P_profiles_from_cloudnet_2025-10-20_v11.py <station> <date> <time> --save-dir <path> [options]

Example:
    python get_T_P_profiles_from_cloudnet_2025-10-20_v11.py Bucharest 15.10.2025 12:00:00 --save-dir ./cache --outcsv-dir ./csv

Arguments:
    station       Required   Cloudnet site name or code (e.g. Bucharest)
    date          Required   Date in "dd.mm.yyyy" format
    time          Required   Time in "hh:mm:ss" UTC format
    --save-dir    Required   Local folder for caching the .nc files
    --src         Optional   Path to an existing .nc file used if download fails
    --nc          Optional   Same as --src (kept for compatibility)
    --outcsv-dir  Optional   Save output profiles as CSV file in the provided folder (defaults to save-dir)

Example output:
    Model time: 2025-10-15T12:00:00
    NC used: C:\data\cache\202510151200_bucharest_ecmwf.nc
    Levels: 137
    Saved CSV: C:\data\profile.csv

LabVIEW Usage
-------------
For LabVIEW integration, use the wrapper function lv_get_profile(), which has the same logic but accepts plain string arguments.

    from get_T_P_profiles_from_cloudnet_2025-10-20_v11 import lv_get_profile

    alt_m, temp_C, press_hPa, model_time_iso, nc_used = lv_get_profile(
        "Bucharest", "15.10.2025", "12:00:00", "./cache"
    )

Functions
---------

read_ecmwf_profile(station, date_ddmmyyyy, time_hhmmss, save_dir, src_path=None, nc_path=None)
    Reads ECMWF profile for a Cloudnet site at the closest model time.
    Returns altitude_m, temperature_C, pressure_hPa, model_time_iso, nc_path_used

lv_get_profile(station, date, time, save_dir, src_path="", nc_path="")
    Wrapper for LabVIEW or simple usage.

_cli()
    Internal function called automatically when running:
        python get_T_P_profiles_from_cloudnet_2025-10-20_v11.py ...
    Parses command-line arguments, calls read_ecmwf_profile(), and optionally exports a CSV file.

File Naming Convention
----------------------
Each downloaded or reused file is saved as:
    YYYYMMDDHHMM_station_ecmwf.nc
Example:
    202510151200_bucharest_ecmwf.nc

Outputs
-------
Variable        Units   Description
altitude_m      m       Altitude levels
temperature_C   °C      Temperature at each level
pressure_hPa    hPa     Pressure at each level
model_time_iso  ISO     Time of the ECMWF model data
nc_path_used    path    Cached NetCDF file used

Requirements
------------
Python >= 3.8
Modules: numpy, pandas, xarray, urllib, ssl, json

Install dependencies:
    pip install numpy pandas xarray

Version
-------
v11 – 20.10.2025
Author: Livio Belegante, INOE
Contact: livio@inoe.ro

License
-------
MIT License (optional — you can change this)
