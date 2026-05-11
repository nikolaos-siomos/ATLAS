#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 13:50:57 2026

@author: nikos, livio
"""
import os
import ssl
import json
import shutil
import tempfile
import urllib.request
from datetime import datetime, date as _date, time as _time

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