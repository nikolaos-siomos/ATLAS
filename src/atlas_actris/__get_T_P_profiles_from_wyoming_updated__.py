#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wyoming radiosonde downloader using the newer University of Wyoming upper-air
endpoint shown at https://weather.uwyo.edu/upperair/sounding.shtml.

The public function keeps the same API expected by find_radiosonde.py:
    _download_wyoming(wmo_id, date, time_utc, save_dir)

It returns an object with:
    ok, path, message

The saved file is a simple comma-separated ASCII file with the column order
expected by read_radiosonde_wyoming():
    PRES,HGHT,TEMP,DWPT,MIXR,RELH
so the existing reader can continue using usecols=[1, 0, 2, 5].
"""

import os
import re
from dataclasses import dataclass
from typing import Optional
from datetime import datetime, timedelta
from io import StringIO

import numpy as np
import pandas as pd
import requests


WYOMING_NEW_ENDPOINT = "https://weather.uwyo.edu/cgi-bin/bufrraob.py"


@dataclass
class DownloadStatus:
    ok: bool
    path: Optional[str] = None
    message: str = ""
    url: Optional[str] = None


def _parse_request_datetime(date: str, time_utc: str) -> datetime:
    """Parse the date/time format currently passed by find_radiosonde."""

    date = str(date).strip()
    time_utc = str(time_utc).strip()

    for fmt in (
        "%d.%m.%Y %H:%M:%S",
        "%d.%m.%Y %H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
    ):
        try:
            return datetime.strptime(f"{date} {time_utc}", fmt)
        except ValueError:
            pass

    raise ValueError(
        f"Could not parse Wyoming radiosonde date/time: date={date}, "
        f"time_utc={time_utc}"
    )


def _nearest_wyoming_hour(dt: datetime) -> datetime:
    """Snap to the nearest nominal Wyoming sounding hour."""

    nominal_hours = (0, 3, 6, 9, 12, 15, 18, 21)

    candidates = []
    for day_offset in (-1, 0, 1):
        base = (dt + timedelta(days=day_offset)).replace(
            hour=0, minute=0, second=0, microsecond=0
        )
        for hour in nominal_hours:
            candidates.append(base + timedelta(hours=hour))

    return min(candidates, key=lambda x: abs(x - dt))


def _build_new_wyoming_url(wmo_id: str, dt: datetime, source: str = "bufr") -> str:
    """Build the newer Wyoming URL."""

    params = {
        "src": source,
        "datetime": dt.strftime("%Y-%m-%d %H:%M:%S"),
        "id": str(wmo_id),
        "type": "TEXT:CSV",
    }

    req = requests.Request("GET", WYOMING_NEW_ENDPOINT, params=params).prepare()
    return req.url


def _strip_html_if_needed(text: str) -> str:
    """Extract preformatted text if the server wraps the response in HTML."""

    if "<pre" not in text.lower():
        return text

    match = re.search(r"<pre[^>]*>(.*?)</pre>", text, flags=re.I | re.S)
    if not match:
        return text

    txt = match.group(1)
    txt = re.sub(r"<[^>]+>", "", txt)
    txt = txt.replace("&nbsp;", " ")
    txt = txt.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
    return txt


def _read_wyoming_csv(text: str) -> pd.DataFrame:
    """Read a Wyoming CSV response and return a dataframe."""

    text = _strip_html_if_needed(text)
    lines = [line.strip() for line in text.splitlines() if line.strip()]

    # Find the real CSV header. The response may contain metadata lines first.
    header_idx = None
    for i, line in enumerate(lines):
        uline = line.upper()
        if "," in line and ("PRES" in uline or "PRESS" in uline) and (
            "HGHT" in uline or "HEIGHT" in uline
        ):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("Could not find a CSV header containing pressure and height")

    csv_text = "\n".join(lines[header_idx:])
    df = pd.read_csv(StringIO(csv_text))

    # Normalize column names.
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _standardize_for_existing_reader(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert Wyoming output to the column order expected by read_radiosonde_wyoming:
        col 0 -> pressure hPa
        col 1 -> height m asl
        col 2 -> temperature C
        col 5 -> relative humidity percent
    """

    aliases = {
        "PRES": ("PRES", "PRESSURE", "P", "PRESS"),
        "HGHT": ("HGHT", "HEIGHT", "HGT", "ALT", "ALTITUDE"),
        "TEMP": ("TEMP", "TEMPERATURE", "T"),
        "DWPT": ("DWPT", "DEWPOINT", "DEW_POINT", "DEWPT", "TD"),
        "MIXR": ("MIXR", "MIXINGRATIO", "MIXING_RATIO"),
        "RELH": ("RELH", "RH", "RELATIVEHUMIDITY", "RELATIVE_HUMIDITY"),
    }

    lookup = {str(c).strip().upper().replace(" ", ""): c for c in df.columns}

    def find_col(name: str, required: bool = True):
        for alias in aliases[name]:
            key = alias.upper().replace(" ", "")
            if key in lookup:
                return lookup[key]
        if required:
            raise ValueError(
                f"Could not find required Wyoming column {name}. "
                f"Available columns: {list(df.columns)}"
            )
        return None

    pres = pd.to_numeric(df[find_col("PRES")], errors="coerce")
    hght = pd.to_numeric(df[find_col("HGHT")], errors="coerce")
    temp = pd.to_numeric(df[find_col("TEMP")], errors="coerce")

    dwpt_col = find_col("DWPT", required=False)
    mixr_col = find_col("MIXR", required=False)
    relh_col = find_col("RELH", required=False)

    dwpt = pd.to_numeric(df[dwpt_col], errors="coerce") if dwpt_col else np.nan
    mixr = pd.to_numeric(df[mixr_col], errors="coerce") if mixr_col else np.nan
    relh = pd.to_numeric(df[relh_col], errors="coerce") if relh_col else np.nan

    out = pd.DataFrame(
        {
            "PRES": pres,
            "HGHT": hght,
            "TEMP": temp,
            "DWPT": dwpt,
            "MIXR": mixr,
            "RELH": relh,
        }
    )

    # Keep only rows with the core columns needed by your reader.
    out = out.dropna(subset=["PRES", "HGHT", "TEMP"], how="any")

    if out.empty:
        raise ValueError("Wyoming response did not contain usable numeric rows")

    return out


def _download_one_wyoming_source(wmo_id: str, dt: datetime, save_dir: str, source: str) -> DownloadStatus:
    """Download one Wyoming source, usually 'bufr' first and then 'temp' as fallback."""

    url = _build_new_wyoming_url(wmo_id=wmo_id, dt=dt, source=source)

    response = requests.get(
        url,
        timeout=40,
        headers={"User-Agent": "atlas-actris-radiosonde-downloader/1.0"},
    )

    if not response.ok:
        return DownloadStatus(
            ok=False,
            message=f"HTTP {response.status_code} from Wyoming",
            url=url,
        )

    text = response.text
    if "no data" in text.lower() or "can't get" in text.lower():
        return DownloadStatus(
            ok=False,
            message="No sounding data returned by Wyoming",
            url=url,
        )

    try:
        df_raw = _read_wyoming_csv(text)
        df_out = _standardize_for_existing_reader(df_raw)
    except Exception as exc:
        return DownloadStatus(
            ok=False,
            message=f"Could not parse Wyoming response: {exc}",
            url=url,
        )

    os.makedirs(save_dir, exist_ok=True)

    # Filename keeps the yyyymmdd_hhmm prefix required by your current checks.
    fname = f"{dt:%Y%m%d_%H%M}_wyoming_{source}_{wmo_id}.txt"
    path = os.path.join(save_dir, fname)

    df_out.to_csv(path, index=False, na_rep="nan")

    return DownloadStatus(
        ok=True,
        path=path,
        message=f"Downloaded Wyoming {source.upper()} sounding",
        url=url,
    )


def _download_wyoming(wmo_id, date, time_utc, save_dir):
    """
    Download a Wyoming radiosonde from the newer endpoint.

    The input date/time can be the measurement midpoint. The downloader snaps it
    to the nearest nominal Wyoming sounding hour before requesting the file.
    BUFR is attempted first because it is the newer endpoint/data source; TEMP is
    attempted as a fallback through the same new endpoint.
    """

    request_dt = _parse_request_datetime(date=date, time_utc=time_utc)
    sounding_dt = _nearest_wyoming_hour(request_dt)

    failures = []

    for source in ("bufr", "temp"):
        status = _download_one_wyoming_source(
            wmo_id=str(wmo_id),
            dt=sounding_dt,
            save_dir=save_dir,
            source=source,
        )

        if status.ok:
            return status

        failures.append(f"{source.upper()}: {status.message}; url={status.url}")

    return DownloadStatus(
        ok=False,
        path=None,
        message=(
            "Wyoming download failed for all attempted sources. "
            f"Requested time={request_dt}, nearest sounding time={sounding_dt}. "
            + " | ".join(failures)
        ),
    )
