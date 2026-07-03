#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wyoming radiosonde downloader using the University of Wyoming upper-air endpoints.

The public function keeps the same API expected by find_radiosonde.py:
    _download_wyoming(wmo_id, date, time_utc, save_dir)

It returns an object with:
    ok, path, message, url

Workflow:
    1) Parse requested measurement date/time.
    2) Build candidate nominal radiosonde times within +/-12 hours.
    3) Reuse an existing local file if already downloaded.
    4) Try the current Wyoming WSGI endpoint:
           https://weather.uwyo.edu/wsgi/sounding
       with src=FM35 and type=TEXT:CSV.
    5) If this fails, try the older fallbacks:
           https://weather.uwyo.edu/cgi-bin/bufrraob.py
           http://weather.uwyo.edu/cgi-bin/sounding
    6) Save the first valid response using the ATLAS filename pattern:
           YYYYMMDD_HHMM_wyoming_<source>_<wmo_id>.txt

Version:
    2026-07-01 WSGI FM35 endpoint + old fallbacks
    contact livio.belegante@actris.ro
"""

import argparse
import html
import json
import os
import re
import tempfile
import urllib.parse
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple

import requests


WYOMING_WSGI_ENDPOINT = "https://weather.uwyo.edu/wsgi/sounding"
WYOMING_NEW_ENDPOINT = "https://weather.uwyo.edu/cgi-bin/bufrraob.py"
WYOMING_OLD_ENDPOINT = "http://weather.uwyo.edu/cgi-bin/sounding"
WYOMING_FAILURE_CACHE = ".wyoming_download_failures.json"
WYOMING_USER_AGENT = "atlas-actris-radiosonde-downloader/1.0"
WYOMING_TIMEOUT = (3.05, 12.0)

WYOMING_DEBUG = False

OLD_START_MARKER = (
    "-----------------------------------------------------------------------------\n"
    "   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV\n"
    "    hPa     m      C      C      %    g/kg    deg   knot     K      K      K \n"
    "-----------------------------------------------------------------------------"
)

OLD_END_MARKER = "</PRE><H3>Station information and sounding indices</H3><PRE>"


@dataclass
class DownloadStatus:
    ok: bool
    path: Optional[str] = None
    message: str = ""
    url: Optional[str] = None


# ----------- _debug ----------
def _debug(message: str) -> None:
    """Print a debug message only when WYOMING_DEBUG is enabled."""
    if WYOMING_DEBUG:
        print(f"[DEBUG] {message}")
# -----------------------------


# ----------- _parse_request_datetime ----------
def _parse_request_datetime(date: str, time_utc: str) -> datetime:
    """Parse the date/time format currently passed by find_radiosonde."""
    _debug(f"Parsing requested date/time: date={date}, time_utc={time_utc}")

    date = str(date).strip()
    time_utc = str(time_utc).strip()

    for fmt in (
        "%d.%m.%Y %H:%M:%S",
        "%d.%m.%Y %H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
    ):
        try:
            parsed = datetime.strptime(f"{date} {time_utc}", fmt)
            _debug(f"Parsed requested datetime as {parsed:%Y-%m-%d %H:%M:%S} UTC using format {fmt}")
            return parsed
        except ValueError:
            pass

    raise ValueError(
        f"Could not parse Wyoming radiosonde date/time: date={date}, "
        f"time_utc={time_utc}"
    )
# -----------------------------


# ----------- _wyoming_candidate_times ----------
def _wyoming_candidate_times(dt: datetime, search_hours: int = 12) -> List[datetime]:
    """Return nominal Wyoming sounding times within +/- search_hours."""
    _debug(f"Building candidate sounding times within +/-{search_hours} h of {dt:%Y-%m-%d %H:%M:%S}")

    nominal_hours = (0, 3, 6, 9, 12, 15, 18, 21)
    start_dt = dt - timedelta(hours=search_hours)
    end_dt = dt + timedelta(hours=search_hours)

    candidates = []
    for day_offset in (-1, 0, 1):
        base = (dt + timedelta(days=day_offset)).replace(
            hour=0, minute=0, second=0, microsecond=0
        )
        for hour in nominal_hours:
            candidate = base + timedelta(hours=hour)
            if start_dt <= candidate <= end_dt:
                candidates.append(candidate)

    candidates = sorted(set(candidates), key=lambda x: (abs(x - dt), x))

    if not candidates:
        raise ValueError(
            f"No Wyoming nominal sounding times found within +/-{search_hours} h "
            f"of {dt}"
        )

    _debug(
        "Candidate order: "
        + ", ".join(candidate.strftime("%Y-%m-%d %H:%M") for candidate in candidates)
    )

    return candidates
# -----------------------------


# ----------- _nearest_wyoming_hour ----------
def _nearest_wyoming_hour(dt: datetime) -> datetime:
    """Snap to the nearest nominal Wyoming sounding hour."""
    candidates = _wyoming_candidate_times(dt=dt, search_hours=12)
    return min(candidates, key=lambda x: abs(x - dt))
# -----------------------------


# ----------- _wyoming_filename ----------
def _wyoming_filename(wmo_id: str, dt: datetime, source: str) -> str:
    """Return the ATLAS-expected local filename for a Wyoming download."""
    return f"{dt:%Y%m%d_%H%M}_wyoming_{source}_{wmo_id}.txt"
# -----------------------------


# ----------- _existing_wyoming_file ----------
def _existing_wyoming_file(wmo_id: str, dt: datetime, save_dir: str) -> Optional[str]:
    """Return an already-downloaded Wyoming file for this WMO/time, if present."""
    if not save_dir or not os.path.isdir(save_dir):
        _debug(f"No existing folder found for cache check: {save_dir}")
        return None

    prefix = f"{dt:%Y%m%d_%H%M}_wyoming_"
    suffix = f"_{wmo_id}.txt"

    _debug(f"Checking existing files with prefix={prefix}, suffix={suffix}")

    for fname in sorted(os.listdir(save_dir)):
        if fname.startswith(prefix) and fname.endswith(suffix):
            path = os.path.join(save_dir, fname)
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                _debug(f"Found existing Wyoming file: {path}")
                return path

    _debug("No existing Wyoming file found for this candidate time")
    return None
# -----------------------------


# ----------- _build_new_wyoming_url ----------
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
# -----------------------------


# ----------- _build_wsgi_wyoming_url ----------
def _build_wsgi_wyoming_url(wmo_id: str, dt: datetime, source: str = "FM35") -> str:
    """Build the current Wyoming WSGI CSV URL."""
    params = {
        "datetime": dt.strftime("%Y-%m-%d %H:%M:%S"),
        "id": str(wmo_id),
        "src": source,
        "type": "TEXT:CSV",
    }

    req = requests.Request("GET", WYOMING_WSGI_ENDPOINT, params=params).prepare()
    return req.url
# -----------------------------


# ----------- _build_old_wyoming_url ----------
def _build_old_wyoming_url(wmo_id: str, dt: datetime) -> str:
    """Build the older Wyoming TEXT:LIST URL."""
    params = {
        "region": "europe",
        "TYPE": "TEXT:LIST",
        "YEAR": f"{dt.year:04d}",
        "MONTH": f"{dt.month:02d}",
        "FROM": f"{dt.day:02d}{dt.hour:02d}",
        "TO": f"{dt.day:02d}{dt.hour:02d}",
        "STNM": str(wmo_id),
    }

    return WYOMING_OLD_ENDPOINT + "?" + urllib.parse.urlencode(params)
# -----------------------------


# ----------- _failure_cache_path ----------
def _failure_cache_path(save_dir: str) -> str:
    """Return the persistent failure cache path."""
    return os.path.join(save_dir, WYOMING_FAILURE_CACHE)
# -----------------------------


# ----------- _load_failure_cache ----------
def _load_failure_cache(save_dir: str) -> Dict[str, Dict[str, str]]:
    """Load cached Wyoming failures. Invalid caches are ignored."""
    path = _failure_cache_path(save_dir)
    _debug(f"Loading failure cache: {path}")

    if not os.path.isfile(path):
        _debug("Failure cache does not exist yet")
        return {}

    try:
        with open(path, "r", encoding="utf-8") as fobj:
            data = json.load(fobj)
    except Exception as exc:
        _debug(f"Could not read failure cache, ignoring it: {exc}")
        return {}

    if not isinstance(data, dict):
        _debug("Failure cache is not a dictionary, ignoring it")
        return {}

    _debug(f"Loaded {len(data)} cached failures")
    return data
# -----------------------------


# ----------- _save_failure_cache ----------
def _save_failure_cache(save_dir: str, cache: Dict[str, Dict[str, str]]) -> None:
    """Atomically save the Wyoming failure cache."""
    _debug(f"Saving failure cache with {len(cache)} entries")

    os.makedirs(save_dir, exist_ok=True)
    path = _failure_cache_path(save_dir)
    fd, tmp_path = tempfile.mkstemp(
        prefix=".wyoming_download_failures.", suffix=".tmp", dir=save_dir
    )

    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fobj:
            json.dump(cache, fobj, indent=2, sort_keys=True)
        os.replace(tmp_path, path)
    except Exception:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise
# -----------------------------


# ----------- _failure_key ----------
def _failure_key(wmo_id: str, dt: datetime, source: str) -> str:
    """Return a stable cache key for one station/time/source request."""
    return f"{wmo_id}|{dt:%Y%m%d%H%M}|{source}"
# -----------------------------


# ----------- _cache_failure ----------
def _cache_failure(
    cache: Dict[str, Dict[str, str]],
    wmo_id: str,
    dt: datetime,
    source: str,
    message: str,
    url: str,
) -> None:
    """Store one failed station/time/source attempt in the cache."""
    _debug(f"Caching failure for {dt:%Y-%m-%d %H:%M} source={source}: {message}")

    cache[_failure_key(wmo_id=wmo_id, dt=dt, source=source)] = {
        "wmo_id": str(wmo_id),
        "datetime_utc": dt.strftime("%Y-%m-%d %H:%M:%S"),
        "source": str(source),
        "message": str(message),
        "url": str(url),
        "cached_at_utc": datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),
    }
# -----------------------------


# ----------- _cached_failure ----------
def _cached_failure(
    cache: Dict[str, Dict[str, str]], wmo_id: str, dt: datetime, source: str
) -> Optional[Dict[str, str]]:
    """Return cached failure metadata for one request, if present."""
    return cache.get(_failure_key(wmo_id=wmo_id, dt=dt, source=source))
# -----------------------------


# ----------- _looks_like_failed_wyoming_response ----------
def _looks_like_failed_wyoming_response(content: bytes, allow_html: bool = False) -> Tuple[bool, str]:
    """Detect obvious Wyoming failure pages without parsing the CSV data."""
    if not content:
        return True, "Empty response from Wyoming"

    sample = content[:4096].decode("utf-8", errors="ignore").lower()
    compact = re.sub(r"\s+", " ", sample)

    failure_markers = (
        "no data",
        "can't get",
        "cannot get",
        "not found",
        "invalid station",
        "error",
        "service unavailable",
        "temporarily unavailable",
    )

    for marker in failure_markers:
        if marker in compact:
            return True, f"Wyoming response indicates failure: {marker}"

    if not allow_html:
        if b"," not in content[:8192] and b"<html" in content[:8192].lower():
            return True, "Wyoming returned an HTML page instead of CSV data"

    return False, ""
# -----------------------------


# ----------- _decode_response_text ----------
def _decode_response_text(content: bytes) -> str:
    """Decode Wyoming response bytes to text."""
    try:
        return content.decode("utf-8")
    except UnicodeDecodeError:
        return content.decode("latin-1", errors="replace")
# -----------------------------


# ----------- _extract_old_wyoming_block ----------
def _extract_old_wyoming_block(raw_html: str) -> List[str]:
    """Extract sounding data lines from the older Wyoming HTML/TEXT:LIST response."""
    text = raw_html.replace("\r\n", "\n").replace("\r", "\n")

    start_idx = text.find(OLD_START_MARKER)
    if start_idx == -1:
        raise RuntimeError("Could not find the old Wyoming header block")

    data_start = start_idx + len(OLD_START_MARKER)
    end_idx = text.find(OLD_END_MARKER, data_start)

    if end_idx == -1:
        raise RuntimeError("Could not find the old Wyoming station information section")

    sub = text[data_start:end_idx]
    lines = sub.split("\n")

    while lines and not lines[0].strip():
        lines.pop(0)

    while lines and not lines[-1].strip():
        lines.pop()

    return lines
# -----------------------------


# ----------- _normalize_old_wyoming_lines ----------
def _normalize_old_wyoming_lines(lines: List[str], sep: str = ",") -> str:
    """Convert old Wyoming whitespace-separated lines to separator-separated lines."""
    out_lines = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()
        out_lines.append(sep.join(parts))

    return "\n".join(out_lines) + "\n"
# -----------------------------


# ----------- _format_old_wyoming_content ----------
def _format_old_wyoming_content(lines: List[str], sep: str = ",") -> bytes:
    """Create ATLAS-readable text content from old Wyoming TEXT:LIST data lines."""
    if sep == ",":
        name_line = "PRES,HGHT,TEMP,DWPT,RELH,MIXR,DRCT,SKNT,THTA,THTE,THTV"
        unit_line = "hPa,m,C,C,%,g/kg,deg,knot,K,K,K"
    else:
        cols = ["PRES", "HGHT", "TEMP", "DWPT", "RELH", "MIXR", "DRCT", "SKNT", "THTA", "THTE", "THTV"]
        units = ["hPa", "m", "C", "C", "%", "g/kg", "deg", "knot", "K", "K", "K"]
        name_line = sep.join(cols)
        unit_line = sep.join(units)

    header_block = (
        "-----------------------------------------------------------------------------\n"
        f"{name_line}\n"
        f"{unit_line}\n"
        "-----------------------------------------------------------------------------\n"
    )

    content = header_block + _normalize_old_wyoming_lines(lines=lines, sep=sep)
    return content.encode("utf-8")
# -----------------------------


# ----------- _write_bytes_atomic ----------
def _write_bytes_atomic(path: str, content: bytes) -> None:
    """Write downloaded content atomically."""
    _debug(f"Writing file atomically: {path}")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(path)}.", suffix=".tmp", dir=os.path.dirname(path)
    )

    try:
        with os.fdopen(fd, "wb") as fobj:
            fobj.write(content)
        os.replace(tmp_path, path)
    except Exception:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise
# -----------------------------


# ----------- _download_one_wyoming_wsgi_source ----------
def _download_one_wyoming_wsgi_source(
    wmo_id: str,
    dt: datetime,
    save_dir: str,
    session: requests.Session,
    source: str = "FM35",
) -> DownloadStatus:
    """Download one Wyoming sounding from the current WSGI CSV endpoint."""
    url = _build_wsgi_wyoming_url(wmo_id=wmo_id, dt=dt, source=source)
    _debug(f"Trying WSGI Wyoming endpoint: source={source}, url={url}")

    try:
        response = session.get(url, timeout=WYOMING_TIMEOUT)
    except requests.RequestException as exc:
        _debug(f"WSGI endpoint request exception: {exc}")
        return DownloadStatus(
            ok=False,
            message=f"Request to WSGI Wyoming endpoint failed: {exc}",
            url=url,
        )

    _debug(f"WSGI endpoint HTTP status: {response.status_code}")

    if not response.ok:
        return DownloadStatus(
            ok=False,
            message=f"HTTP {response.status_code} from WSGI Wyoming endpoint",
            url=url,
        )

    failed, failure_message = _looks_like_failed_wyoming_response(response.content, allow_html=False)
    if failed:
        _debug(f"WSGI endpoint response rejected: {failure_message}")
        return DownloadStatus(ok=False, message=failure_message, url=url)

    fname = _wyoming_filename(wmo_id=wmo_id, dt=dt, source=source.lower())
    path = os.path.join(save_dir, fname)
    _write_bytes_atomic(path=path, content=response.content)

    _debug(f"WSGI endpoint success. Saved to {path}")

    return DownloadStatus(
        ok=True,
        path=path,
        message=(
            f"Downloaded raw Wyoming {source} CSV sounding "
            f"for {dt:%Y-%m-%d %H:%M} UTC using WSGI endpoint"
        ),
        url=url,
    )
# -----------------------------


# ----------- _download_one_wyoming_new_source ----------
def _download_one_wyoming_new_source(
    wmo_id: str,
    dt: datetime,
    save_dir: str,
    source: str,
    session: requests.Session,
) -> DownloadStatus:
    """Download one Wyoming source from the newer BUFR/RAOB endpoint."""
    url = _build_new_wyoming_url(wmo_id=wmo_id, dt=dt, source=source)
    _debug(f"Trying NEW Wyoming endpoint: source={source.upper()}, url={url}")

    try:
        response = session.get(url, timeout=WYOMING_TIMEOUT)
    except requests.RequestException as exc:
        _debug(f"NEW endpoint request exception: {exc}")
        return DownloadStatus(
            ok=False,
            message=f"Request to new Wyoming endpoint failed: {exc}",
            url=url,
        )

    _debug(f"NEW endpoint HTTP status: {response.status_code}")

    if not response.ok:
        return DownloadStatus(
            ok=False,
            message=f"HTTP {response.status_code} from new Wyoming endpoint",
            url=url,
        )

    failed, failure_message = _looks_like_failed_wyoming_response(response.content, allow_html=False)
    if failed:
        _debug(f"NEW endpoint response rejected: {failure_message}")
        return DownloadStatus(ok=False, message=failure_message, url=url)

    fname = _wyoming_filename(wmo_id=wmo_id, dt=dt, source=source)
    path = os.path.join(save_dir, fname)
    _write_bytes_atomic(path=path, content=response.content)

    _debug(f"NEW endpoint success. Saved to {path}")

    return DownloadStatus(
        ok=True,
        path=path,
        message=(
            f"Downloaded raw Wyoming {source.upper()} CSV sounding "
            f"for {dt:%Y-%m-%d %H:%M} UTC using new endpoint"
        ),
        url=url,
    )
# -----------------------------


# ----------- _download_one_wyoming_old_source ----------
def _download_one_wyoming_old_source(
    wmo_id: str,
    dt: datetime,
    save_dir: str,
    session: requests.Session,
) -> DownloadStatus:
    """Download one Wyoming sounding from the older TEXT:LIST endpoint."""
    source = "old"
    url = _build_old_wyoming_url(wmo_id=wmo_id, dt=dt)
    _debug(f"Trying OLD Wyoming endpoint: url={url}")

    try:
        response = session.get(url, timeout=WYOMING_TIMEOUT)
    except requests.RequestException as exc:
        _debug(f"OLD endpoint request exception: {exc}")
        return DownloadStatus(
            ok=False,
            message=f"Request to old Wyoming endpoint failed: {exc}",
            url=url,
        )

    _debug(f"OLD endpoint HTTP status: {response.status_code}")

    if not response.ok:
        return DownloadStatus(
            ok=False,
            message=f"HTTP {response.status_code} from old Wyoming endpoint",
            url=url,
        )

    failed, failure_message = _looks_like_failed_wyoming_response(response.content, allow_html=True)
    if failed:
        _debug(f"OLD endpoint response rejected before parsing: {failure_message}")
        return DownloadStatus(ok=False, message=failure_message, url=url)

    try:
        raw_html = _decode_response_text(response.content)
        raw_html = html.unescape(raw_html)
        lines = _extract_old_wyoming_block(raw_html=raw_html)
        content = _format_old_wyoming_content(lines=lines, sep=",")
    except Exception as exc:
        _debug(f"OLD endpoint parsing failed: {exc}")
        return DownloadStatus(
            ok=False,
            message=f"Could not parse old Wyoming endpoint response: {exc}",
            url=url,
        )

    fname = _wyoming_filename(wmo_id=wmo_id, dt=dt, source=source)
    path = os.path.join(save_dir, fname)
    _write_bytes_atomic(path=path, content=content)

    _debug(f"OLD endpoint success. Saved to {path}")

    return DownloadStatus(
        ok=True,
        path=path,
        message=(
            f"Downloaded and reformatted Wyoming TEXT:LIST sounding "
            f"for {dt:%Y-%m-%d %H:%M} UTC using old endpoint"
        ),
        url=url,
    )
# -----------------------------


# ----------- _download_wyoming ----------
def _download_wyoming(wmo_id, date, time_utc, save_dir):
    """
    Download the closest available Wyoming radiosonde within +/-12 hours.

    The current WSGI endpoint is tried first:
        wsgi/sounding with source=FM35 and type=TEXT:CSV

    If it fails for a candidate time, the old fallbacks are used:
        bufrraob.py with source=bufr
        bufrraob.py with source=temp
        sounding with TYPE=TEXT:LIST
    """
    _debug("Starting Wyoming download workflow")

    request_dt = _parse_request_datetime(date=date, time_utc=time_utc)
    candidate_times = _wyoming_candidate_times(dt=request_dt, search_hours=12)
    wmo_id = str(wmo_id).strip()
    os.makedirs(save_dir, exist_ok=True)

    _debug(f"WMO ID: {wmo_id}")
    _debug(f"Output directory: {os.path.abspath(save_dir)}")

    cache = _load_failure_cache(save_dir=save_dir)
    failures = []
    attempted_network = False
    cache_changed = False

    session = requests.Session()
    session.headers.update({"User-Agent": WYOMING_USER_AGENT})

    try:
        for sounding_dt in candidate_times:
            _debug(f"Processing candidate time: {sounding_dt:%Y-%m-%d %H:%M} UTC")

            existing_path = _existing_wyoming_file(
                wmo_id=wmo_id,
                dt=sounding_dt,
                save_dir=save_dir,
            )

            if existing_path is not None:
                _debug("Returning existing local file")
                return DownloadStatus(
                    ok=True,
                    path=existing_path,
                    message=(
                        "Using already-downloaded Wyoming sounding "
                        f"for {sounding_dt:%Y-%m-%d %H:%M} UTC"
                    ),
                )

            cached_wsgi = _cached_failure(
                cache=cache,
                wmo_id=wmo_id,
                dt=sounding_dt,
                source="fm35",
            )

            if cached_wsgi is not None:
                _debug(
                    f"Skipping cached WSGI failure: {sounding_dt:%Y-%m-%d %H:%M} "
                    f"FM35 - {cached_wsgi.get('message', 'unknown failure')}"
                )
                failures.append(
                    f"{sounding_dt:%Y-%m-%d %H:%M} UTC WSGI/FM35: "
                    f"cached failure: {cached_wsgi.get('message', 'unknown failure')}"
                )
            else:
                attempted_network = True
                status = _download_one_wyoming_wsgi_source(
                    wmo_id=wmo_id,
                    dt=sounding_dt,
                    save_dir=save_dir,
                    session=session,
                    source="FM35",
                )

                if status.ok:
                    _debug("Workflow finished successfully using WSGI endpoint")
                    return status

                failures.append(
                    f"{sounding_dt:%Y-%m-%d %H:%M} UTC WSGI/FM35: "
                    f"{status.message}; url={status.url}"
                )
                _cache_failure(
                    cache=cache,
                    wmo_id=wmo_id,
                    dt=sounding_dt,
                    source="fm35",
                    message=status.message,
                    url=status.url or "",
                )
                cache_changed = True

            for source in ("bufr", "temp"):
                cached = _cached_failure(
                    cache=cache,
                    wmo_id=wmo_id,
                    dt=sounding_dt,
                    source=source,
                )

                if cached is not None:
                    _debug(
                        f"Skipping cached failure: {sounding_dt:%Y-%m-%d %H:%M} "
                        f"{source.upper()} - {cached.get('message', 'unknown failure')}"
                    )
                    failures.append(
                        f"{sounding_dt:%Y-%m-%d %H:%M} UTC {source.upper()}: "
                        f"cached failure: {cached.get('message', 'unknown failure')}"
                    )
                    continue

                attempted_network = True
                status = _download_one_wyoming_new_source(
                    wmo_id=wmo_id,
                    dt=sounding_dt,
                    save_dir=save_dir,
                    source=source,
                    session=session,
                )

                if status.ok:
                    _debug("Workflow finished successfully using NEW endpoint")
                    return status

                failures.append(
                    f"{sounding_dt:%Y-%m-%d %H:%M} UTC {source.upper()}: "
                    f"{status.message}; url={status.url}"
                )
                _cache_failure(
                    cache=cache,
                    wmo_id=wmo_id,
                    dt=sounding_dt,
                    source=source,
                    message=status.message,
                    url=status.url or "",
                )
                cache_changed = True

            cached_old = _cached_failure(
                cache=cache,
                wmo_id=wmo_id,
                dt=sounding_dt,
                source="old",
            )

            if cached_old is not None:
                _debug(
                    f"Skipping cached OLD failure: {sounding_dt:%Y-%m-%d %H:%M} - "
                    f"{cached_old.get('message', 'unknown failure')}"
                )
                failures.append(
                    f"{sounding_dt:%Y-%m-%d %H:%M} UTC OLD: "
                    f"cached failure: {cached_old.get('message', 'unknown failure')}"
                )
                continue

            attempted_network = True
            _debug("Both NEW sources failed or were skipped; trying OLD endpoint fallback")

            old_status = _download_one_wyoming_old_source(
                wmo_id=wmo_id,
                dt=sounding_dt,
                save_dir=save_dir,
                session=session,
            )

            if old_status.ok:
                _debug("Workflow finished successfully using OLD endpoint fallback")
                return old_status

            failures.append(
                f"{sounding_dt:%Y-%m-%d %H:%M} UTC OLD: "
                f"{old_status.message}; url={old_status.url}"
            )
            _cache_failure(
                cache=cache,
                wmo_id=wmo_id,
                dt=sounding_dt,
                source="old",
                message=old_status.message,
                url=old_status.url or "",
            )
            cache_changed = True

    finally:
        session.close()
        _debug("Closed HTTP session")

        if cache_changed:
            _save_failure_cache(save_dir=save_dir, cache=cache)

    if attempted_network:
        prefix = "Wyoming download failed for all attempted sources within +/-12 h."
    else:
        prefix = (
            "Wyoming download aborted because all candidate sources have cached "
            "failures within +/-12 h."
        )

    _debug("Workflow failed for all candidates")

    return DownloadStatus(
        ok=False,
        path=None,
        message=(
            f"{prefix} Requested time={request_dt}. Attempted sounding times="
            f"{', '.join(dt.strftime('%Y-%m-%d %H:%M') for dt in candidate_times)}. "
            + " | ".join(failures)
        ),
    )
# -----------------------------


# ----------- _main ----------
def _main() -> int:
    """Command line interface for debugging the Wyoming downloader."""
    global WYOMING_DEBUG

    parser = argparse.ArgumentParser(
        description="Debug and test the ATLAS Wyoming radiosonde downloader."
    )

    parser.add_argument("wmo", help="WMO station ID, e.g. 15420")
    parser.add_argument("date", help="Date, e.g. 27.05.2026 or 2026-05-27")
    parser.add_argument("time", help="UTC time, e.g. 08:30:00")
    parser.add_argument("--outdir", required=True, help="Output folder for downloaded radiosonde files")
    parser.add_argument("--debug", action="store_true", help="Print each workflow step in the CMD window")
    parser.add_argument("--clear-failure-cache", action="store_true", help="Delete .wyoming_download_failures.json before running")

    args = parser.parse_args()
    WYOMING_DEBUG = bool(args.debug)

    print("[CMD] Wyoming downloader test started")
    print(f"[CMD] WMO     : {args.wmo}")
    print(f"[CMD] Date    : {args.date}")
    print(f"[CMD] Time UTC: {args.time}")
    print(f"[CMD] Outdir  : {os.path.abspath(args.outdir)}")
    print(f"[CMD] Debug   : {WYOMING_DEBUG}")

    os.makedirs(args.outdir, exist_ok=True)

    if args.clear_failure_cache:
        cache_path = _failure_cache_path(args.outdir)
        print(f"[CMD] Clearing failure cache: {cache_path}")

        if os.path.isfile(cache_path):
            os.remove(cache_path)
            print("[CMD] Failure cache removed")
        else:
            print("[CMD] Failure cache did not exist")

    result = _download_wyoming(
        wmo_id=args.wmo,
        date=args.date,
        time_utc=args.time,
        save_dir=args.outdir,
    )

    print("[CMD] Wyoming downloader test finished")
    print("OK      :", result.ok)
    print("PATH    :", result.path)
    print("MESSAGE :", result.message)
    print("URL     :", result.url)

    return 0 if result.ok else 1
# -----------------------------


if __name__ == "__main__":
    raise SystemExit(_main())
