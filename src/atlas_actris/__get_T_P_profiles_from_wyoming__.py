#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wyoming radiosonde downloader using the University of Wyoming upper-air
endpoint shown at https://weather.uwyo.edu/upperair/sounding.shtml.

The public function keeps the same API expected by find_radiosonde.py:
    _download_wyoming(wmo_id, date, time_utc, save_dir)

It returns an object with:
    ok, path, message, url

This module only downloads the closest available Wyoming CSV response within
+/-12 hours and saves it with the filename pattern already used by ATLAS:
    YYYYMMDD_HHMM_wyoming_<source>_<wmo_id>.txt

It does not parse, standardize, or re-export the radiosonde data. Your existing
reader can skip/handle the header later.
"""

import json
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple

import requests


WYOMING_NEW_ENDPOINT = "https://weather.uwyo.edu/cgi-bin/bufrraob.py"
WYOMING_FAILURE_CACHE = ".wyoming_download_failures.json"
WYOMING_USER_AGENT = "atlas-actris-radiosonde-downloader/1.0"

# Fast-fail timeout: separate connect/read values avoid long hangs when the
# Wyoming server is temporarily unavailable.
WYOMING_TIMEOUT = (3.05, 12.0)


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


def _wyoming_candidate_times(dt: datetime, search_hours: int = 12) -> List[datetime]:
    """Return nominal Wyoming sounding times within +/- search_hours.

    Wyoming offers the nominal hours visible on the upper-air page:
    00, 03, 06, 09, 12, 15, 18, and 21 UTC. The candidates are sorted by
    temporal distance from the requested measurement time.
    """

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

    return candidates


def _nearest_wyoming_hour(dt: datetime) -> datetime:
    """Snap to the nearest nominal Wyoming sounding hour."""

    candidates = _wyoming_candidate_times(dt=dt, search_hours=12)
    return min(candidates, key=lambda x: abs(x - dt))


def _wyoming_filename(wmo_id: str, dt: datetime, source: str) -> str:
    """Return the ATLAS-expected local filename for a Wyoming download."""

    return f"{dt:%Y%m%d_%H%M}_wyoming_{source}_{wmo_id}.txt"


def _existing_wyoming_file(wmo_id: str, dt: datetime, save_dir: str) -> Optional[str]:
    """Return an already-downloaded Wyoming file for this WMO/time, if present."""

    if not save_dir or not os.path.isdir(save_dir):
        return None

    prefix = f"{dt:%Y%m%d_%H%M}_wyoming_"
    suffix = f"_{wmo_id}.txt"

    for fname in sorted(os.listdir(save_dir)):
        if fname.startswith(prefix) and fname.endswith(suffix):
            path = os.path.join(save_dir, fname)
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                return path

    return None


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


def _failure_cache_path(save_dir: str) -> str:
    """Return the persistent failure cache path."""

    return os.path.join(save_dir, WYOMING_FAILURE_CACHE)


def _load_failure_cache(save_dir: str) -> Dict[str, Dict[str, str]]:
    """Load cached Wyoming failures. Invalid caches are ignored."""

    path = _failure_cache_path(save_dir)
    if not os.path.isfile(path):
        return {}

    try:
        with open(path, "r", encoding="utf-8") as fobj:
            data = json.load(fobj)
    except Exception:
        return {}

    if not isinstance(data, dict):
        return {}

    return data


def _save_failure_cache(save_dir: str, cache: Dict[str, Dict[str, str]]) -> None:
    """Atomically save the Wyoming failure cache."""

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


def _failure_key(wmo_id: str, dt: datetime, source: str) -> str:
    """Return a stable cache key for one station/time/source request."""

    return f"{wmo_id}|{dt:%Y%m%d%H%M}|{source}"


def _cache_failure(
    cache: Dict[str, Dict[str, str]],
    wmo_id: str,
    dt: datetime,
    source: str,
    message: str,
    url: str,
) -> None:
    """Store one failed station/time/source attempt in the cache."""

    cache[_failure_key(wmo_id=wmo_id, dt=dt, source=source)] = {
        "wmo_id": str(wmo_id),
        "datetime_utc": dt.strftime("%Y-%m-%d %H:%M:%S"),
        "source": str(source),
        "message": str(message),
        "url": str(url),
        "cached_at_utc": datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),
    }


def _cached_failure(
    cache: Dict[str, Dict[str, str]], wmo_id: str, dt: datetime, source: str
) -> Optional[Dict[str, str]]:
    """Return cached failure metadata for one request, if present."""

    return cache.get(_failure_key(wmo_id=wmo_id, dt=dt, source=source))


def _looks_like_failed_wyoming_response(content: bytes) -> Tuple[bool, str]:
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

    # Successful CSV responses from this endpoint should contain commas. If the
    # server returns an HTML page without CSV-like text, do not save it as data.
    if b"," not in content[:8192] and b"<html" in content[:8192].lower():
        return True, "Wyoming returned an HTML page instead of CSV data"

    return False, ""


def _write_bytes_atomic(path: str, content: bytes) -> None:
    """Write downloaded content atomically."""

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


def _download_one_wyoming_source(
    wmo_id: str,
    dt: datetime,
    save_dir: str,
    source: str,
    session: requests.Session,
) -> DownloadStatus:
    """Download one Wyoming source and save the raw CSV response."""

    url = _build_new_wyoming_url(wmo_id=wmo_id, dt=dt, source=source)

    try:
        response = session.get(url, timeout=WYOMING_TIMEOUT)
    except requests.RequestException as exc:
        return DownloadStatus(
            ok=False,
            message=f"Request to Wyoming failed: {exc}",
            url=url,
        )

    if not response.ok:
        return DownloadStatus(
            ok=False,
            message=f"HTTP {response.status_code} from Wyoming",
            url=url,
        )

    failed, failure_message = _looks_like_failed_wyoming_response(response.content)
    if failed:
        return DownloadStatus(ok=False, message=failure_message, url=url)

    fname = _wyoming_filename(wmo_id=wmo_id, dt=dt, source=source)
    path = os.path.join(save_dir, fname)
    _write_bytes_atomic(path=path, content=response.content)

    return DownloadStatus(
        ok=True,
        path=path,
        message=(
            f"Downloaded raw Wyoming {source.upper()} CSV sounding "
            f"for {dt:%Y-%m-%d %H:%M} UTC"
        ),
        url=url,
    )


def _download_wyoming(wmo_id, date, time_utc, save_dir):
    """
    Download the closest available Wyoming radiosonde CSV within +/-12 hours.

    Existing local files are reused first. Failed station/time/source attempts
    are cached in save_dir/.wyoming_download_failures.json; a later call skips
    cached failures and aborts immediately when every candidate source is known
    to have failed before.
    """

    request_dt = _parse_request_datetime(date=date, time_utc=time_utc)
    candidate_times = _wyoming_candidate_times(dt=request_dt, search_hours=12)
    wmo_id = str(wmo_id).strip()
    os.makedirs(save_dir, exist_ok=True)

    cache = _load_failure_cache(save_dir=save_dir)
    failures = []
    attempted_network = False
    cache_changed = False

    session = requests.Session()
    session.headers.update({"User-Agent": WYOMING_USER_AGENT})

    try:
        for sounding_dt in candidate_times:
            existing_path = _existing_wyoming_file(
                wmo_id=wmo_id,
                dt=sounding_dt,
                save_dir=save_dir,
            )

            if existing_path is not None:
                return DownloadStatus(
                    ok=True,
                    path=existing_path,
                    message=(
                        "Using already-downloaded Wyoming sounding "
                        f"for {sounding_dt:%Y-%m-%d %H:%M} UTC"
                    ),
                )

            for source in ("bufr", "temp"):
                cached = _cached_failure(
                    cache=cache,
                    wmo_id=wmo_id,
                    dt=sounding_dt,
                    source=source,
                )

                if cached is not None:
                    failures.append(
                        f"{sounding_dt:%Y-%m-%d %H:%M} UTC {source.upper()}: "
                        f"cached failure: {cached.get('message', 'unknown failure')}"
                    )
                    continue

                attempted_network = True
                status = _download_one_wyoming_source(
                    wmo_id=wmo_id,
                    dt=sounding_dt,
                    save_dir=save_dir,
                    source=source,
                    session=session,
                )

                if status.ok:
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
    finally:
        session.close()
        if cache_changed:
            _save_failure_cache(save_dir=save_dir, cache=cache)

    if attempted_network:
        prefix = "Wyoming download failed for all attempted sources within +/-12 h."
    else:
        prefix = (
            "Wyoming download aborted because all candidate sources have cached "
            "failures within +/-12 h."
        )

    return DownloadStatus(
        ok=False,
        path=None,
        message=(
            f"{prefix} Requested time={request_dt}. Attempted sounding times="
            f"{', '.join(dt.strftime('%Y-%m-%d %H:%M') for dt in candidate_times)}. "
            + " | ".join(failures)
        ),
    )
