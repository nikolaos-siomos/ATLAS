#!/usr/bin/env python3
"""
Fast downloader for one University of Wyoming radiosonde sounding.

CLI usage, exactly four positional arguments:

    python wyoming_radiosonde_downloader.py WMO_ID date time save-dir

Example:

    python wyoming_radiosonde_downloader.py 15420 01.01.2024 03:10:00 data

Import usage:

    from wyoming_radiosonde_downloader import _download_cloudnet_ecmwf
    result = _download_cloudnet_ecmwf("15420", "01.01.2024", "03:10:00", "data")

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import json
import re
import socket
import time as time_module
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


WYOMING_ENDPOINT = "https://weather.uwyo.edu/wsgi/sounding"
USER_AGENT = "wyoming-radiosonde-downloader/1.2"

SEARCH_WINDOW_HOURS = 6
STANDARD_SOUNDING_HOURS = (0, 6, 12, 18)
DATA_SOURCES = ("FM35", "BUFR")

# Main speed/safety tradeoff.
# Lower = faster failure when Wyoming hangs.
# Higher = more tolerant of slow Wyoming responses.
REQUEST_TIMEOUT_SECONDS = 5

# Exact filename stem requested:
#   <yyyymmdd>_<hhmm>_wyoming_<wmo_id>
# Set to ".csv" if you want a CSV extension.
OUTPUT_EXTENSION = ".dat"


@dataclass(frozen=True)
class DownloadResult:
    """Structured result returned by _download_cloudnet_ecmwf()."""

    ok: bool
    status: str
    message: str
    station: str
    requested_datetime_utc: Optional[str] = None
    sounding_datetime_utc: Optional[str] = None
    source: Optional[str] = None
    path: Optional[str] = None
    url: Optional[str] = None
    checked_urls: int = 0
    elapsed_seconds: Optional[float] = None
    exit_code: int = 0


@dataclass(frozen=True)
class FetchResult:
    """Result of one archive URL fetch."""

    content: Optional[bytes]
    error: Optional[str]
    classification: str


def parse_wmo_id(wmo_id: str) -> str:
    """Validate and normalize a WMO station ID."""
    value = str(wmo_id).strip()

    if not re.fullmatch(r"\d{5}", value):
        raise ValueError(
            f"Invalid WMO station id {wmo_id!r}. "
            "Use a five-digit numeric WMO id, e.g. 15420."
        )

    return value


def parse_requested_datetime(date: str, time_utc: str) -> datetime:
    """Parse dd.mm.yyyy and hh:mm:ss as a UTC datetime."""
    try:
        dt = datetime.strptime(f"{date} {time_utc}", "%d.%m.%Y %H:%M:%S")
    except ValueError as exc:
        raise ValueError(
            "Invalid date/time. Expected date as dd.mm.yyyy and time as hh:mm:ss UTC, "
            "for example: 01.01.2024 00:00:00"
        ) from exc

    return dt.replace(tzinfo=timezone.utc)


def candidate_datetimes(requested: datetime) -> list[datetime]:
    """
    Build candidate radiosonde datetimes within +/- 6 hours.

    Candidates use 00, 06, 12, and 18 UTC and are sorted by closeness to the
    requested datetime. Earlier times win exact ties.
    """
    window_start = requested - timedelta(hours=SEARCH_WINDOW_HOURS)
    window_end = requested + timedelta(hours=SEARCH_WINDOW_HOURS)

    candidates: list[datetime] = []
    current_date = window_start.date() - timedelta(days=1)
    end_date = window_end.date() + timedelta(days=1)

    while current_date <= end_date:
        for hour in STANDARD_SOUNDING_HOURS:
            candidate = datetime(
                current_date.year,
                current_date.month,
                current_date.day,
                hour,
                0,
                0,
                tzinfo=timezone.utc,
            )

            if window_start <= candidate <= window_end:
                candidates.append(candidate)

        current_date += timedelta(days=1)

    return sorted(candidates, key=lambda x: (abs(x - requested), x))


def build_url(wmo_id: str, sounding_dt: datetime, source: str) -> str:
    """Create a University of Wyoming CSV query URL."""
    params = {
        "datetime": sounding_dt.strftime("%Y-%m-%d %H:%M:%S"),
        "id": wmo_id,
        "src": source,
        "type": "TEXT:CSV",
    }

    return f"{WYOMING_ENDPOINT}?{urlencode(params)}"


def classify_sounding_content(content: Optional[bytes]) -> str:
    """
    Classify a Wyoming response body.

    Returns one of:
        valid
        station_not_found
        no_data
        invalid
    """
    if not content:
        return "invalid"

    text = content.decode("utf-8", errors="ignore").strip()

    if len(text) < 40:
        return "invalid"

    lowered = text.lower()

    station_failure_phrases = (
        "invalid station",
        "station not found",
        "unknown station",
        "invalid station number",
        "invalid wmo",
    )
    if any(phrase in lowered for phrase in station_failure_phrases):
        return "station_not_found"

    no_data_phrases = (
        "no data",
        "no observations",
        "can't get",
        "cannot get",
        "not available",
        "sorry",
    )
    if any(phrase in lowered for phrase in no_data_phrases):
        return "no_data"

    if "<html" in lowered or "</html>" in lowered:
        return "invalid"

    if "pressure" in lowered and "height" in lowered and "temperature" in lowered:
        return "valid"

    numeric_lines = 0
    for line in text.splitlines():
        if line.count(",") >= 4 and re.search(r"\d", line):
            numeric_lines += 1

    if numeric_lines >= 3:
        return "valid"

    return "invalid"


def fetch_url(url: str, timeout: int = REQUEST_TIMEOUT_SECONDS) -> FetchResult:
    """
    Fetch a URL and classify the response.

    There are no automatic retries. This keeps failure cases fast.
    """
    try:
        request = Request(url, headers={"User-Agent": USER_AGENT})

        with urlopen(request, timeout=timeout) as response:
            content = response.read()

        return FetchResult(
            content=content,
            error=None,
            classification=classify_sounding_content(content),
        )

    except HTTPError as exc:
        body = b""
        try:
            body = exc.read() or b""
        except Exception:
            body = b""

        classification = classify_sounding_content(body)
        if classification == "invalid":
            classification = "no_data"

        return FetchResult(
            content=body,
            error=f"http_error:{exc.code}",
            classification=classification,
        )

    except (URLError, socket.timeout, TimeoutError, OSError) as exc:
        return FetchResult(
            content=None,
            error=f"connection_error:{exc}",
            classification="invalid",
        )


def output_filename(sounding_dt: datetime, wmo_id: str) -> str:
    """
    Create output filename using the downloaded radiosonde datetime.

    Format:
        <yyyymmdd>_<hhmm>_wyoming_<wmo_id>
    """
    return f"{sounding_dt:%Y%m%d}_{sounding_dt:%H%M}_wyoming_{wmo_id}{OUTPUT_EXTENSION}"


def make_result(
    *,
    ok: bool,
    status: str,
    message: str,
    station: str,
    started: float,
    requested_dt: Optional[datetime] = None,
    sounding_dt: Optional[datetime] = None,
    source: Optional[str] = None,
    path: Optional[str] = None,
    url: Optional[str] = None,
    checked_urls: int = 0,
    exit_code: int = 0,
) -> DownloadResult:
    """Build a DownloadResult with elapsed time filled in."""
    return DownloadResult(
        ok=ok,
        status=status,
        message=message,
        station=station,
        requested_datetime_utc=requested_dt.isoformat() if requested_dt else None,
        sounding_datetime_utc=sounding_dt.isoformat() if sounding_dt else None,
        source=source,
        path=path,
        url=url,
        checked_urls=checked_urls,
        elapsed_seconds=round(time_module.perf_counter() - started, 3),
        exit_code=exit_code,
    )


def download_closest_sounding(
    wmo_id: str,
    requested_dt: datetime,
    save_dir: str | Path,
    *,
    timeout: int = REQUEST_TIMEOUT_SECONDS,
) -> DownloadResult:
    """
    Download the closest available sounding within +/- 6 hours.

    This version is staged rather than fully parallel:
      - try the closest candidate time first;
      - for that time, try FM35 and then BUFR;
      - return immediately after a valid sounding is found;
      - do not write a manifest file.

    This usually beats the fully parallel version because it does not wait for
    slow, farther-away candidate URLs after a closer valid sounding has already
    been found.
    """
    started = time_module.perf_counter()
    station = parse_wmo_id(wmo_id)
    output_root = Path(save_dir)
    candidates = candidate_datetimes(requested_dt)

    if not candidates:
        return make_result(
            ok=False,
            status="no_radiosonde",
            message="No standard radiosonde time exists within +/- 6 hours.",
            station=station,
            requested_dt=requested_dt,
            started=started,
            exit_code=3,
        )

    checked_urls = 0
    saw_connection_error = False
    saw_any_non_connection_response = False

    for sounding_dt in candidates:
        for source in DATA_SOURCES:
            url = build_url(station, sounding_dt, source)
            checked_urls += 1

            fetched = fetch_url(url, timeout=timeout)

            if fetched.error and fetched.error.startswith("connection_error"):
                saw_connection_error = True
                continue

            saw_any_non_connection_response = True

            if fetched.classification == "station_not_found":
                return make_result(
                    ok=False,
                    status="station_not_found",
                    message=(
                        "The archive response indicates that this WMO station was not found. "
                        "Check that the WMO ID is correct and available in the Wyoming archive."
                    ),
                    station=station,
                    requested_dt=requested_dt,
                    started=started,
                    url=url,
                    checked_urls=checked_urls,
                    exit_code=3,
                )

            if fetched.classification == "valid" and fetched.content:
                try:
                    output_root.mkdir(parents=True, exist_ok=True)
                    filename = output_filename(sounding_dt, station)
                    path = output_root / filename
                    path.write_bytes(fetched.content)

                except OSError as exc:
                    return make_result(
                        ok=False,
                        status="error",
                        message=f"File-system error: {exc}",
                        station=station,
                        requested_dt=requested_dt,
                        started=started,
                        url=url,
                        checked_urls=checked_urls,
                        exit_code=4,
                    )

                return make_result(
                    ok=True,
                    status="downloaded",
                    message="Downloaded the closest available radiosonde sounding.",
                    station=station,
                    requested_dt=requested_dt,
                    sounding_dt=sounding_dt,
                    source=source,
                    path=str(path),
                    url=url,
                    checked_urls=checked_urls,
                    started=started,
                    exit_code=0,
                )

    if saw_connection_error and not saw_any_non_connection_response:
        return make_result(
            ok=False,
            status="connection_error",
            message=(
                "Could not contact the University of Wyoming archive. Check the internet "
                "connection, DNS, proxy/firewall settings, or try again later."
            ),
            station=station,
            requested_dt=requested_dt,
            checked_urls=checked_urls,
            started=started,
            exit_code=2,
        )

    return make_result(
        ok=False,
        status="no_radiosonde",
        message=(
            "No radiosonde was found for this WMO id within +/- 6 hours of the "
            "requested date/time. The station may not exist in the Wyoming archive, "
            "or it may have no sounding near the requested time."
        ),
        station=station,
        requested_dt=requested_dt,
        checked_urls=checked_urls,
        started=started,
        exit_code=3,
    )


def _download_wyoming(
    wmo_id: str,
    date: str,
    time_utc: str,
    save_dir: str | Path,
) -> DownloadResult:
    """
    Callable API for other Python code.

    Parameters
    ----------
    wmo_id:
        Five-digit WMO station id, e.g. "15420".
    date:
        Date as dd.mm.yyyy, e.g. "01.01.2024".
    time_utc:
        Time as hh:mm:ss UTC, e.g. "03:10:00".
    save_dir:
        Directory where the downloaded file is written.

    Returns
    -------
    DownloadResult
        Structured success/failure result. Expected errors are returned, not raised.
    """
    started = time_module.perf_counter()

    try:
        station = parse_wmo_id(wmo_id)
        requested_dt = parse_requested_datetime(date, time_utc)

        return download_closest_sounding(
            wmo_id=station,
            requested_dt=requested_dt,
            save_dir=save_dir,
        )

    except ValueError as exc:
        return make_result(
            ok=False,
            status="invalid_input",
            message=str(exc),
            station=str(wmo_id),
            started=started,
            exit_code=1,
        )

    except OSError as exc:
        return make_result(
            ok=False,
            status="error",
            message=f"File-system error: {exc}",
            station=str(wmo_id),
            started=started,
            exit_code=4,
        )

    except Exception as exc:
        return make_result(
            ok=False,
            status="error",
            message=f"Unexpected error: {exc}",
            station=str(wmo_id),
            started=started,
            exit_code=4,
        )


def build_arg_parser() -> argparse.ArgumentParser:
    """Create CLI parser outside _download_cloudnet_ecmwf(), so _download_cloudnet_ecmwf() remains callable."""
    parser = argparse.ArgumentParser(
        description="Download one closest University of Wyoming radiosonde sounding."
    )

    parser.add_argument("wmo_id", help="five-digit WMO station id, e.g. 15420")
    parser.add_argument("date", help='date as "dd.mm.yyyy", e.g. 01.01.2024')
    parser.add_argument("time", help='time as "hh:mm:ss" UTC, e.g. 03:10:00')
    parser.add_argument("save_dir", help="directory where the downloaded file is saved")

    return parser


def cli(argv: Optional[list[str]] = None) -> int:
    """CLI entry point."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    result = _download_cloudnet_ecmwf(
        wmo_id=args.wmo_id,
        date=args.date,
        time_utc=args.time,
        save_dir=args.save_dir,
    )

    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))
    return result.exit_code


if __name__ == "__main__":
    raise SystemExit(cli())
