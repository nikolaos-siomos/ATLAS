#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:13:42 2026

@author: nikos
"""

"""Find common measurement intervals in the IR campaign ODS protocol."""

import os
import warnings
import subprocess

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"

def common_intervals(path, test_type, period, date, *, systems=None,
                     exclude_systems=None, max_gap_minutes=15,
                     min_duration_minutes=15):
    """Return [(start, end), ...] as UTC strings in YYYYMMDD_HHMM format.

    Uses only Python's standard library. Reads the IR_campaign sheet with date
    labels in A. Detects three-column (start time, stop time, type) or five-column
    (start date, start time, stop date, stop time, type) groups from the headers.
    Matches the exact test label (ignoring surrounding whitespace).
    period is 'Daytime' or 'Nighttime' (case-insensitive). Only rows under that
    label in column A are included, until the next label or date header.
    date is the required block date as DD.MM.YYYY, e.g. '17.09.2026'. Only that
    block is searched; matching overnight intervals can end on the next day.
    Invalid measurements in other date blocks are not evaluated.

    Compares only systems with entries for the requested test, date, and period.
    Systems without matching entries are ignored, even if explicitly selected.
    Set systems to a sequence of header names to compare specific systems, e.g.
    ("POLIS1064", "M3DUSA") in the newer file. Names must match the file's headers.
    One remaining system returns its own intervals; no matching systems returns [].
    Matching entries with missing/invalid dates or times still raise ValueError.
    exclude_systems is an optional sequence of header names to leave out. If
    systems is omitted, considers all recorded systems except those excluded.
    If both are supplied, exclusions are removed from systems (exclusion wins).
    Unknown names or an empty final selection raise ValueError.

    After finding the common intervals, merges gaps strictly shorter than
    max_gap_minutes (default 15). A gap of exactly 15 minutes is not merged.
    Then removes intervals shorter than min_duration_minutes (default 15);
    exactly 15 minutes is retained. Duration is the merged start-to-end span,
    including bridged gaps, not just simultaneous measurement time.
    Set either parameter to 0 to disable that behavior. Both must be finite,
    nonnegative numbers of minutes. Overlapping/touching intervals always merge.
    Merging is applied to final common intervals, not individual system records.
    Daytime and Nighttime sections
    belong to their preceding date. Explicit start/stop dates take precedence.
    A blank start date uses the block date. A blank stop date uses the start date,
    adding one day if the stop time is earlier than the start time. Equal times
    on the same date mean zero duration. Each date block is evaluated separately.
    Dates may be ODS date values or DD.MM.YYYY, YYYY.MM.DD, or YYYY-MM-DD text.
    If both date cells contain times and both time cells are blank, reads those
    misplaced times with a warning and uses the block date. Does not edit the ODS.

    Raises ValueError for missing/invalid dates or times on matching rows, unknown system
    names, or an unrecognized sheet layout. Returns [] if no overlap exists or
    the requested date/period is absent.
    Does not infer interruptions from rows of other test types or mixed labels
    such as 'ray+drk sequence in pcb'. Input times are assumed to be UTC.
    """
    from datetime import datetime, timedelta
    import re
    import warnings
    import math
    from xml.etree import ElementTree as ET
    from zipfile import ZipFile

    for name, value in (("max_gap_minutes", max_gap_minutes),
                        ("min_duration_minutes", min_duration_minutes)):
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or value < 0):
            raise ValueError(f"{name} must be a finite, nonnegative number of minutes")
    if not isinstance(date, str) or not re.fullmatch(r"\d{2}\.\d{2}\.\d{4}", date.strip()):
        raise ValueError("date must use DD.MM.YYYY, e.g. '17.09.2026'")
    try:
        target_date = datetime.strptime(date.strip(), "%d.%m.%Y")
    except ValueError as exc:
        raise ValueError(f"Invalid date: {date!r}; use DD.MM.YYYY") from exc
    if not isinstance(test_type, str) or not test_type.strip():
        raise ValueError("test_type must be a nonempty string")
    test_type = test_type.strip()
    if not isinstance(period, str) or period.strip().lower() not in ("daytime", "nighttime"):
        raise ValueError("period must be 'Daytime' or 'Nighttime'")
    period = period.strip().lower()
    ns = {
        "t": "urn:oasis:names:tc:opendocument:xmlns:table:1.0",
        "x": "urn:oasis:names:tc:opendocument:xmlns:text:1.0",
        "o": "urn:oasis:names:tc:opendocument:xmlns:office:1.0",
    }
    with ZipFile(path) as ods:
        root = ET.fromstring(ods.read("content.xml"))
    sheet = next((s for s in root.findall(".//t:table", ns)
                  if s.get(f"{{{ns['t']}}}name") == "IR_campaign"), None)
    if sheet is None:
        raise ValueError("Sheet 'IR_campaign' was not found")

    def text(cell):
        return " ".join("".join(p.itertext())
                        for p in cell.findall("x:p", ns)).strip()

    def has_value(cell):
        return bool(text(cell) or cell.get(f"{{{ns['o']}}}time-value")
                    or cell.get(f"{{{ns['o']}}}date-value"))

    def calendar_date(cell):
        value = cell.get(f"{{{ns['o']}}}date-value")
        value = value[:10] if value else text(cell)
        if not value:
            return None
        for fmt in ("%Y-%m-%d", "%d.%m.%Y", "%Y.%m.%d"):
            try:
                return datetime.strptime(value, fmt)
            except ValueError:
                pass
        raise ValueError(f"invalid date: {value!r}")

    def merge(intervals, gap_minutes=0):
        result = []
        for start, end in sorted(intervals):
            if start >= end:
                continue
            if result and (start <= result[-1][1]
                           or (start - result[-1][1]).total_seconds() / 60 < gap_minutes):
                result[-1] = (result[-1][0], max(result[-1][1], end))
            else:
                result.append((start, end))
        return result

    def clock(cell):
        # Prefer the underlying ODS time value over its displayed formatting.
        value = cell.get(f"{{{ns['o']}}}time-value")
        if value:
            match = re.fullmatch(
                r"PT(?:(\d+)H)?(?:(\d+)M)?(?:(\d+(?:\.\d+)?)S)?", value)
            if match:
                h, m, s = (float(v or 0) for v in match.groups())
                if h < 24 and m < 60 and s < 60:
                    return timedelta(hours=h, minutes=m, seconds=s)
        else:
            for fmt in ("%H:%M", "%H:%M:%S"):
                try:
                    t = datetime.strptime(text(cell), fmt)
                    return timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)
                except ValueError:
                    pass
        raise ValueError("missing or invalid time")

    groups = None
    header = None
    selected = None
    blocks = []
    block = None
    current_period = None
    row_number = 0
    empty = ET.Element("empty")
    for row in sheet.findall("t:table-row", ns):
        row_number += 1
        cells = {}
        col = 0
        for cell in row:
            if cell.tag not in (f"{{{ns['t']}}}table-cell",
                                f"{{{ns['t']}}}covered-table-cell"):
                continue
            repeat = int(cell.get(f"{{{ns['t']}}}number-columns-repeated", "1"))
            if has_value(cell):
                for index in range(col, col + repeat):
                    cells[index] = cell
            col += repeat
        first = cells.get(0, empty)
        if header is None:
            if text(first).lower() == "date":
                header = cells
        elif groups is None:
            labels = {i: text(c).lower() for i, c in cells.items()}
            starts = sorted(i for i, label in labels.items()
                            if label == "start time (utc)")
            if not starts:
                row_number += int(row.get(f"{{{ns['t']}}}number-rows-repeated", "1")) - 1
                continue
            groups = []
            for start_col in starts:
                dated = labels.get(start_col - 1) == "start date"
                first_col = start_col - 1 if dated else start_col
                end_col = start_col + (2 if dated else 1)
                type_col = end_col + 1
                if (labels.get(end_col) != "stop time (utc)"
                        or labels.get(type_col) != "type"
                        or (dated and labels.get(start_col + 1) != "stop date")):
                    raise ValueError(f"Unrecognized measurement columns near column {start_col + 1}")
                # Some older headers split a system name over adjacent cells.
                name = "".join(dict.fromkeys(text(header.get(i, empty))
                                            for i in range(first_col, type_col + 1)))
                if not name:
                    raise ValueError(f"Missing system name at column {first_col + 1}")
                groups.append((name, start_col, end_col, type_col,
                               first_col if dated else None,
                               start_col + 1 if dated else None))
            names = [group[0] for group in groups]
            if len(set(names)) != len(names):
                raise ValueError("Duplicate system names in header")
            def system_names(value, argument, default):
                if value is None:
                    return set(default)
                if isinstance(value, str):
                    raise ValueError(f"{argument} must be a sequence of names, not a string")
                try:
                    values = list(value)
                except TypeError as exc:
                    raise ValueError(f"{argument} must be a sequence of names") from exc
                if any(not isinstance(name, str) for name in values):
                    raise ValueError(f"{argument} must contain only system names as strings")
                requested = set(values)
                unknown = requested - set(names)
                if unknown:
                    raise ValueError(
                        f"Unknown names in {argument}: {', '.join(sorted(unknown))}. "
                        f"Choose from: {', '.join(names)}"
                    )
                return requested

            included = system_names(systems, "systems", names)
            excluded = system_names(exclude_systems, "exclude_systems", ())
            selected = included - excluded
            if not selected:
                raise ValueError("No systems remain after applying systems and exclude_systems")
        else:
            label = text(first).lower()
            date = calendar_date(first) if has_value(first) and label not in ("daytime", "nighttime") else None
            if date is not None:
                block = (date, {}) if date == target_date else None
                if block is not None:
                    blocks.append(block)
                current_period = None
            if label in ("daytime", "nighttime"):
                current_period = label
            if block is not None and current_period == period:
                date, measurements = block
                for name, start_col, end_col, type_col, start_date_col, end_date_col in groups:
                    start_cell = cells.get(start_col, empty)
                    end_cell = cells.get(end_col, empty)
                    type_cell = cells.get(type_col, empty)
                    start_date_cell = cells.get(start_date_col, empty)
                    end_date_cell = cells.get(end_date_col, empty)
                    if name not in selected or not any(has_value(c) for c in
                        (start_cell, end_cell, type_cell, start_date_cell, end_date_cell)):
                        continue
                    if text(type_cell) == test_type:
                        intervals = measurements.setdefault(name, [])
                        try:
                            if (start_date_col is not None
                                    and not has_value(start_cell) and not has_value(end_cell)
                                    and has_value(start_date_cell) and has_value(end_date_cell)
                                    and not start_date_cell.get(f"{{{ns['o']}}}date-value")
                                    and not end_date_cell.get(f"{{{ns['o']}}}date-value")):
                                # Recognize only a complete pair of misplaced time values.
                                clock(start_date_cell)
                                clock(end_date_cell)
                                warnings.warn(
                                    f"{date:%d.%m.%Y}, {name}, row {row_number}: "
                                    "times are in Start Date/Stop Date columns; "
                                    "using these as times with the block date.",
                                    UserWarning, stacklevel=2,
                                )
                                start_cell, end_cell = start_date_cell, end_date_cell
                                start_date_cell = end_date_cell = empty
                            start_date = calendar_date(start_date_cell) or date
                            explicit_end_date = calendar_date(end_date_cell)
                            start = start_date + clock(start_cell)
                            end = (explicit_end_date or start_date) + clock(end_cell)
                            if explicit_end_date is None and end < start:
                                end += timedelta(days=1)
                            if end < start:
                                raise ValueError("explicit stop date/time is earlier than start date/time")
                        except ValueError as exc:
                            raise ValueError(
                                f"{date:%d.%m.%Y}, {name}, row {row_number}: {exc}"
                            ) from exc
                        intervals.append((start, end))
        row_number += int(row.get(f"{{{ns['t']}}}number-rows-repeated", "1")) - 1

    if groups is None:
        raise ValueError("Could not find system headers and Start Time (UTC) columns")
    result = []
    for date, measurements in blocks:
        required = set(measurements)
        if not required:
            continue
        common = None
        for name in sorted(required):
            intervals = merge(measurements[name])
            common = intervals if common is None else merge(
                (max(a, c), min(b, d))
                for a, b in common for c, d in intervals
                if max(a, c) < min(b, d)
            )
        common = merge(common, gap_minutes=max_gap_minutes)
        result.extend((start, end) for start, end in common
                      if (end - start).total_seconds() / 60 >= min_duration_minutes)
    return [(start.strftime("%Y%m%d_%H%M"), end.strftime("%Y%m%d_%H%M"))
            for start, end in sorted(result)]


#------------------------------------------------------------------------------
#MAIN
#------------------------------------------------------------------------------
dates = common_intervals(
    path = '/home/nikos/Nextcloud5/CARS-Room-Lidar/Campaigns/2026_IR_Campaign/data/IR_campaign_intercomparison_protocol.ods',
    test_type='ray_pcb',
    date = '17.09.2026',
    period= 'Nighttime',
    exclude_systems=["ALPHA"],
    )

init_files = [
    '/home/nikos/Nextcloud5/CARS-Room-Lidar/Campaigns/2026_IR_Campaign/data/POLIS/call_atlas_POL_9_N_20260917.ini',
    '/home/nikos/Nextcloud5/CARS-Room-Lidar/Campaigns/2026_IR_Campaign/data/PollyXT/call_atlas_PXT_N_20260917.ini',
    '/home/nikos/Nextcloud5/CARS-Room-Lidar/Campaigns/2026_IR_Campaign/data/RALI/call_atlas_RAL_N_1136_20260917.ini',
    '/home/nikos/Nextcloud5/CARS-Room-Lidar/Campaigns/2026_IR_Campaign/data/Emoral/call_atlas_EMO_N_1028_20260917.ini',
    ]

dir_out = '/home/nikos/Big_Data/Intercomparison_IR_1064/intercomparison'

for dt in dates:
    for init_file in init_files:
        subprocess.run([
            "atlas-prepare-intercomparison", "-i", init_file,
            '-s','ray_pcb', dt[0], dt[1],
            '-o', dir_out
            ], check=True)
