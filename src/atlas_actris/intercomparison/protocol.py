import subprocess

def common_intervals(path, test_type, period, date, *, systems=None,
                     exclude_systems=None, min_duration_minutes=0, verbose=True):
    """Return [(common_start, common_end), (gap_start, gap_end), ...] in UTC.

    date: DD.MM.YYYY block label; period: Daytime or Nighttime.
    ray_pcb uses each system's ray_pcb entries if present, otherwise its ray
    entries. If no included system has ray_pcb, prints an exit message and
    returns [] normally (also in IPython). ray uses only ray entries.
    Other test types match exactly. The caller should skip processing when [].
    Systems without eligible entries are ignored. systems includes names;
    exclude_systems removes names, taking precedence. Names match the headers.

    The common envelope is latest system start to earliest system end. Each
    system's envelope spans all its selected entries. Gaps are found between
    sorted, overlapping/touching-merged entries of each system, clipped to the
    common envelope and returned chronologically, WITHOUT union/deduplication
    across systems. The printed report also shows full, unclipped system gaps.
    No gaps are filled or filtered by duration. min_duration_minutes (default 0)
    optionally filters the common envelope by elapsed duration, including gaps.
    The old max_gap_minutes argument is no longer used/supported.

    Prints each system's choice, source intervals, envelope and gaps by default;
    verbose=False suppresses the report. Missing times in selected entries raise
    ValueError; unused ray entries do not affect a system using ray_pcb.
    No positive common envelope returns []. Overlapping gaps may cover the whole
    envelope: callers must apply their union to determine actual usable time.

    Reads both three- and five-column ODS layouts. Explicit measurement dates
    take priority. Blank start dates use the block date; blank stop dates use
    the start date, plus one day if stop time is earlier. Dates accept ODS values,
    DD.MM.YYYY, YYYY.MM.DD, YYYY-MM-DD. Misplaced pairs of times in date columns
    are recognized only when both time columns are blank, with a warning.
    Uses only the standard library; does not modify the workbook or run atlas.
    """
    from datetime import datetime, timedelta
    import re
    import warnings
    import math
    from xml.etree import ElementTree as ET
    from zipfile import ZipFile

    for name, value in (("min_duration_minutes", min_duration_minutes),):
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

    def merge(intervals):
        result = []
        for start, end in sorted(intervals):
            if start >= end:
                continue
            if result and start <= result[-1][1]:
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
                    entry_type = text(type_cell)
                    candidates = ("ray_pcb", "ray") if test_type == "ray_pcb" else (test_type,)
                    if entry_type in candidates:
                        measurements.setdefault(name, {}).setdefault(entry_type, []).append(
                            (row_number, start_cell, end_cell, start_date_cell,
                             end_date_cell, start_date_col)
                        )
        row_number += int(row.get(f"{{{ns['t']}}}number-rows-repeated", "1")) - 1

    if groups is None:
        raise ValueError("Could not find system headers and Start Time (UTC) columns")
    # Aggregate matching date sections before selecting each system's type.
    entries = {}
    for block_date, measurements in blocks:
        for name, types in measurements.items():
            for kind, rows in types.items():
                entries.setdefault(name, {}).setdefault(kind, []).extend(rows)

    def report(message):
        if verbose:
            print(message, flush=True)

    def stamp(value):
        return value.strftime("%Y%m%d_%H%M")

    def spans(values):
        return "; ".join(f"{stamp(a)} -> {stamp(b)}" for a, b in values) or "none"

    report(f"{target_date:%d.%m.%Y} {period.title()} | requested={test_type} | UTC")
    if test_type == "ray_pcb" and not any(v.get("ray_pcb") for v in entries.values()):
        print(f"No ray_pcb entries for any included system in "
              f"{target_date:%d.%m.%Y} {period.title()}. Nothing to process.", flush=True)
        return []

    resolved = {}
    date = target_date
    for name in names:
        if name not in selected:
            report(f"{name}: excluded by system selection")
            continue
        types = entries.get(name, {})
        kind = ("ray_pcb" if types.get("ray_pcb") else "ray") if test_type == "ray_pcb" else test_type
        if not types.get(kind):
            report(f"{name}: excluded (no eligible entries)")
            continue
        intervals = []
        for row_number, start_cell, end_cell, start_date_cell, end_date_cell, start_date_col in types[kind]:
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
            if end <= start:
                raise ValueError(f"{name}, row {row_number}: selected measurement has zero duration")
            intervals.append((start, end))
        merged = merge(intervals)
        gaps = [(a[1], b[0]) for a, b in zip(merged, merged[1:])]
        resolved[name] = (merged[0][0], merged[-1][1], gaps)
        fallback = " (fallback from ray_pcb)" if kind != test_type else ""
        report(f"{name}: {kind}{fallback}")
        report(f"  Entries: {spans(sorted(intervals))}")
        report(f"  Start/end: {stamp(merged[0][0])} -> {stamp(merged[-1][1])}")
        report(f"  Gaps: {spans(gaps)}")

    if not resolved:
        report("No eligible systems; returning [].")
        return []
    start = max(v[0] for v in resolved.values())
    end = min(v[1] for v in resolved.values())
    if start >= end:
        report("No common outer interval; returning [].")
        return []
    if (end - start).total_seconds() / 60 < min_duration_minutes:
        report("Common outer interval is shorter than min_duration_minutes; returning [].")
        return []
    gaps = sorted((max(a, start), min(b, end))
                  for _, _, system_gaps in resolved.values() for a, b in system_gaps
                  if max(a, start) < min(b, end))
    report(f"Common start/end: {stamp(start)} -> {stamp(end)}")
    report(f"Returned gaps (clipped, not unioned): {spans(gaps)}")
    return [(stamp(start), stamp(end))] + [(stamp(a), stamp(b)) for a, b in gaps]

def iterate_by_protocol(dates, init_files, dir_out: str):
    """Run preparation for each init file, or return normally if dates is empty.

    dates[0] is the common start/end pair; subsequent pairs are gaps.
    ATLAS expects each exclusion as a test/start/end triplet.
    """
    if not dates:
        print("No common dates. Nothing to process.", flush=True)
        return

    start, end = dates[0]
    exclude_args = []
    if len(dates) > 1:
        exclude_args.append("-e")
        for test_type in ("ray", "ray_pcb"):
            for gap_start, gap_end in dates[1:]:
                exclude_args.extend([test_type, gap_start, gap_end])

    for init_file in init_files:
        subprocess.run([
            "atlas-prepare-intercomparison", "-i", init_file,
            '-s','ray_pcb', start, end,'ray', start, end
            ] + exclude_args + [
            '-o', dir_out, '-q', 'ray', 'ray_pcb' 
            ], check=True)
