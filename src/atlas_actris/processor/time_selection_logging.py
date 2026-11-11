"""Concise time-selection diagnostics; one line per original INI triplet."""
from utils.printouts import print_entry


def _format_time(value):
    return value.strftime("%Y-%m-%d %H:%M:%S" if value.second else "%Y-%m-%d %H:%M")


def print_time_selections(action, events, requests=None):
    """Summarize expanded sector rules without hiding resolved calendar dates.

    Event counts concern coordinates only. Slice matches describe each interval,
    not the final union of several slice rules. Exclusion counts reflect rows
    actually removed from the running mask, after slicing and prior exclusions.
    """
    if requests is not None:
        expected = [(member, request["start"], request["stop"])
                    for request in requests for member in request["members"]]
        actual = [(e["key"], e["start_text"], e["stop_text"]) for e in events]
        if expected != actual:
            requests = None  # caller_info was edited after parsing.
    if requests is None:
        requests = [dict(key=e["key"], start=e["start_text"], stop=e["stop_text"],
                         members=[e["key"]]) for e in events]
    offset = 0
    for request in requests:
        count = len(request["members"])
        selected = events[offset:offset + count]
        offset += count
        key = request["key"]
        heading = f"{action} {key}"
        if not selected:
            print_entry(f"{heading}: {request['start']} -> {request['stop']}; skipped (overridden).")
            continue
        available = [event for event in selected if event["source"] is not None]
        if not available:
            print_entry(f"{heading}: {request['start']} -> {request['stop']}; skipped (no data).")
            continue
        # Several logical aliases can point to the same physical data.
        unique = {}
        for event in available:
            unique.setdefault((event["source"], event["start"], event["stop"]), event)
        available = list(unique.values())
        ranges = list(dict.fromkeys((_format_time(e["start"]), _format_time(e["stop"])) for e in available))
        if len(ranges) == 1:
            interval = f"{ranges[0][0]} -> {ranges[0][1]} UTC"
        else:
            interval = "; ".join(
                f"{e['key']}: {_format_time(e['start'])} -> {_format_time(e['stop'])} UTC"
                for e in available
            )
        if len(selected) == 1 and available[0]["source"] != key:
            heading += f" ({available[0]['source']})"
        total = sum(e["total"] for e in available)
        if action == "Slice":
            matched = sum(e["matched"] for e in available)
            outcome = f"{matched}/{total} profiles in interval"
            if matched == total:
                outcome += " (unchanged)"
        else:
            removed = sum(e["removed"] for e in available)
            outcome = f"{removed} {'profile' if removed == 1 else 'profiles'} removed"
            if removed == 0:
                outcome += " (unchanged)"
        absent = [e["key"] for e in selected if e["source"] is None]
        if absent:
            outcome += f"; no data: {', '.join(absent)}"
        print_entry(f"{heading}: {interval}; {outcome}.")
