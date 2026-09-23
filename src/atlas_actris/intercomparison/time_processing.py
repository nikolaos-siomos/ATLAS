#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Temporal filtering and averaging for ATLAS intercomparison bundles.

The implementation follows the user-facing ATLAS slice/exclude semantics while
operating directly on the already prepared intercomparison bundle arrays.

Processing order for every channel/pair dataset:
1. If ``time`` is absent, leave the array unchanged.
2. If ``time`` is singleton, squeeze/drop it; slicing is intentionally skipped.
3. Otherwise apply global ``slice_measurement`` intervals, then
   ``exclude_measurement`` intervals.
4. Average the remaining signal/ratio over ``time``.
5. Average its uncertainty following the existing ATLAS convention:
   mean(error) / sqrt(N), where N is the number of finite time samples.
"""

from __future__ import annotations

from copy import copy
from datetime import datetime
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import xarray as xr


_ABSOLUTE_TIME_FORMATS = (
    "%Y%m%d",
    "%Y%m%d_%H",
    "%Y%m%d_%H%M",
    "%Y%m%d_%H%M%S",
)


def _copy_nested(value: Any) -> Any:
    """Copy dict containers without copying the underlying xarray/dask data."""

    if isinstance(value, dict):
        return {key: _copy_nested(item) for key, item in value.items()}
    return value


def _is_hhmm(value: str) -> bool:
    text = str(value).strip()
    return len(text) == 4 and text.isdigit()


def _parse_absolute_time(value: str) -> np.datetime64:
    text = str(value).strip()
    for fmt in _ABSOLUTE_TIME_FORMATS:
        try:
            return np.datetime64(datetime.strptime(text, fmt))
        except ValueError:
            continue
    raise ValueError(f"Unsupported absolute time value: {value!r}")


def _hhmm_minutes(value: str) -> int:
    text = str(value).strip()
    hour = int(text[:2])
    minute = int(text[2:4])
    return hour * 60 + minute


def _time_of_day_minutes(time: xr.DataArray) -> xr.DataArray:
    """Return minute-of-day for a datetime64 time coordinate."""

    return time.dt.hour * 60 + time.dt.minute


def _interval_mask(time: xr.DataArray, start: str, stop: str) -> xr.DataArray:
    """Build one inclusive ATLAS-style temporal interval mask.

    HHMM intervals are interpreted against the time-of-day of the measurements.
    If stop is earlier than start, the interval crosses midnight. Absolute
    values are compared directly as datetime64 values.
    """

    start = str(start).strip()
    stop = str(stop).strip()

    if _is_hhmm(start) and _is_hhmm(stop):
        start_minute = _hhmm_minutes(start)
        stop_minute = _hhmm_minutes(stop)
        minute_of_day = _time_of_day_minutes(time)

        if stop_minute >= start_minute:
            return (minute_of_day >= start_minute) & (minute_of_day <= stop_minute)

        # Crossing midnight, e.g. 2300 -> 0135.
        return (minute_of_day >= start_minute) | (minute_of_day <= stop_minute)

    if _is_hhmm(start) != _is_hhmm(stop):
        raise ValueError(
            "A slice/exclude interval must use HHMM for both boundaries or "
            "absolute date/time formats for both boundaries. "
            f"Got {start!r}, {stop!r}."
        )

    start_time = _parse_absolute_time(start)
    stop_time = _parse_absolute_time(stop)
    if stop_time < start_time:
        raise ValueError(
            f"Absolute time interval stop {stop!r} is earlier than start {start!r}. "
            "For absolute intervals, provide the next date explicitly."
        )

    return (time >= start_time) & (time <= stop_time)


def _combine_interval_masks(
    time: xr.DataArray,
    intervals: Sequence[tuple[str, str]],
) -> xr.DataArray:
    """Return the union of all supplied intervals."""

    mask = xr.zeros_like(time, dtype=bool)
    for start, stop in intervals:
        mask = mask | _interval_mask(time, start, stop)
    return mask


def _filter_time(
    value: xr.DataArray,
    *,
    slice_intervals: Sequence[tuple[str, str]],
    exclude_intervals: Sequence[tuple[str, str]],
    location: str,
) -> xr.DataArray:
    """Apply slice first and exclusions second along ``time``."""

    if "time" not in value.dims:
        return value

    if value.sizes["time"] <= 1:
        return value

    if "time" not in value.coords:
        raise ValueError(f"{location}: array has a 'time' dimension but no time coordinate")

    time = value["time"]
    if not np.issubdtype(time.dtype, np.datetime64):
        raise ValueError(
            f"{location}: time coordinate must be datetime64 for temporal slicing; "
            f"got dtype {time.dtype}."
        )

    keep = xr.ones_like(time, dtype=bool)

    if slice_intervals:
        keep = _combine_interval_masks(time, slice_intervals)

    if exclude_intervals:
        exclude = _combine_interval_masks(time, exclude_intervals)
        keep = keep & ~exclude

    filtered = value.where(keep, drop=True)
    if filtered.sizes.get("time", 0) == 0:
        raise ValueError(
            f"{location}: slice_measurement/exclude_measurement removed all time samples."
        )

    return filtered




def _diagnose_signal_stage(value: Any, *, label: str) -> None:
    """Print finite-value coverage without changing the processing result."""
    if not isinstance(value, xr.DataArray):
        print(f"        {label}: not an xarray.DataArray ({type(value).__name__})")
        return

    finite = value.notnull()
    total = int(value.size)
    finite_total = int(finite.sum().compute().item())
    print(f"        {label}: dims={value.dims}, shape={value.shape}, dtype={value.dtype}")
    print(f"        {label}: finite values={finite_total}/{total}")

    if "time" in value.dims:
        other_dims = [dim for dim in value.dims if dim != "time"]
        if other_dims:
            finite_per_time = finite.any(dim=other_dims)
            n_times_with_data = int(finite_per_time.sum().compute().item())
            print(
                f"        {label}: time samples containing at least one finite value="
                f"{n_times_with_data}/{value.sizes['time']}"
            )
        try:
            t0 = value.time.values[0]
            t1 = value.time.values[-1]
            print(f"        {label}: time extent={t0} -> {t1}")
        except Exception:
            pass
    else:
        if value.ndim:
            finite_any = finite
            n_finite = int(finite_any.sum().compute().item())
            print(f"        {label}: finite profile bins={n_finite}/{value.size}")


def _mean_signal(
    value: Any,
    *,
    slice_intervals: Sequence[tuple[str, str]],
    exclude_intervals: Sequence[tuple[str, str]],
    location: str,
) -> tuple[Any, dict[str, Any]]:
    """Filter and collapse a signal/ratio to a time-independent array."""

    if not isinstance(value, xr.DataArray) or "time" not in value.dims:
        return value, {"mode": "no_time", "n_before": None, "n_after": None}

    n_before = int(value.sizes["time"])
    if n_before == 1:
        return value.squeeze("time", drop=True), {
            "mode": "singleton",
            "n_before": 1,
            "n_after": 1,
        }

    print(f"\n    TIME-PROCESSING DIAGNOSTIC: {location}")
    _diagnose_signal_stage(value, label="before filtering")

    filtered = _filter_time(
        value,
        slice_intervals=slice_intervals,
        exclude_intervals=exclude_intervals,
        location=location,
    )
    n_after = int(filtered.sizes["time"])
    _diagnose_signal_stage(filtered, label="after filtering")

    averaged = filtered.mean(dim="time", skipna=True)
    _diagnose_signal_stage(averaged, label="after averaging")

    return averaged, {
        "mode": "filtered_mean",
        "n_before": n_before,
        "n_after": n_after,
    }


def _mean_error(
    value: Any,
    *,
    slice_intervals: Sequence[tuple[str, str]],
    exclude_intervals: Sequence[tuple[str, str]],
    location: str,
) -> Any:
    """Filter and collapse uncertainty using the existing ATLAS mean-error rule."""

    if not isinstance(value, xr.DataArray) or "time" not in value.dims:
        return value

    if value.sizes["time"] == 1:
        return value.squeeze("time", drop=True)

    filtered = _filter_time(
        value,
        slice_intervals=slice_intervals,
        exclude_intervals=exclude_intervals,
        location=location,
    )

    n = filtered.notnull().sum(dim="time")
    mean_error = filtered.mean(dim="time", skipna=True) / np.sqrt(n)
    return mean_error.where(n > 0)


def _process_group_collection(
    groups: Mapping[str, Mapping[str, Any]],
    *,
    group_kind: str,
    slice_intervals: Sequence[tuple[str, str]],
    exclude_intervals: Sequence[tuple[str, str]],
) -> dict[str, dict[str, Any]]:
    """Process channel or pair groups while retaining bundle metadata."""

    out: dict[str, dict[str, Any]] = {}

    for group_id, group in groups.items():
        group_out = _copy_nested(group)

        for dataset_id, dataset in group_out.get("entries", {}).items():
            base = f"[{group_kind}:{group_id}] entry {dataset_id!r}"

            dataset["signal"], time_status = _mean_signal(
                dataset.get("signal"),
                slice_intervals=slice_intervals,
                exclude_intervals=exclude_intervals,
                location=f"{base} signal",
            )
            dataset["error"] = _mean_error(
                dataset.get("error"),
                slice_intervals=slice_intervals,
                exclude_intervals=exclude_intervals,
                location=f"{base} error",
            )

            mode = time_status["mode"]
            if mode == "singleton":
                print(f"    - {dataset_id}: singleton time dimension -> dropped time; slicing skipped")
            elif mode == "filtered_mean":
                print(
                    f"    - {dataset_id}: {time_status['n_before']} -> "
                    f"{time_status['n_after']} time samples after filtering -> averaged"
                )
            else:
                print(f"    - {dataset_id}: no time dimension -> left unchanged")

        out[group_id] = group_out

    return out


def filter_and_average_intercomparison_bundles(
    intercomparison_info: Mapping[str, Any],
    intercomparison_bundles: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply global temporal slice/exclude rules and collapse time.

    ``slice_measurement`` and ``exclude_measurement`` are already validated and
    converted by ``parse_intercomparison_ini`` to lists of ``(start, stop)``
    tuples. The same global windows are applied to every participating entry.

    A singleton ``time`` dimension is interpreted as an already averaged
    exported product: no slicing is attempted and ``time`` is simply dropped.
    """

    general = intercomparison_info["general"]
    slice_intervals = general.get("slice_measurement") or []
    exclude_intervals = general.get("exclude_measurement") or []

    print()
    print("-----------------------------------------------")
    print("Temporal filtering and averaging")
    print("-----------------------------------------------")
    print(f"Slice intervals: {slice_intervals if slice_intervals else 'none'}")
    print(f"Exclude intervals: {exclude_intervals if exclude_intervals else 'none'}")

    channel_groups = intercomparison_bundles.get("channel_groups", {})
    pair_groups = intercomparison_bundles.get("pair_groups", {})
    if channel_groups:
        print("Channel groups:")
    else:
        print("Channel groups: none")

    out = {
        "datasets": _copy_nested(intercomparison_bundles.get("datasets", {})),
        "channel_groups": _process_group_collection(
            channel_groups,
            group_kind="channel_group",
            slice_intervals=slice_intervals,
            exclude_intervals=exclude_intervals,
        ),
        "pair_groups": {},
    }

    if pair_groups:
        print("Pair groups:")
        out["pair_groups"] = _process_group_collection(
            pair_groups,
            group_kind="pair_group",
            slice_intervals=slice_intervals,
            exclude_intervals=exclude_intervals,
        )
    else:
        print("Pair groups: none")

    print("Temporal filtering and averaging complete.")
    print()
    return out
