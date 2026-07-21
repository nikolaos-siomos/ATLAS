#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 17:55:19 2026

@author: nikos
"""
import numpy as np
import copy
import re
import xarray as xr
import pandas as pd

from typing import Any, Dict, Optional
from utils.printouts import endpoint
from utils.error_classes import CustomWarning
from utils.dataarray_utils import shallow_copy
from utils.time_conversions import iso_to_datetimes

from utils.printouts import print_header, print_subsection, print_entry

from processor.definitions import (
    assign_drk, 
    profile_instances, 
    )

def hhmm_to_datetime(hhmm: str, base: pd.Timestamp) -> pd.Timestamp:
    """Convert HHMM string to datetime on the date of base timestamp."""

    value = str(hhmm).strip()

    if not value.isdigit() or len(value) != 4:
        raise ValueError(
            "HHMM times must contain exactly 4 digits, for example '2130'. "
            f"Got {hhmm!r}."
        )

    return pd.Timestamp(
        year=base.year,
        month=base.month,
        day=base.day,
        hour=int(value[:2]),
        minute=int(value[2:]),
    )


def _parse_slice_time(value: str, base: pd.Timestamp) -> tuple[pd.Timestamp, bool]:
    """Parse a slice/exclude time string.

    Accepted formats are:
    - HHMM              -> time on the date of ``base``
    - YYYYMMDD          -> absolute date at 00:00:00
    - YYYYMMDD_HH       -> absolute date and hour
    - YYYYMMDD_HHMM     -> absolute date, hour, and minute
    - YYYYMMDD_HHMMSS   -> absolute date, hour, minute, and second

    Returns
    -------
    timestamp, is_relative_hhmm
        ``is_relative_hhmm`` is True only for the legacy HHMM format.
        It is used by ``make_interval`` to preserve the old overnight behavior.
    """

    text = str(value).strip()

    if text.isdigit() and len(text) == 4:
        return hhmm_to_datetime(text, base), True

    match = re.fullmatch(r"(\d{8})(?:_(\d{2}|\d{4}|\d{6}))?", text)

    if match is None:
        raise ValueError(
            "Invalid slice/exclude time format. Accepted formats are: "
            "HHMM, YYYYMMDD, YYYYMMDD_HH, YYYYMMDD_HHMM, "
            f"YYYYMMDD_HHMMSS. Got {value!r}."
        )

    date_part, time_part = match.groups()
    year = int(date_part[:4])
    month = int(date_part[4:6])
    day = int(date_part[6:8])

    hour = 0
    minute = 0
    second = 0

    if time_part is not None:
        hour = int(time_part[:2])

        if len(time_part) >= 4:
            minute = int(time_part[2:4])

        if len(time_part) == 6:
            second = int(time_part[4:6])

    return pd.Timestamp(
        year=year,
        month=month,
        day=day,
        hour=hour,
        minute=minute,
        second=second,
    ), False


def make_interval(start_str: str, stop_str: str, base: pd.Timestamp):
    """Create a start/stop interval from supported slice/exclude strings.

    The legacy HHMM format keeps the old behavior: if the stop time is not
    after the start time, the stop time is assumed to be on the next day.
    Absolute YYYYMMDD-based formats are not auto-shifted; use the next date
    explicitly when the interval crosses midnight.
    """

    start_dt, start_is_hhmm = _parse_slice_time(start_str, base)
    stop_dt, stop_is_hhmm = _parse_slice_time(stop_str, base)

    if stop_dt <= start_dt:
        if start_is_hhmm and stop_is_hhmm:
            stop_dt += pd.Timedelta(days=1)
        else:
            raise ValueError(
                "The slice/exclude stop time must be after the start time "
                "for YYYYMMDD-based formats. For overnight absolute "
                "intervals, write the stop time with the next date. "
                f"Got start={start_str!r}, stop={stop_str!r}."
            )

    return start_dt, stop_dt

def apply_time_mask(qa_tests, mask_time, output_data):
    
    del_keys = []
    for key in qa_tests:
        if not mask_time[key].any().item():
            del_keys.append(key)

    for key in qa_tests:
        if key in del_keys:
            CustomWarning(
                f"{key} measurement will not be processed because masking "
                "removes all profiles"
            )
            del output_data["profile"][key]
            if key in output_data["time_mask"]:
                del output_data["time_mask"][key]
            del output_data["time_info"][key]
            del output_data["shots"][key]

        else:
            output_data["profile"][key] = \
                output_data["profile"][key].where(mask_time[key]).copy()
            
            output_data["time_mask"][key] = mask_time[key]

            output_data["time_info"][key] = \
                output_data["time_info"][key].where(mask_time[key]).copy()
            
            output_data["shots"][key] = \
                output_data["shots"][key].where(mask_time[key]).copy()
    
    return output_data

def select_by_time_mask(qa_tests, mask_time, output_data):
    """
    Apply a 1D time mask to all non-empty dictionaries inside output_data.

    True  = keep
    False = drop

    RAM-safe behavior:
    - compute only the small 1D mask once per measurement
    - slice large xarray objects with isel()
    - do not use where(..., drop=True)
    - do not create large NaN-filled arrays
    - do not add extra output_data sections
    """

    output_data.setdefault("time_mask", {})

    keep_indices = {}
    del_keys = []

    # Compute only the small 1D masks.
    for key in qa_tests:
        mask = mask_time[key]

        if hasattr(mask.data, "compute"):
            mask_values = mask.compute().values
        else:
            mask_values = mask.values

        keep_idx = np.flatnonzero(mask_values)

        if keep_idx.size == 0:
            del_keys.append(key)
        else:
            keep_indices[key] = keep_idx

    # Remove fully masked measurements from all non-empty output sections.
    for key in del_keys:
        CustomWarning(
            f"{key} measurement will not be processed because masking "
            "removes all profiles"
        )

        for _, subdict in output_data.items():
            if isinstance(subdict, dict) and subdict and key in subdict:
                del subdict[key]

    # Lazily slice every xarray object that has a time dimension.
    for key in qa_tests:
        if key in del_keys:
            continue

        keep_idx = keep_indices[key]

        for _, subdict in output_data.items():

            if not isinstance(subdict, dict):
                continue

            if not subdict:
                continue

            if key not in subdict:
                continue

            obj = subdict[key]

            if isinstance(obj, (xr.DataArray, xr.Dataset)) and "time" in obj.dims:
                subdict[key] = obj.isel(time=keep_idx)

        output_data["time_mask"][key] = mask_time[key].isel(time=keep_idx)

    return output_data

def _time_mask_from_time_info(time_info_key: xr.DataArray) -> xr.DataArray:
    """
    Create a full-True boolean mask using time_info as reference.
    """

    if "time" not in time_info_key.dims:
        raise ValueError(
            "time_info entry must contain a 'time' dimension in order "
            "to build the slice/exclude mask."
        )

    return xr.ones_like(time_info_key["time"], dtype=bool)


def _start_time_array_from_time_info(time_info_key: xr.DataArray) -> xr.DataArray:
    """
    Return start_time values from time_info as datetime values with the same
    time coordinates as time_info.
    """

    start_time = time_info_key.sel({"parameters": "start_time"})

    start_time_values = iso_to_datetimes(start_time.values)

    return xr.DataArray(
        start_time_values,
        dims=start_time.dims,
        coords=start_time.coords,
        name="start_time",
    )


def _resolve_loading_map_key(
    key: str,
    available_keys,
    loading_map,
    context: Optional[str] = None,
    warn: bool = True,
):
    """Resolve a QA-test key or loading-map alias to an available data key.

    Examples
    --------
    If ``key`` is ``"drk_ray"`` and ``loading_map["drk_ray"] == "drk"``,
    this returns ``"drk"`` when ``"drk"`` exists in ``available_keys``.
    Existing data keys are returned unchanged. Unresolvable keys return None.
    """

    if key in available_keys:
        return key

    if loading_map is None:
        loading_map = {}

    context_prefix = f"{context}: " if context else ""

    if key not in loading_map:
        if warn:
            CustomWarning(
                f"{context_prefix}key {key!r} was not found in available "
                "measurements or loading-map aliases. Entry skipped."
            )
        return None

    mapped_key = loading_map[key]

    if mapped_key in available_keys:
        return mapped_key

    if warn:
        CustomWarning(
            f"{context_prefix}key {key!r} maps to {mapped_key!r}, but the "
            "mapped key is not available. Entry skipped."
        )

    return None


def compute_slice_and_exclude(
        processing_info: Dict[str, Any],
        input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = copy.deepcopy(input_data)
    output_data.setdefault("time_mask", {})

    time_info = output_data["time_info"]

    qa_tests = list(time_info.keys())
    loading_map = processing_info["caller_info"].get("loading_map", {})

    # Start with full-True masks for all measurements.
    # The mask is based on the time dimension/coordinate of time_info.
    mask_time = {
        key: xr.ones_like(time_info[key]["time"], dtype=bool)
        for key in qa_tests
    }

    slicer = processing_info["caller_info"]["slice_measurement"]

    if slicer is not None:
        slice_meas = [slicer[i] for i in range(0, len(slicer), 3)]
        slice_start_time = [slicer[i] for i in range(1, len(slicer), 3)]
        slice_stop_time = [slicer[i] for i in range(2, len(slicer), 3)]

        # If slicing is specified, start from all-False for mentioned keys.
        # Mentioned keys may be direct output_data keys or loading-map aliases.
        for key in slice_meas:
            resolved_key = _resolve_loading_map_key(
                key, qa_tests, loading_map, context="slice_measurement"
            )

            if resolved_key is not None:
                mask_time[resolved_key] = xr.zeros_like(
                    time_info[resolved_key]["time"],
                    dtype=bool,
                )

        for i, key in enumerate(slice_meas):
            resolved_key = _resolve_loading_map_key(
                key, qa_tests, loading_map, context="slice_measurement"
            )

            if resolved_key is None:
                continue

            time = time_info[resolved_key]["time"]

            start_times = iso_to_datetimes(
                time_info[resolved_key].sel({"parameters": "start_time"}).values
            )

            start_dt, stop_dt = make_interval(
                start_str=slice_start_time[i],
                stop_str=slice_stop_time[i],
                base=start_times[0],
            )

            mask_interval = (time >= start_dt) & (time <= stop_dt)

            mask_time[resolved_key] = mask_time[resolved_key] | mask_interval

    slicer = processing_info["caller_info"]["exclude_measurement"]

    if slicer is not None:
        slice_meas = [slicer[i] for i in range(0, len(slicer), 3)]
        slice_start_time = [slicer[i] for i in range(1, len(slicer), 3)]
        slice_stop_time = [slicer[i] for i in range(2, len(slicer), 3)]

        for i, key in enumerate(slice_meas):
            resolved_key = _resolve_loading_map_key(
                key, qa_tests, loading_map, context="exclude_measurement"
            )

            if resolved_key is None:
                continue

            time = time_info[resolved_key]["time"]

            start_times = iso_to_datetimes(
                time_info[resolved_key].sel({"parameters": "start_time"}).values
            )

            start_dt, stop_dt = make_interval(
                start_str=slice_start_time[i],
                stop_str=slice_stop_time[i],
                base=start_times[0],
            )

            mask_interval = (time < start_dt) | (time > stop_dt)

            mask_time[resolved_key] = mask_time[resolved_key] & mask_interval

    output_data = select_by_time_mask(qa_tests, mask_time, output_data)

    if "profile" in output_data and not output_data["profile"]:
        endpoint(3)

    return output_data

def compute_screen_low_shots(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Replaces values in all bins with NaNs if the number of shots is
    lower than `shot_limit` of the maximum shots for that channel.

    Behavior:
    1. Prints a warning for channels that contain low-shot measurements
       across the time dimension.
    2. Creates a time mask that is True where any channel has low shots.
    3. Passes the inverse of that mask to `apply_time_mask`, assuming
       `apply_time_mask` expects True = keep.
    """

    low_shot_threshold = processing_info['caller_info']['low_shot_threshold']
    low_shot_threshold_percent = round(100.0 * low_shot_threshold)

    output_data = copy.deepcopy(input_data)
    output_data.setdefault("time_mask", {})

    profiles = output_data["profile"]
    time_info = output_data["time_info"]
    shots = output_data["shots"]

    qa_tests = list(profiles.keys())

    # True where any channel has low shots at that time
    bad_time_mask = {}

    for key in qa_tests:
        filename = time_info[key].sel({"parameters": "filename"})

        # True where shots are low for a given time/channel
        mask_low_shots = xr.where(
            (shots[key] / shots[key].max("time") < low_shot_threshold) &
            (shots[key].max("time") > 20),
            True,
            False
        ).compute()

        # True for channels that have at least one low-shot time
        any_low_shots_by_channel = mask_low_shots.any("time")

        # Get channel labels safely
        if "channel" in any_low_shots_by_channel.coords:
            channel_values = any_low_shots_by_channel.coords["channel"].values
        elif "channel" in any_low_shots_by_channel.dims:
            channel_values = np.arange(any_low_shots_by_channel.sizes["channel"])
        else:
            raise KeyError(
                f"'channel' not found in dims/coords for dataset '{key}'. "
                f"dims={any_low_shots_by_channel.dims}, "
                f"coords={list(any_low_shots_by_channel.coords)}"
            )

        low_channel_values = channel_values[any_low_shots_by_channel.values]

        if len(low_channel_values) > 0:
            print(f"{key} dataset")
            print(
                f"-- Warning: The following files have less shots than "
                f"{low_shot_threshold_percent}% of the maximum number of shots "
                f"encountered and will be screened out:"
            )

            for ch in low_channel_values:
                # Boolean mask over time for this channel only
                mask_time_ch = mask_low_shots.sel(channel=ch)

                # Times where this channel has low shots
                time_low_shots = mask_low_shots["time"].values[mask_time_ch.values]

                for t in time_low_shots:
                    print(f"    Filename: {filename.loc[t].values} | Channel: {ch}")

            print("")

        # True over time where ANY channel has low shots
        bad_time_mask[key] = mask_low_shots.any("channel")

    # If apply_time_mask expects True = keep, invert the bad-time mask
    keep_time_mask = {key: ~bad_time_mask[key] for key in qa_tests}

    output_data = select_by_time_mask(qa_tests, keep_time_mask, output_data)
    
    if not output_data["profile"]:
        endpoint(4)

    return output_data

def compute_asign_dark_blocks(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    loading_map = processing_info["caller_info"]["loading_map"]

    output_data = shallow_copy(input_data)

    profiles = output_data["profile"]

    for key, drk_key_alias in assign_drk.items():

        if not profiles:
            continue

        if key not in profiles:
            continue

        if drk_key_alias in profiles:
            drk_key = drk_key_alias
        elif drk_key_alias in loading_map:
            drk_key = loading_map[drk_key_alias]
        else:
            continue

        if drk_key not in profiles:
            continue

        drk = profiles[drk_key]
        sig = profiles[key]

        best_start_time, best_end_time, is_block = (
            select_dark_block_closest_to_signal(sig, drk)
        )

        if not is_block:
            continue

        if drk_key_alias in profiles:
            continue

        t_slice = slice(best_start_time, best_end_time)

        print_subsection(key)

        start_time = pd.Timestamp(best_start_time).strftime("%Y%m%d %H:%M:%S")
        end_time = pd.Timestamp(best_end_time).strftime("%Y%m%d %H:%M:%S")

        print_entry(
            f"Assigned part of {drk_key} to {drk_key_alias}:\n"
            f"  --start: {start_time}\n"
            f"  --end: {end_time}"
        )

        for db_key, db in output_data.items():

            if not isinstance(db, dict):
                continue

            if drk_key in db:
                src_key = drk_key
            elif key in db:
                src_key = key
            else:
                continue

            arr = db[src_key]

            if isinstance(arr, xr.DataArray) and "time" in arr.dims:
                arr = arr.sel({"time": t_slice})

            output_data[db_key][drk_key_alias] = arr

    print_entry("Dark measurement blocks successfully assigned!")

    return output_data

def select_dark_block_closest_to_signal(
    sig: xr.DataArray,
    drk: xr.DataArray,
    max_gap: str = "1h",
    time_dim: str = "time",
):
    """
    Select the dark time block closest to the mean signal time only if
    the dark measurement contains separated time blocks.

    Returns
    -------
    best_start_time : np.datetime64
        Start time of the selected dark interval.

    best_end_time : np.datetime64
        End time of the selected dark interval.

    is_block : bool
        False if the full/original dark measurement should be used.
        True if the dark measurement should be sliced to one selected time block.
    """

    if time_dim not in sig.dims:
        raise ValueError(f"sig has no '{time_dim}' dimension")

    if time_dim not in drk.dims:
        raise ValueError(f"drk has no '{time_dim}' dimension")

    if drk.sizes[time_dim] == 0:
        raise ValueError("drk has no time entries")

    # Keep these only if you want to guarantee correct ordering.
    # If you truly need the exact original object returned when no slicing happens,
    # remove these two lines.
    drk = drk.sortby(time_dim)
    sig = sig.sortby(time_dim)

    drk_time = drk[time_dim]

    # Single dark profile: no block selection possible.
    # Return False because no separated dark block was selected.
    if drk.sizes[time_dim] == 1:
        start_time = drk_time.values[0]
        end_time = drk_time.values[0]
        return start_time, end_time, False

    gaps = drk_time.diff(time_dim)
    max_gap_td = np.timedelta64(pd.Timedelta(max_gap).value, "ns")

    split_after = np.where(gaps.values > max_gap_td)[0]

    # No gap larger than max_gap.
    # This means drk is one continuous block, so return the full drk.
    if len(split_after) == 0:
        start_time = drk_time.values[0]
        end_time = drk_time.values[-1]
        return start_time, end_time, False

    # Here, separated blocks exist, so select only one block.
    starts = np.r_[0, split_after + 1]
    ends = np.r_[split_after + 1, drk.sizes[time_dim]]

    sig_mean_time = sig[time_dim].mean().values

    best_i = np.argmin([
        abs(
            drk_time.isel({time_dim: slice(start, end)}).mean().values
            - sig_mean_time
        )
        for start, end in zip(starts, ends)
    ])

    best_start = starts[best_i]
    best_end = ends[best_i]

    best_start_time = drk_time.values[best_start]
    best_end_time = drk_time.values[best_end - 1]

    return best_start_time, best_end_time, True