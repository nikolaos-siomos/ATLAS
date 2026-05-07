#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 17:55:19 2026

@author: nikos
"""
import numpy as np
import copy
import xarray as xr
import pandas as pd
from helper_functions.printouts import endpoint
from typing import Any, Dict
from utils.time_conversions import iso_to_datetimes
from utils.error_classes import CustomWarning

def hhmm_to_datetime(hhmm: str, base: pd.Timestamp) -> pd.Timestamp:
    """Convert hhmm string to datetime on the date of base timestamp."""
    return pd.Timestamp(year=base.year,
                        month=base.month,
                        day=base.day,
                        hour=int(hhmm[:2]),
                        minute=int(hhmm[2:]))

def make_interval(start_str: str, stop_str: str, base: pd.Timestamp):
    start_dt = hhmm_to_datetime(start_str, base)
    stop_dt = hhmm_to_datetime(stop_str, base)

    # If stop is "before" start → assume it's the next day
    if stop_dt <= start_dt:
        stop_dt += pd.Timedelta(days=1)

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

def compute_slice_and_exclude(
        processing_info: Dict[str, Any],input_data: Dict[str, Dict[str, Any]]) -> \
    Dict[str, Dict[str, Any]]:
    
    output_data = copy.deepcopy(input_data)
    # output_data.setdefault("time_mask", {})

    profiles = output_data["profile"]
    time_info = output_data["time_info"]

    qa_tests = list(profiles.keys())

    # Start with full-True masks for all measurements
    mask_time = {
        key: xr.ones_like(profiles[key]["time"], dtype=bool)
        for key in qa_tests
    }

    slicer = processing_info['caller_info']["slice_measurement"]
    if slicer is not None:
        slice_meas = [slicer[i] for i in range(0, len(slicer), 3)]
        slice_start_time = [slicer[i] for i in range(1, len(slicer), 3)]
        slice_stop_time = [slicer[i] for i in range(2, len(slicer), 3)]

        # If slicing is specified, start from all-False for mentioned keys
        for key in slice_meas:
            if key in qa_tests:
                mask_time[key] = xr.zeros_like(profiles[key]["time"], dtype=bool)

        for i, key in enumerate(slice_meas):
            if key in qa_tests:
                start_times = iso_to_datetimes(
                    time_info[key].sel({"parameters": "start_time"}).values
                )

                start_dt, stop_dt = make_interval(
                    start_str=slice_start_time[i],
                    stop_str=slice_stop_time[i],
                    base=start_times[0]
                )

                time = profiles[key]["time"]
                mask_interval = (time >= start_dt) & (time <= stop_dt)
                mask_time[key] = mask_time[key] | mask_interval

    slicer = processing_info['caller_info']["exclude_measurement"]
    if slicer is not None:
        slice_meas = [slicer[i] for i in range(0, len(slicer), 3)]
        slice_start_time = [slicer[i] for i in range(1, len(slicer), 3)]
        slice_stop_time = [slicer[i] for i in range(2, len(slicer), 3)]

        for i, key in enumerate(slice_meas):
            if key in qa_tests:
                start_times = iso_to_datetimes(
                    time_info[key].sel({"parameters": "start_time"}).values
                )

                start_dt, stop_dt = make_interval(
                    start_str=slice_start_time[i],
                    stop_str=slice_stop_time[i],
                    base=start_times[0]
                )

                time = profiles[key]["time"]
                mask_interval = (time < start_dt) | (time > stop_dt)
                mask_time[key] = mask_time[key] & mask_interval

    output_data = apply_time_mask(qa_tests, mask_time, output_data)
    
    if not output_data["profile"]:
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

    output_data = apply_time_mask(qa_tests, keep_time_mask, output_data)
    
    if not output_data["profile"]:
        endpoint(4)

    return output_data