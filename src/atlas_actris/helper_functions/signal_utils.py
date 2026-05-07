#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 20:13:56 2025

@author: nikos
"""

import numpy as np
import xarray as xr
import pandas as pd
from typing import Tuple
from utils.error_classes import CustomWarning


def expected_profiles_per_window(freq: str,
                                 sample_dt: np.timedelta64,
                                 include_end: bool = False) -> int:
    """
    Return the expected number of profiles within one resampling window.

    freq        : pandas offset alias (e.g. '3min', '1H', '90s', '1H30min')
    sample_dt   : np.timedelta64 (e.g. np.timedelta64(15_000_000_000, 'ns'))
    include_end : if True, count includes the right edge when it lands exactly on a sample
                  (i.e., closed='both' / right-closed windows)

    Assumes windows are left-closed, right-open by default (include_end=False),
    matching xarray/pandas resample defaults.
    """
    # Convert both to integer nanoseconds to avoid unit mismatches
    W_ns = pd.to_timedelta(freq).value                  # window length in ns (int)
    dt_ns = pd.to_timedelta(sample_dt).value            # sample spacing in ns (int)

    if W_ns <= 0 or dt_ns <= 0:
        raise ValueError("freq and sample_dt must be positive durations.")

    count = W_ns // dt_ns  # floor(W/dt) for [left-closed, right-open) bins

    # If you want to include the right edge when exactly aligned:
    if include_end and (W_ns % dt_ns == 0):
        count += 1

    return int(count)


def temporal_averaging(sig: xr.DataArray, averaging_rate: str, 
                       averaging_threshold: float) -> Tuple[xr.DataArray, xr.DataArray]:
    
    delta_t_min = np.min(sig.time[1:].values-sig.time[:-1].values)
    
    expected_profiles = expected_profiles_per_window(averaging_rate, delta_t_min)
    
    if expected_profiles < 1:
        CustomWarning(f"The provided averaging_rate ({averaging_rate}) is smaller than the temporal resolution. Averaging is not possible, the raw temporal resolution will be used")
        sig_avg = sig
    
    else:    
        
        origin = pd.Timestamp(sig.time.values[0])

        actual_profiles = sig.resample({"time" : averaging_rate},
                                       label = "left", 
                                       origin = origin).count()
        
        mask_incomplete = actual_profiles / expected_profiles < averaging_threshold
        
        sig_avg = sig.resample({"time" : averaging_rate},
                               label = "left", 
                               origin = origin).mean()
    
        sig_avg = sig_avg.where(~mask_incomplete, np.nan)
        
    return(sig_avg, mask_incomplete)
    