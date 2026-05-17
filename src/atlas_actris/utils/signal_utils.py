#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 20:13:56 2025

@author: nikos
"""

import numpy as np
import xarray as xr
import pandas as pd
from scipy.signal import savgol_coeffs

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

def temporal_averaging_error(
        sig_err: xr.DataArray, averaging_rate: str, 
        averaging_threshold: float) -> Tuple[xr.DataArray, xr.DataArray]:
    
    delta_t_min = np.min(sig_err.time[1:].values-sig_err.time[:-1].values)
    
    expected_profiles = expected_profiles_per_window(averaging_rate, delta_t_min)
    
    if expected_profiles < 1:
        CustomWarning(f"The provided averaging_rate ({averaging_rate}) is smaller than the temporal resolution. Averaging is not possible, the raw temporal resolution will be used")
        sig_avg_err = sig_err
    
    else:    
        
        origin = pd.Timestamp(sig_err.time.values[0])

        actual_profiles = sig_err.resample({"time" : averaging_rate},
                                           label = "left", 
                                           origin = origin).count()
        
        mask_incomplete = actual_profiles / expected_profiles < averaging_threshold
        
        sig_mean_err = sig_err.resample({"time" : averaging_rate},
                                       label = "left", 
                                       origin = origin).mean()

        sig_avg_err = sig_mean_err / np.sqrt(actual_profiles)
        
        sig_avg_err = sig_avg_err.where(~mask_incomplete, np.nan)
        
    return(sig_avg_err, mask_incomplete)

def rolling_noise(
    sig: xr.DataArray,
    smooth_window: int = 5,
    noise_window: int = 31,
    dim: str = "bins",
):
    # 1. Local smooth signal estimate
    sig_sm = (
        sig
        .rolling({dim: smooth_window}, center=True, min_periods=smooth_window)
        .mean()
    )

    # 2. Residual = high-frequency component
    residual = sig - sig_sm

    # 3. Local rolling noise estimate
    noise = (
        residual
        .rolling({dim: noise_window}, center=True, min_periods=noise_window)
        .std()
    )

    return noise

def rolling_noise_savgol(
    sig: xr.DataArray,
    polyorder: int = 3,
    noise_window: int = 31,
    dim: str = "bins",
):
    """
    Estimate local noise using:
      1. Savitzky-Golay low-pass smoothing along `dim`
      2. residual = signal - smoothed signal
      3. rolling std of residual along `dim`

    Good when real features, such as aerosol layers, should be less distorted
    than with a simple rolling mean.
    """

    if noise_window % 2 == 0:
        raise ValueError("noise_window should usually be odd when center=True.")
    if polyorder >= noise_window:
        raise ValueError("polyorder must be smaller than smooth_window.")

    win_dim = "_savgol_window"

    coeffs = xr.DataArray(
        savgol_coeffs(
            window_length=noise_window,
            polyorder=polyorder,
            use="dot",
        ),
        dims=[win_dim],
    )

    rolled = (
        sig
        .rolling({dim: noise_window}, center=True, min_periods=noise_window)
        .construct(win_dim)
    )

    sig_sm = (rolled * coeffs).sum(win_dim, skipna=False)

    residual = sig - sig_sm

    noise = (
        residual
        .rolling({dim: noise_window}, center=True, min_periods=noise_window)
        .std()
    )

    return noise
    