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
from scipy.ndimage import uniform_filter1d

from typing import Tuple, Optional, Set, Dict, Any
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


def temporal_averaging(
    sig: xr.DataArray,
    averaging_period: str,
    averaging_threshold: float,
) -> Tuple[Optional[xr.DataArray], Optional[xr.DataArray]]:
    """Average along time, or return ``(None, None)`` when averaging is impossible.

    Returning ``None`` allows the calling processing stage to leave its already
    initialized output arrays unchanged when ``averaging_period`` is shorter than
    the native measurement resolution.
    """

    if sig.sizes.get("time", 0) < 2:
        CustomWarning(
            "Temporal averaging was skipped because fewer than two time "
            "samples are available."
        )
        return None, None

    delta_t_min = np.min(sig.time[1:].values - sig.time[:-1].values)
    expected_profiles = expected_profiles_per_window(averaging_period, delta_t_min)

    if expected_profiles < 1:
        CustomWarning(
            f"The provided averaging_period ({averaging_period}) is smaller than "
            "the temporal resolution. Averaging was skipped."
        )
        return None, None

    origin = pd.Timestamp(sig.time.values[0])

    actual_profiles = sig.resample(
        {"time": averaging_period},
        label="left",
        origin=origin,
    ).count()

    mask_incomplete = (
        actual_profiles / expected_profiles < averaging_threshold
    )

    sig_avg = sig.resample(
        {"time": averaging_period},
        label="left",
        origin=origin,
    ).mean()

    sig_avg = sig_avg.where(~mask_incomplete, np.nan)

    return sig_avg, mask_incomplete


def temporal_averaging_error(
    sig_err: xr.DataArray,
    averaging_period: str,
    averaging_threshold: float,
) -> Tuple[Optional[xr.DataArray], Optional[xr.DataArray]]:
    """Average temporal errors, or return ``(None, None)`` when skipped."""

    if sig_err.sizes.get("time", 0) < 2:
        CustomWarning(
            "Temporal error averaging was skipped because fewer than two "
            "time samples are available."
        )
        return None, None

    delta_t_min = np.min(
        sig_err.time[1:].values - sig_err.time[:-1].values
    )
    expected_profiles = expected_profiles_per_window(
        averaging_period,
        delta_t_min,
    )

    if expected_profiles < 1:
        CustomWarning(
            f"The provided averaging_period ({averaging_period}) is smaller than "
            "the temporal resolution. Error averaging was skipped."
        )
        return None, None

    origin = pd.Timestamp(sig_err.time.values[0])

    actual_profiles = sig_err.resample(
        {"time": averaging_period},
        label="left",
        origin=origin,
    ).count()

    mask_incomplete = (
        actual_profiles / expected_profiles < averaging_threshold
    )

    sig_mean_err = sig_err.resample(
        {"time": averaging_period},
        label="left",
        origin=origin,
    ).mean()

    sig_avg_err = sig_mean_err / np.sqrt(actual_profiles)
    sig_avg_err = sig_avg_err.where(~mask_incomplete, np.nan)

    return sig_avg_err, mask_incomplete

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

def fast_rolling_mean(
    da: xr.DataArray,
    window: int,
    dim: str = "bins",
    min_periods: int = None,
) -> xr.DataArray:

    if min_periods is None:
        min_periods = window

    axis = da.get_axis_num(dim)

    values = da.values.astype(float, copy=False)
    valid = np.isfinite(values)

    values_filled = np.where(valid, values, 0.0)

    summed = (
        uniform_filter1d(
            values_filled,
            size=window,
            axis=axis,
            mode="constant",
            cval=0.0,
        )
        * window
    )

    counts = (
        uniform_filter1d(
            valid.astype(float),
            size=window,
            axis=axis,
            mode="constant",
            cval=0.0,
        )
        * window
    )

    mean = summed / counts
    mean = np.where(counts >= min_periods, mean, np.nan)

    return xr.DataArray(
        mean.astype(da.dtype, copy=False),
        dims=da.dims,
        coords=da.coords,
        attrs=da.attrs,
        name=da.name,
    )

def fast_rolling_std(
    da: xr.DataArray,
    window: int,
    dim: str = "bins",
    min_periods: int = None,
) -> xr.DataArray:

    if min_periods is None:
        min_periods = window

    axis = da.get_axis_num(dim)

    values = da.values.astype(float, copy=False)
    valid = np.isfinite(values)

    values_filled = np.where(valid, values, 0.0)
    values2_filled = np.where(valid, values**2, 0.0)

    summed = (
        uniform_filter1d(
            values_filled,
            size=window,
            axis=axis,
            mode="constant",
            cval=0.0,
        )
        * window
    )

    summed2 = (
        uniform_filter1d(
            values2_filled,
            size=window,
            axis=axis,
            mode="constant",
            cval=0.0,
        )
        * window
    )

    counts = (
        uniform_filter1d(
            valid.astype(float),
            size=window,
            axis=axis,
            mode="constant",
            cval=0.0,
        )
        * window
    )

    mean = summed / counts
    mean2 = summed2 / counts

    var = mean2 - mean**2
    var = np.maximum(var, 0.0)

    std = np.sqrt(var)
    std = np.where(counts >= min_periods, std, np.nan)

    return xr.DataArray(
        std.astype(da.dtype, copy=False),
        dims=da.dims,
        coords=da.coords,
        attrs=da.attrs,
        name=da.name,
    )

def fast_rolling_noise(
    sig: xr.DataArray,
    smooth_window: int = 5,
    noise_window: int = 31,
    dim: str = "bins",
):
    # 1. Local smooth signal estimate
    sig_sm = fast_rolling_mean(
        da=sig,
        window=smooth_window,
        dim=dim,
        min_periods=smooth_window,
    )

    # 2. Residual = high-frequency component
    residual = sig - sig_sm

    # 3. Local rolling noise estimate
    noise = fast_rolling_std(
        da=residual,
        window=noise_window,
        dim=dim,
        min_periods=noise_window,
    )

    return noise

def fast_rolling_mean_range(
    da: xr.DataArray,
    ranges: xr.DataArray,
    window: int = 1000,
    smooth_above: float = 2000.0,
    dim: str = "bins",
) -> xr.DataArray:
    """
    Smooth da along bins and replace values above smooth_above range
    with smoothed values.

    The original dimensions are preserved, including singleton time.
    """

    axis = da.get_axis_num(dim)

    values = da.values.astype(float, copy=False)
    valid = np.isfinite(values)

    values_filled = np.where(valid, values, 0.0)

    summed = (
        uniform_filter1d(
            values_filled,
            size=window,
            axis=axis,
            mode="constant",
            cval=0.0,
        )
        * window
    )

    counts = (
        uniform_filter1d(
            valid.astype(float),
            size=window,
            axis=axis,
            mode="constant",
            cval=0.0,
        )
        * window
    )

    smooth_values = summed / counts
    smooth_values = np.where(counts > 0, smooth_values, np.nan)

    smoothed = xr.DataArray(
        smooth_values.astype(da.dtype, copy=False),
        dims=da.dims,
        coords=da.coords,
        attrs=da.attrs,
        name=da.name,
    )

    out = xr.where(ranges > smooth_above, smoothed, da)

    # Important: restore original dimension order and keep singleton time.
    out = out.transpose(*da.dims)

    return out
    
def _true_bin_grid(
    zero_bin: xr.DataArray,
    n_old_bins: int,
    resolution: xr.DataArray,
    zenith_angle_rad: float,
    max_height_agl: Optional[float] = None,
    ) -> np.ndarray:
    """
    Build common true-bin left-edge grid.

    Convention:
        true_bin = file_bin + zero_bin[channel]

    The returned grid:
        - starts at true bin 0
        - covers the union of all channel-shifted bins
        - is optionally limited by max_height_agl
    """

    if "channels" in zero_bin.dims and "channel" not in zero_bin.dims:
        zero_bin = zero_bin.rename({"channels": "channel"})

    if "channels" in resolution.dims and "channel" not in resolution.dims:
        resolution = resolution.rename({"channels": "channel"})

    zero_bin = zero_bin.astype("float64")
    resolution = resolution.astype("float64")

    first_valid = zero_bin
    last_valid = zero_bin + n_old_bins

    target_min = max(0, int(np.floor(float(first_valid.min().values))))
    target_max = int(np.ceil(float(last_valid.max().values)) - 1.0)

    if max_height_agl is not None:
        max_bin = (
            max_height_agl
            / (resolution * np.cos(zenith_angle_rad))
            - 0.5
        )

        # Use max across channels because each channel has a different resolution.
        target_max = min(
            target_max,
            int(np.floor(float(max_bin.max().values))),
        )

    if target_max < target_min:
        raise ValueError(
            "No valid bins remain after applying zero-bin and height limits."
        )

    return np.arange(target_min, target_max + 1, dtype="float64")

def _rebin_to_true_bins(
    da: xr.DataArray,
    zero_bin: xr.DataArray,
    target_bins: np.ndarray,
    output_dtype: str = "float32",
) -> xr.DataArray:
    """
    Lazily rebin a DataArray from file bins to true bins.

    Fast path
    ---------
    If all channel shifts are integers, use simple per-channel positional
    slicing plus padding. This avoids the expensive interpolation graph.

    Fractional path
    ---------------
    If any shift has a fractional component, use the conservative
    fractional-overlap implementation.
    """

    if "bins" not in da.dims or "channel" not in da.dims:
        return da

    if "channels" in zero_bin.dims and "channel" not in zero_bin.dims:
        zero_bin = zero_bin.rename({"channels": "channel"})

    zero_bin = zero_bin.sel(channel=da.channel)
    shifts = np.asarray(zero_bin.values, dtype="float64")
    shift_int = np.rint(shifts).astype("int64")
    is_integer_shift = np.all(np.isclose(shifts, shift_int, atol=1e-10, rtol=0.0))

    n_old_bins = da.sizes["bins"]
    target_bins = np.asarray(target_bins, dtype="float64")

    if is_integer_shift:
        out_channels = []
        target_start = int(target_bins[0])
        target_stop = int(target_bins[-1]) + 1

        for ch, shift in zip(da.channel.values, shift_int):
            # true_bin = file_bin + shift
            old_start = target_start - int(shift)
            old_stop = target_stop - int(shift)

            src_start = max(old_start, 0)
            src_stop = min(old_stop, n_old_bins)

            left_pad = max(0, -old_start)
            right_pad = max(0, old_stop - n_old_bins)

            da_ch = da.sel(channel=ch, drop=False)
            sliced = da_ch.isel(bins=slice(src_start, src_stop))

            pieces = []
            if left_pad:
                left = xr.full_like(
                    da_ch.isel(bins=slice(0, left_pad)),
                    False if np.issubdtype(da.dtype, np.bool_) else np.nan,
                )
                pieces.append(left)

            pieces.append(sliced)

            if right_pad:
                right = xr.full_like(
                    da_ch.isel(bins=slice(0, right_pad)),
                    False if np.issubdtype(da.dtype, np.bool_) else np.nan,
                )
                pieces.append(right)

            out_ch = xr.concat(pieces, dim="bins") if len(pieces) > 1 else pieces[0]
            out_ch = out_ch.assign_coords(bins=target_bins)

            if not np.issubdtype(da.dtype, np.bool_):
                out_ch = out_ch.astype(output_dtype)

            out_channels.append(out_ch)

        out = xr.concat(out_channels, dim="channel")
        out = out.assign_coords(channel=da.channel)
        return out.transpose(*da.dims)

    # Fractional-overlap fallback.
    out_channels = []

    for ch in da.channel.values:
        zb = float(zero_bin.sel(channel=ch).values)

        shift_int_ch = int(np.floor(zb))
        frac = float(zb - shift_int_ch)

        weight_main = 1.0 - frac
        weight_prev = frac

        old_i_main = xr.DataArray(
            target_bins - shift_int_ch,
            dims=["bins_new"],
            coords={"bins_new": target_bins},
        )

        old_i_prev = old_i_main - 1

        valid_main = (
            (old_i_main >= 0)
            & (old_i_main < n_old_bins)
            & (weight_main > 0.0)
        )

        valid_prev = (
            (old_i_prev >= 0)
            & (old_i_prev < n_old_bins)
            & (weight_prev > 0.0)
        )

        valid_out = valid_main | valid_prev

        da_ch = da.sel(channel=ch)
        work = da_ch.rename({"bins": "bins_old"})
        work = work.assign_coords(bins_old=np.arange(n_old_bins))

        old_i_main_clip = old_i_main.clip(
            min=0,
            max=n_old_bins - 1,
        ).astype("int64")

        old_i_prev_clip = old_i_prev.clip(
            min=0,
            max=n_old_bins - 1,
        ).astype("int64")

        main = work.isel(bins_old=old_i_main_clip)
        prev = work.isel(bins_old=old_i_prev_clip)

        if np.issubdtype(da.dtype, np.bool_):
            out_ch = main.where(valid_main, False) | prev.where(valid_prev, False)
            out_ch = out_ch.where(valid_out, False)
        else:
            out_ch = (
                main.where(valid_main, 0.0) * weight_main
                + prev.where(valid_prev, 0.0) * weight_prev
            )
            out_ch = out_ch.where(valid_out, np.nan)
            out_ch = out_ch.astype(output_dtype)

        out_ch = out_ch.rename({"bins_new": "bins"})
        out_ch = out_ch.assign_coords(bins=target_bins)
        out_channels.append(out_ch)

    out = xr.concat(out_channels, dim=da.channel)
    out = out.assign_coords(channel=da.channel)

    return out.transpose(*da.dims)


def _rebin_and_trim_all_binned_arrays(
    output_data: Dict[str, Any],
    key: str,
    zero_bin: xr.DataArray,
    target_bins: np.ndarray,
    mask_bins: xr.DataArray,
    skip_stores: set[str] = {"range", "height_agl", "height_asl"},
) -> None:
    """
    Rebin and mask all channel/bin arrays for one QA key.

    Repeated references to the same DataArray are transformed only once.
    ``drop=True`` is intentionally avoided: ``target_bins`` already defines
    the shared retained grid, while the channel-specific height mask only
    needs to set invalid channel/bin cells to NaN (or False).
    """

    transformed: Dict[int, xr.DataArray] = {}

    for store_name, store in output_data.items():

        if store_name in skip_stores:
            continue
        if not isinstance(store, dict) or key not in store:
            continue

        arr = store[key]

        if not isinstance(arr, xr.DataArray):
            continue
        if "channel" not in arr.dims or "bins" not in arr.dims:
            continue

        arr_id = id(arr)
        if arr_id not in transformed:
            rebinned = _rebin_to_true_bins(
                da=arr,
                zero_bin=zero_bin,
                target_bins=target_bins,
                output_dtype="float32",
            )

            channel_mask = mask_bins
            if np.issubdtype(rebinned.dtype, np.bool_):
                rebinned = rebinned.where(channel_mask, False)
            else:
                rebinned = rebinned.where(channel_mask)

            transformed[arr_id] = rebinned.reset_coords(drop=True)

        store[key] = transformed[arr_id]
