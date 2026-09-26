#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare and generate ATLAS intercomparison group plots.

This module plays the role that ``qa_test_rayleigh_fit.py`` plays for the
Rayleigh-fit plot: it extracts one comparison group, applies optional plotting
smoothing / local-STD estimation, determines plot limits, prepares relative
comparisons to the reference entry, creates text, and calls
``visualizer.plot_intercomparison``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

import numpy as np
import xarray as xr

from visualizer.make_text import (
    GenerateIntercomparisonText,
    IntercomparisonLibraries,
)
from visualizer.plot_utils import smoothing
from visualizer import plot_intercomparison


PHYSICAL_VERTICAL_SCALES = {"range", "height_agl", "height_asl"}


def _to_plot_units(values: xr.DataArray, vertical_scale: str) -> np.ndarray:
    arr = np.asarray(values.values, dtype=float)
    if vertical_scale in PHYSICAL_VERTICAL_SCALES:
        return 1.0e-3 * arr
    return arr


def _array_values(value: Optional[xr.DataArray]) -> Optional[np.ndarray]:
    if not isinstance(value, xr.DataArray):
        return None
    return np.asarray(value.values, dtype=float)


def _entry_label(
    intercomparison_info: Mapping[str, Any],
    entry_id: str,
    entry: Mapping[str, Any],
) -> str:
    if entry.get("entry_label"):
        return str(entry["entry_label"])

    dataset_id = str(entry.get("dataset_id", ""))
    dataset = intercomparison_info["datasets"].get(dataset_id, {})
    dataset_label = dataset.get("dataset_label") or dataset.get("system_label") or dataset_id
    product_id = entry.get("atlas_channel_id") or entry.get("atlas_pair_id")
    if product_id:
        return f"{dataset_label} - {product_id}"
    return dataset_label or entry_id


def _figure_title(
    intercomparison_info: Mapping[str, Any],
    group_id: str,
    group_kind: str,
    group: Mapping[str, Any],
) -> str:
    """Return a compact, parameter-oriented intercomparison figure title."""

    general = intercomparison_info.get("general", {})
    label = group.get("label") or group_id

    reference_dataset = intercomparison_info.get("reference_dataset", "reference")
    reference_entry = group.get("reference_entry")
    entries = group.get("entries", {})
    if reference_entry in entries:
        reference_dataset = entries[reference_entry].get(
            "dataset_id", reference_dataset
        )

    dataset = intercomparison_info.get("datasets", {}).get(reference_dataset, {})
    reference_label = (
        dataset.get("dataset_label")
        or dataset.get("system_label")
        or reference_dataset
    )

    if group_kind == "channel_group":
        plot_type = "Channel Intercomparison"
    elif group_kind == "pair_group":
        plot_type = "Pair Intercomparison"
    else:
        plot_type = "Intercomparison"

    vertical_scale = general.get("vertical_scale", "height_asl")
    if bool(general.get("plot_native_scale", False)):
        vertical_part = (
            f"vertical_scale={vertical_scale} | plot_native_scale=True"
        )
    else:
        method = general.get("vertical_method", "")
        vertical_part = (
            f"vertical_scale={vertical_scale} | vertical_method={method}"
        )
        if method == "vertical_binning":
            bin_width = group.get("vertical_bin_width")
            if bin_width is None:
                bin_width = general.get("vertical_bin_width")
            if bin_width is not None:
                vertical_part += f" | vertical_bin_width={bin_width}"

    return (
        f"ATLAS {plot_type} - {label}\n"
        f"Reference: {reference_label} | {vertical_part}"
    )

def _native_group_arrays(
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    X: Dict[str, np.ndarray] = {}
    Y: Dict[str, np.ndarray] = {}
    YE: Dict[str, Optional[np.ndarray]] = {}

    for entry_id, entry in group.get("entries", {}).items():
        signal = entry.get("signal")
        vertical = entry.get("metadata", {}).get(vertical_scale)
        if not isinstance(signal, xr.DataArray) or not isinstance(vertical, xr.DataArray):
            continue

        X[entry_id] = _to_plot_units(vertical, vertical_scale)
        Y[entry_id] = np.asarray(signal.values, dtype=float)
        YE[entry_id] = _array_values(entry.get("error"))

    return X, Y, YE


def _harmonized_group_arrays(
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    signals = group.get("signals")
    if not isinstance(signals, xr.DataArray) or "entry" not in signals.dims:
        raise ValueError(
            "Harmonized group array 'signals' is missing. Run vertical harmonization "
            "before plotting, or set plot_native_scale=True."
        )

    errors = group.get("errors")
    vertical = group.get("vertical_grid")
    if not isinstance(vertical, xr.DataArray):
        raise ValueError("Harmonized group vertical_grid is missing")

    x_common = _to_plot_units(vertical, vertical_scale)
    X: Dict[str, np.ndarray] = {}
    Y: Dict[str, np.ndarray] = {}
    YE: Dict[str, Optional[np.ndarray]] = {}

    entry_ids = [str(v) for v in signals["entry"].values]
    for entry_id in entry_ids:
        X[entry_id] = x_common.copy()
        Y[entry_id] = np.asarray(signals.sel(entry=entry_id).values, dtype=float)
        if isinstance(errors, xr.DataArray) and entry_id in [str(v) for v in errors["entry"].values]:
            YE[entry_id] = np.asarray(errors.sel(entry=entry_id).values, dtype=float)
        else:
            YE[entry_id] = None

    return X, Y, YE


def _reference_molecular(
    group: Mapping[str, Any],
    *,
    reference_entry: str,
    vertical_scale: str,
    plot_native_scale: bool,
) -> Optional[Dict[str, np.ndarray]]:
    if not group.get("plot_molecular", False):
        return None

    if plot_native_scale:
        entry = group.get("entries", {}).get(reference_entry, {})
        profile = entry.get("molecular", {}).get("profile")
        vertical = entry.get("metadata", {}).get(vertical_scale)
        if not isinstance(profile, xr.DataArray) or not isinstance(vertical, xr.DataArray):
            return None
        return {
            "x": _to_plot_units(vertical, vertical_scale),
            "y": np.asarray(profile.values, dtype=float),
            "label": "reference molecular",
        }

    molecular = group.get("molecular_profiles")
    vertical = group.get("vertical_grid")
    if not isinstance(molecular, xr.DataArray) or not isinstance(vertical, xr.DataArray):
        return None
    entry_values = [str(v) for v in molecular["entry"].values]
    if reference_entry not in entry_values:
        return None

    return {
        "x": _to_plot_units(vertical, vertical_scale),
        "y": np.asarray(molecular.sel(entry=reference_entry).values, dtype=float),
        "label": "reference molecular",
    }


def _smooth_profiles(
    X: Mapping[str, np.ndarray],
    Y: Mapping[str, np.ndarray],
    YE: Mapping[str, Optional[np.ndarray]],
    group: Mapping[str, Any],
    *,
    apply_smoothing: bool,
) -> Tuple[Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    Y_out: Dict[str, np.ndarray] = {}
    YE_out: Dict[str, Optional[np.ndarray]] = {}

    settings = {
        "smooth": group.get("smooth", False),
        "smoothing_range": group.get("smoothing_range", []),
        "smoothing_window": group.get("smoothing_window"),
    }

    for dataset_id, y in Y.items():
        if apply_smoothing and settings["smooth"] and settings["smoothing_window"]:
            y_sm, y_std = smoothing(
                args=settings,
                x_vals=np.asarray(X[dataset_id], dtype=float),
                y_vals=np.asarray(y, dtype=float),
                err_type="std",
            )
            Y_out[dataset_id] = np.asarray(y_sm, dtype=float)
            YE_out[dataset_id] = np.asarray(y_std, dtype=float)
        else:
            Y_out[dataset_id] = np.asarray(y, dtype=float).copy()
            err = YE.get(dataset_id)
            YE_out[dataset_id] = None if err is None else np.asarray(err, dtype=float).copy()

    return Y_out, YE_out


def _profiles_for_axis_limits(
    X: Mapping[str, np.ndarray],
    Y: Mapping[str, np.ndarray],
    stage_errors: Mapping[str, Optional[np.ndarray]],
    group: Mapping[str, Any],
    *,
    plotted_data_are_binned: bool,
) -> Tuple[Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    """Return signal/error pairs used only for automatic axis limits.

    The signal used for limit estimation is always the plotted processed signal:
    conservatively binned values are used directly, while non-binned values may
    be smoothed (including a temporary smoothing pass when plotting smoothing is
    disabled).  SNR, however, is *always* evaluated with the uncertainty loaded
    from the stage and propagated through intercomparison preprocessing / vertical
    harmonization.  The local standard deviation calculated by plotting smoothing
    is deliberately not used as an SNR uncertainty.
    """
    stored_errors = {
        key: None if stage_errors.get(key) is None else np.asarray(stage_errors[key], dtype=float)
        for key in Y
    }

    if plotted_data_are_binned or group.get("smooth", False):
        return (
            {key: np.asarray(value, dtype=float) for key, value in Y.items()},
            stored_errors,
        )

    window = group.get("smoothing_window")
    if not window:
        return (
            {key: np.asarray(value, dtype=float) for key, value in Y.items()},
            stored_errors,
        )

    settings = {
        "smooth": True,
        "smoothing_range": group.get("smoothing_range", []),
        "smoothing_window": window,
    }

    out_y: Dict[str, np.ndarray] = {}
    for dataset_id, y in Y.items():
        y_sm, _ = smoothing(
            args=settings,
            x_vals=np.asarray(X[dataset_id], dtype=float),
            y_vals=np.asarray(y, dtype=float),
            err_type="std",
        )
        out_y[dataset_id] = np.asarray(y_sm, dtype=float)

    return out_y, stored_errors


def _smooth_molecular(
    molecular: Optional[Dict[str, np.ndarray]],
    group: Mapping[str, Any],
    *,
    apply_smoothing: bool,
):
    if molecular is None:
        return None
    if (
        not apply_smoothing
        or not group.get("smooth", False)
        or not group.get("smoothing_window")
    ):
        return molecular

    settings = {
        "smooth": True,
        "smoothing_range": group.get("smoothing_range", []),
        "smoothing_window": group.get("smoothing_window"),
    }
    y_sm, _ = smoothing(
        args=settings,
        x_vals=np.asarray(molecular["x"], dtype=float),
        y_vals=np.asarray(molecular["y"], dtype=float),
        err_type="sem",
    )
    out = dict(molecular)
    out["y"] = np.asarray(y_sm, dtype=float)
    return out



def _slice_profile_dicts(
    X: Mapping[str, np.ndarray],
    Y: Mapping[str, np.ndarray],
    YE: Mapping[str, Optional[np.ndarray]],
    x_lims,
):
    """Slice plotting dictionaries to one horizontal interval."""
    X_out: Dict[str, np.ndarray] = {}
    Y_out: Dict[str, np.ndarray] = {}
    YE_out: Dict[str, Optional[np.ndarray]] = {}

    low, high = float(x_lims[0]), float(x_lims[1])
    for entry_id, y in Y.items():
        x = np.asarray(X[entry_id], dtype=float)
        y = np.asarray(y, dtype=float)
        mask = np.isfinite(x) & (x >= low) & (x <= high)
        X_out[entry_id] = x[mask]
        Y_out[entry_id] = y[mask]

        err = YE.get(entry_id)
        if err is None:
            YE_out[entry_id] = None
        else:
            err = np.asarray(err, dtype=float)
            YE_out[entry_id] = err[mask]

    return X_out, Y_out, YE_out


def _slice_molecular(molecular, x_lims):
    if molecular is None:
        return None
    x = np.asarray(molecular["x"], dtype=float)
    y = np.asarray(molecular["y"], dtype=float)
    low, high = float(x_lims[0]), float(x_lims[1])
    mask = np.isfinite(x) & (x >= low) & (x <= high)
    out = dict(molecular)
    out["x"] = x[mask]
    out["y"] = y[mask]
    return out

def _finite_concat(values):
    pieces = []
    for value in values:
        arr = np.asarray(value, dtype=float).ravel()
        arr = arr[np.isfinite(arr)]
        if arr.size:
            pieces.append(arr)
    if not pieces:
        return np.asarray([], dtype=float)
    return np.concatenate(pieces)


def _auto_x_lims(X: Mapping[str, np.ndarray]) -> list:
    vals = _finite_concat(X.values())
    if vals.size == 0:
        return [0.0, 1.0]
    low = float(np.nanmin(vals))
    high = float(np.nanmax(vals))
    if high <= low:
        return [low, low + 1.0]
    return [low, high]


def _auto_y_lims(
    Y,
    YE,
    molecular,
    *,
    use_log_y_scale: bool,
    dynamic_range_floor: float = 1.0e-5,
    snr_threshold: float = 1.0,
) -> list:
    """Determine automatic intercomparison y-axis limits.

    The logic intentionally follows the Rayleigh-fit convention when a
    molecular profile is plotted: the upper limit is controlled by the
    measured-profile peak, while the lower limit is controlled by the
    molecular profile.  Without molecular data, logarithmic plots are
    protected against isolated low/noisy values by a fixed dynamic-range
    floor relative to the measured peak.
    """

    signal_vals = _finite_concat(Y.values())
    if signal_vals.size == 0:
        return [1.0e-6, 1.0] if use_log_y_scale else [0.0, 1.0]

    # Use only statistically meaningful signal points to determine the upper
    # automatic limit. This prevents isolated noisy spikes from controlling
    # the full plot range. SNR is defined as abs(signal) / uncertainty.
    snr_signal_pieces = []
    for dataset_id, y in Y.items():
        y = np.asarray(y, dtype=float)
        err = YE.get(dataset_id)
        if err is None:
            continue
        err = np.asarray(err, dtype=float)
        if err.shape != y.shape:
            continue
        valid = (
            np.isfinite(y)
            & np.isfinite(err)
            & (err > 0.0)
            & (np.abs(y) / err > snr_threshold)
        )
        if np.any(valid):
            snr_signal_pieces.append(y[valid])

    if snr_signal_pieces:
        signal_for_max = np.concatenate(snr_signal_pieces)
        signal_max = float(np.nanmax(signal_for_max))
    else:
        print(
            f"      Warning: no finite signal points with SNR > {snr_threshold:g}; "
            "falling back to all finite signal values for automatic y limits"
        )
        signal_max = float(np.nanmax(signal_vals))

    if not np.isfinite(signal_max):
        return [1.0e-6, 1.0] if use_log_y_scale else [0.0, 1.0]

    # Rayleigh-fit-style upper margin.  This deliberately uses the measured
    # profiles, not the molecular profile, so a molecular tail cannot control
    # the upper plotting limit.
    upper_margin_factor = 2.5
    if signal_max > 0.0:
        high = upper_margin_factor * signal_max
    else:
        high = 1.0

    molecular_vals = np.asarray([], dtype=float)
    if molecular is not None:
        molecular_vals = _finite_concat([molecular.get("y", [])])

    if molecular_vals.size:
        # Match the Rayleigh-fit convention: the molecular profile determines
        # the lower limit, with a factor-of-two margin.
        if use_log_y_scale:
            molecular_positive = molecular_vals[molecular_vals > 0.0]
            if molecular_positive.size:
                low = float(np.nanmin(molecular_positive)) / 2.0
            else:
                low = max(high * dynamic_range_floor, 1.0e-12)
        else:
            low = float(np.nanmin(molecular_vals)) / 2.0

    elif use_log_y_scale:
        # No molecular reference is available.  Ignore non-positive samples
        # and prevent a few tiny positive noise excursions from creating an
        # excessive logarithmic dynamic range.
        positive = signal_vals[signal_vals > 0.0]
        if positive.size == 0:
            print(
                "      Warning: no positive values available for logarithmic "
                "y scale; using fallback limits"
            )
            return [1.0e-6, 1.0]

        low_from_data = float(np.nanmin(positive)) / 2.0
        low_from_dynamic_range = high * dynamic_range_floor
        low = max(low_from_data, low_from_dynamic_range)

    else:
        # Without molecular data, negative excursions are normally noise for
        # the intercomparison products, so keep the linear scale anchored at 0.
        low = 0.0

    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        if use_log_y_scale:
            high = max(high, 1.0)
            low = min(max(high * dynamic_range_floor, 1.0e-12), high / 10.0)
        else:
            low = 0.0
            high = max(signal_max * upper_margin_factor, 1.0)

    return [low, high]


def _nice_tick_spacing(span: float, target_intervals: int = 5) -> float:
    """Return a readable major-tick spacing for a displayed x interval.

    The value is chosen from 1, 2, 2.5, 5, or 10 times a power of ten,
    targeting roughly ``target_intervals`` major intervals.
    """
    span = float(span)
    if not np.isfinite(span) or span <= 0.0:
        return 1.0

    raw = span / max(int(target_intervals), 1)
    exponent = np.floor(np.log10(raw))
    scale = 10.0 ** exponent
    fraction = raw / scale

    for candidate in (1.0, 2.0, 2.5, 5.0, 10.0):
        if fraction <= candidate:
            return float(candidate * scale)

    return float(10.0 * scale)


def _near_range_y_lims(Y: Mapping[str, np.ndarray]) -> list:
    """Return linear y limits for the near-range signal panel.

    If every finite measured signal value in the displayed near-range interval
    is non-negative, anchor the axis at zero and add 10 percent headroom above
    the maximum.  Otherwise extend both extrema outward by 10 percent of their
    absolute magnitudes.
    """
    vals = _finite_concat(Y.values())
    if vals.size == 0:
        return [0.0, 1.0]

    low = float(np.nanmin(vals))
    high = float(np.nanmax(vals))

    if low >= 0.0:
        if high <= 0.0 or not np.isfinite(high):
            return [0.0, 1.0]
        return [0.0, 1.10 * high]

    low_lim = low - 0.10 * abs(low)
    high_lim = high + 0.10 * abs(high)

    if not np.isfinite(low_lim) or not np.isfinite(high_lim) or high_lim <= low_lim:
        extent = max(abs(low), abs(high), 1.0)
        return [-1.10 * extent, 1.10 * extent]

    return [low_lim, high_lim]


def _differences_to_reference(
    Y: Mapping[str, np.ndarray],
    YE: Mapping[str, Optional[np.ndarray]],
    *,
    reference_entry: str,
    difference_mode: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    """Calculate entry differences to the reference.

    Channel groups use relative differences, while pair groups use absolute
    differences. Uncertainties from the entry and reference are treated as
    independent and propagated in quadrature.
    """
    if reference_entry not in Y:
        raise ValueError(
            f"Reference entry {reference_entry!r} is missing from the plotted group"
        )

    y_ref = np.asarray(Y[reference_entry], dtype=float)
    e_ref = YE.get(reference_entry)
    e_ref_arr = None if e_ref is None else np.asarray(e_ref, dtype=float)

    differences: Dict[str, np.ndarray] = {}
    difference_error: Dict[str, Optional[np.ndarray]] = {}

    for entry_id, y in Y.items():
        if entry_id == reference_entry:
            continue

        y = np.asarray(y, dtype=float)
        if y.shape != y_ref.shape:
            raise ValueError(
                f"Cannot calculate differences for {entry_id!r}: its shape {y.shape} "
                f"does not match the reference shape {y_ref.shape}."
            )

        diff = np.full_like(y, np.nan, dtype=float)

        if difference_mode == "absolute":
            valid = np.isfinite(y) & np.isfinite(y_ref)
            diff[valid] = y[valid] - y_ref[valid]
        elif difference_mode == "relative":
            valid = np.isfinite(y) & np.isfinite(y_ref) & (y_ref != 0)
            diff[valid] = (y[valid] - y_ref[valid]) / y_ref[valid]
        else:
            raise ValueError(
                f"Unsupported difference_mode {difference_mode!r}; "
                "expected 'relative' or 'absolute'"
            )

        differences[entry_id] = diff

        err = YE.get(entry_id)
        if err is None and e_ref_arr is None:
            difference_error[entry_id] = None
            continue

        err_arr = (
            np.zeros_like(y, dtype=float)
            if err is None
            else np.asarray(err, dtype=float)
        )
        ref_err_arr = (
            np.zeros_like(y_ref, dtype=float)
            if e_ref_arr is None
            else e_ref_arr
        )

        sigma = np.full_like(y, np.nan, dtype=float)
        valid_err = valid & np.isfinite(err_arr) & np.isfinite(ref_err_arr)

        if difference_mode == "absolute":
            sigma[valid_err] = np.sqrt(
                err_arr[valid_err] ** 2 + ref_err_arr[valid_err] ** 2
            )
        else:
            sigma[valid_err] = np.sqrt(
                (err_arr[valid_err] / y_ref[valid_err]) ** 2
                + (
                    y[valid_err] * ref_err_arr[valid_err]
                    / (y_ref[valid_err] ** 2)
                ) ** 2
            )

        difference_error[entry_id] = sigma

    return differences, difference_error


def _auto_absolute_difference_lims(
    differences: Mapping[str, np.ndarray],
    Y: Mapping[str, np.ndarray],
    YE: Mapping[str, Optional[np.ndarray]],
    *,
    reference_entry: str,
    snr_threshold: float = 1.0,
) -> list:
    """Return robust symmetric limits for pair absolute differences.

    Only locations where both the compared entry and reference have SNR above
    the threshold are allowed to control the automatic limit. The complete
    difference curve is still plotted; this filter is only for axis scaling.
    """
    y_ref = np.asarray(Y[reference_entry], dtype=float)
    e_ref = YE.get(reference_entry)
    e_ref = None if e_ref is None else np.asarray(e_ref, dtype=float)

    pieces = []
    for entry_id, diff in differences.items():
        y = np.asarray(Y[entry_id], dtype=float)
        e = YE.get(entry_id)
        e = None if e is None else np.asarray(e, dtype=float)
        diff = np.asarray(diff, dtype=float)

        if e is None or e_ref is None or e.shape != y.shape or e_ref.shape != y_ref.shape:
            continue

        valid = (
            np.isfinite(diff)
            & np.isfinite(y)
            & np.isfinite(y_ref)
            & np.isfinite(e)
            & np.isfinite(e_ref)
            & (e > 0.0)
            & (e_ref > 0.0)
            & (np.abs(y) / e > snr_threshold)
            & (np.abs(y_ref) / e_ref > snr_threshold)
        )
        if np.any(valid):
            pieces.append(diff[valid])

    if pieces:
        vals = np.concatenate(pieces)
    else:
        print(
            f"      Warning: no pair-difference points with both SNR > {snr_threshold:g}; "
            "falling back to all finite differences for automatic limits"
        )
        vals = _finite_concat(differences.values())

    if vals.size == 0:
        return [-1.0, 1.0]

    extent = float(np.nanmax(np.abs(vals)))
    if not np.isfinite(extent) or extent == 0.0:
        extent = 1.0
    extent *= 1.10
    return [-extent, extent]


def _plot_one_group(
    intercomparison_info: Mapping[str, Any],
    group_id: str,
    group: Mapping[str, Any],
    *,
    group_kind: str,
) -> str:
    general = intercomparison_info["general"]
    vertical_scale = general["vertical_scale"]
    plot_native_scale = bool(general.get("plot_native_scale", False))
    reference_entry = group.get("reference_entry")

    # Plotting/processing configuration belongs to intercomparison_info, while
    # the bundle group contains the loaded/processed data.  Merge the two for
    # this plotting call so newly added plotting settings do not have to be
    # duplicated inside prepare_intercomparison_bundles().
    info_group_key = "channel_groups" if group_kind == "channel_group" else "pair_groups"
    configured_group = intercomparison_info.get(info_group_key, {}).get(group_id, {})
    group_settings = dict(configured_group)
    group_settings.update(group)

    if plot_native_scale:
        X, Y, YE = _native_group_arrays(group, vertical_scale=vertical_scale)
    else:
        X, Y, YE = _harmonized_group_arrays(group, vertical_scale=vertical_scale)

    if not Y:
        raise ValueError(f"[{group_kind}:{group_id}] contains no plottable signals")

    # Keep the stage-derived / preprocessing-propagated uncertainty separate
    # from any local STD later calculated only for plotting smoothing.  Automatic
    # SNR-based axis limits must use this stored uncertainty.
    YE_stage = {
        key: None if YE.get(key) is None else np.asarray(YE[key], dtype=float).copy()
        for key in Y
    }

    # Preserve the processed, unsmoothed plotting arrays.  The near-range row
    # starts from these arrays and applies its own independent smoothing.
    X_base = {key: np.asarray(value, dtype=float).copy() for key, value in X.items()}
    Y_base = {key: np.asarray(value, dtype=float).copy() for key, value in Y.items()}
    YE_base = {
        key: None if YE.get(key) is None else np.asarray(YE[key], dtype=float).copy()
        for key in Y
    }

    # Smoothing/local-STD estimation is only useful when the plotted profiles
    # have not already been conservatively binned. Binning already provides a
    # weighted mean and propagated uncertainty. Native-scale plotting is also
    # considered unbinned, even if vertical_method=vertical_binning, because
    # it deliberately bypasses the harmonized binned arrays.
    vertical_method = general["vertical_method"]
    plotted_data_are_binned = (
        not plot_native_scale
        and vertical_method == "vertical_binning"
    )
    apply_smoothing = not plotted_data_are_binned

    Y, YE = _smooth_profiles(
        X, Y, YE, group_settings, apply_smoothing=apply_smoothing
    )
    molecular = _reference_molecular(
        group,
        reference_entry=reference_entry,
        vertical_scale=vertical_scale,
        plot_native_scale=plot_native_scale,
    )
    molecular = _smooth_molecular(
        molecular, group_settings, apply_smoothing=apply_smoothing
    )

    use_log_y_scale = bool(group_settings.get("use_log_y_scale", False))

    axis_limit_Y, axis_limit_YE = _profiles_for_axis_limits(
        X,
        Y,
        YE_stage,
        group_settings,
        plotted_data_are_binned=plotted_data_are_binned,
    )

    x_lims = list(group_settings.get("x_lims") or _auto_x_lims(X))
    auto_y_limits = not bool(group_settings.get("y_lims"))
    y_lims = list(
        group_settings.get("y_lims")
        or _auto_y_lims(
            axis_limit_Y,
            axis_limit_YE,
            molecular,
            use_log_y_scale=use_log_y_scale,
        )
    )

    # Pair ratios are linear by default and their physically meaningful lower
    # plotting bound is zero. Keep explicit user y_lims untouched.
    if group_kind == "pair_group" and auto_y_limits and not use_log_y_scale:
        y_lims[0] = 0.0

    difference_mode = (
        "relative" if group_kind == "channel_group" else "absolute"
    )

    if plot_native_scale:
        differences = {}
        difference_error = {}
    else:
        differences, difference_error = _differences_to_reference(
            Y,
            YE,
            reference_entry=reference_entry,
            difference_mode=difference_mode,
        )

    lib = IntercomparisonLibraries(
        intercomparison_info=dict(intercomparison_info),
        group_id=group_id,
        group_kind=group_kind,
        group=dict(group_settings),
    )
    text_generator = GenerateIntercomparisonText(lib)

    entry_labels = {
        entry_id: _entry_label(intercomparison_info, entry_id, entry)
        for entry_id, entry in group.get("entries", {}).items()
    }
    entry_color_indices = {
        entry_id: entry.get("color_index")
        for entry_id, entry in group.get("entries", {}).items()
        if entry.get("color_index") is not None
    }

    normalisation_region = None
    if vertical_scale in PHYSICAL_VERTICAL_SCALES and group_settings.get("normalisation"):
        normalisation_region = group_settings.get("normalisation_region")

    configured_difference_lims = group_settings.get("difference_y_lims")

    if configured_difference_lims:
        difference_lims = list(configured_difference_lims)
    elif difference_mode == "absolute":
        difference_lims = _auto_absolute_difference_lims(
            differences,
            Y,
            YE_stage,
            reference_entry=reference_entry,
        )
    else:
        # Channel relative-difference defaults are normally resolved by the
        # parser, but retain an automatic fallback for robustness.
        difference_lims = [-0.4, 0.4]

    # ------------------------------------------------------------------
    # Second row: near-range comparison with independent plotting controls.
    # The near row mirrors the full-range controls (x/y limits, tick spacing,
    # smoothing, and difference limits) while remaining linear by design.
    # ------------------------------------------------------------------
    near_x_lims = list(group_settings.get("near_x_lims") or [0.0, 2.0])
    if float(near_x_lims[1]) <= float(near_x_lims[0]):
        raise ValueError(
            f"[{group_kind}:{group_id}] near_x_lims must be increasing, got {near_x_lims}"
        )

    near_settings = dict(group_settings)
    near_smooth = group_settings.get("near_smooth")
    if near_smooth is None:
        # Normally resolved by config.py, but keep this fallback explicit.
        near_smooth = False
    near_settings["smooth"] = bool(near_smooth)
    near_settings["smoothing_range"] = list(
        group_settings.get("near_smoothing_range") or near_x_lims
    )
    near_settings["smoothing_window"] = float(
        group_settings.get("near_smoothing_window", 0.1)
    )

    # Keep near-range smoothing completely independent from the full-range row.
    # In particular, pair_near_smooth/channel_near_smooth=False must bypass the
    # smoothing routine entirely, rather than relying on the helper to inspect
    # the flag internally.
    if near_settings["smooth"]:
        # Smooth before clipping to the displayed near-range interval. If the
        # smoothing implementation cannot evaluate the full window at an edge,
        # keep the original finite sample there rather than shortening the line.
        near_Y_full, near_YE_full = _smooth_profiles(
            X_base, Y_base, YE_base, near_settings, apply_smoothing=True
        )
        for key in near_Y_full:
            smoothed = np.asarray(near_Y_full[key], dtype=float)
            original = np.asarray(Y_base[key], dtype=float)
            missing = ~np.isfinite(smoothed) & np.isfinite(original)
            if np.any(missing):
                smoothed = smoothed.copy()
                smoothed[missing] = original[missing]
                near_Y_full[key] = smoothed

                smoothed_err = near_YE_full.get(key)
                original_err = YE_base.get(key)
                if smoothed_err is not None and original_err is not None:
                    smoothed_err = np.asarray(smoothed_err, dtype=float).copy()
                    original_err = np.asarray(original_err, dtype=float)
                    fill_err = (
                        missing
                        & ~np.isfinite(smoothed_err)
                        & np.isfinite(original_err)
                    )
                    smoothed_err[fill_err] = original_err[fill_err]
                    near_YE_full[key] = smoothed_err
    else:
        # No smoothing means exactly the processed unsmoothed arrays and the
        # stage/preprocessing-propagated errors are used in the near-range row.
        near_Y_full = {
            key: np.asarray(value, dtype=float).copy()
            for key, value in Y_base.items()
        }
        near_YE_full = {
            key: (
                None
                if YE_base.get(key) is None
                else np.asarray(YE_base[key], dtype=float).copy()
            )
            for key in Y_base
        }

    near_X, near_Y, near_YE = _slice_profile_dicts(
        X_base, near_Y_full, near_YE_full, near_x_lims
    )

    near_molecular_full = _reference_molecular(
        group,
        reference_entry=reference_entry,
        vertical_scale=vertical_scale,
        plot_native_scale=plot_native_scale,
    )
    if near_settings["smooth"]:
        near_molecular_full = _smooth_molecular(
            near_molecular_full, near_settings, apply_smoothing=True
        )
    near_molecular = _slice_molecular(near_molecular_full, near_x_lims)

    # The near-range signal panel stays linear, independently of the full-range
    # channel-group scale.  Empty near_y_lims retains the automatic near-range
    # rule; explicit limits override it.
    configured_near_y_lims = group_settings.get("near_y_lims")
    near_y_lims = (
        list(configured_near_y_lims)
        if configured_near_y_lims
        else _near_range_y_lims(near_Y)
    )

    if plot_native_scale:
        near_differences = {}
        near_difference_error = {}
    else:
        near_differences, near_difference_error = _differences_to_reference(
            near_Y,
            near_YE,
            reference_entry=reference_entry,
            difference_mode=difference_mode,
        )

    configured_near_difference_lims = group_settings.get("near_difference_y_lims")
    if configured_near_difference_lims:
        near_difference_lims = list(configured_near_difference_lims)
    elif difference_mode == "absolute":
        near_stage_errors = {}
        for key in near_Y:
            if YE_stage.get(key) is None:
                near_stage_errors[key] = None
                continue
            x_base = np.asarray(X_base[key], dtype=float)
            mask = (
                np.isfinite(x_base)
                & (x_base >= near_x_lims[0])
                & (x_base <= near_x_lims[1])
            )
            near_stage_errors[key] = np.asarray(YE_stage[key], dtype=float)[mask]
        near_difference_lims = _auto_absolute_difference_lims(
            near_differences,
            near_Y,
            near_stage_errors,
            reference_entry=reference_entry,
        )
    else:
        near_difference_lims = list(difference_lims)

    args = {
        "title": _figure_title(
            intercomparison_info, group_id, group_kind, group_settings
        ),
        "filename": text_generator.make_filename(),
        "plot_folder": str(Path(general["output_folder"]) / "plots"),
        "dpi": general["dpi"],
        "vertical_scale": vertical_scale,
        "plot_native_scale": plot_native_scale,
        "x_lims": x_lims,
        "x_tick": float(group_settings["x_tick"]),
        "y_lims": y_lims,
        "difference_lims": difference_lims,
        "difference_mode": difference_mode,
        "use_log_y_scale": use_log_y_scale,
        "normalisation_region": normalisation_region,
        "entry_labels": entry_labels,
        "entry_color_indices": entry_color_indices,
        "left_y_label": str(
            group_settings.get("quantity_label")
            or ("Signal" if group_kind == "channel_group" else "VLDR")
        ),
    }

    near_x_tick = float(
        group_settings.get("near_x_tick")
        or _nice_tick_spacing(near_x_lims[1] - near_x_lims[0])
    )

    near_args = dict(args)
    near_args.update({
        "x_lims": near_x_lims,
        "x_tick": near_x_tick,
        "y_lims": near_y_lims,
        "difference_lims": near_difference_lims,
        "use_log_y_scale": False,
    })

    print(f"    - {group_id}")
    if plotted_data_are_binned:
        print(
            "        smoothing: skipped (conservative binning already applied; "
            "using binned values and propagated errors)"
        )
    else:
        print(f"        smoothing: {'on' if group_settings.get('smooth') else 'off'}")
        if group_settings.get("smooth"):
            print(
                f"        smoothing range: {group_settings.get('smoothing_range')}; "
                f"window: {group_settings.get('smoothing_window')}"
            )
    print(f"        x limits: {args['x_lims']}")
    print(f"        y limits: {args['y_lims']}")
    print(
        f"        near range: {near_x_lims}; x tick: {near_x_tick:g}; "
        f"smoothing: {'on' if near_settings['smooth'] else 'off'}; "
        f"range: {near_settings['smoothing_range']}; "
        f"window: {near_settings['smoothing_window']}; y scale: linear; "
        f"y limits: {near_y_lims}"
    )
    if not group_settings.get("y_lims"):
        print(
            "        axis-limit signal basis: "
            + ("binned values" if plotted_data_are_binned else "smoothed values")
            + " with SNR > 1 using stage-derived processed errors"
        )
    if group_kind == "channel_group":
        print(f"        y scale: {'log' if use_log_y_scale else 'linear'}")
    else:
        print("        y scale: linear")
    if plot_native_scale:
        print("        right panel: empty (native scales are not aligned)")
    else:
        print(
            f"        right panel: {difference_mode} differences to "
            f"{reference_entry!r}"
        )

    return plot_intercomparison.generate_plot(
        X=X,
        Y=Y,
        YE=YE,
        differences=differences,
        difference_error=difference_error,
        molecular=molecular,
        args=args,
        near_X=near_X,
        near_Y=near_Y,
        near_YE=near_YE,
        near_differences=near_differences,
        near_difference_error=near_difference_error,
        near_molecular=near_molecular,
        near_args=near_args,
    )


def generate_intercomparison(
    intercomparison_info: Mapping[str, Any],
    intercomparison_bundles: Mapping[str, Any],
) -> Dict[str, Dict[str, str]]:
    """Generate one two-panel figure for every active channel/pair group."""
    output: Dict[str, Dict[str, str]] = {
        "channel_groups": {},
        "pair_groups": {},
    }

    print()
    print("-----------------------------------------------")
    print("Intercomparison plots")
    print("-----------------------------------------------")

    for group_id, group in intercomparison_bundles.get("channel_groups", {}).items():
        output["channel_groups"][group_id] = _plot_one_group(
            intercomparison_info,
            group_id,
            group,
            group_kind="channel_group",
        )

    for group_id, group in intercomparison_bundles.get("pair_groups", {}).items():
        output["pair_groups"][group_id] = _plot_one_group(
            intercomparison_info,
            group_id,
            group,
            group_kind="pair_group",
        )

    print("Intercomparison plotting complete.")
    return output
