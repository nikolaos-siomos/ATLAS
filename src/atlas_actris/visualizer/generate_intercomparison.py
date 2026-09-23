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

    args = {
        "title": text_generator.make_title(),
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
        "left_y_label": (
            "Signal"
            if group_kind == "channel_group"
            else "Pair ratio"
        ),
    }

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
