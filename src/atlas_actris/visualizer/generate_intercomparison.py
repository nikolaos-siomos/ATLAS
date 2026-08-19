#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare and generate ATLAS intercomparison group plots.

This module plays the role that ``qa_test_rayleigh_fit.py`` plays for the
Rayleigh-fit plot: it extracts one comparison group, applies optional plotting
smoothing / local-STD estimation, determines plot limits, prepares relative
comparisons to the reference dataset, creates text, and calls
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


def _dataset_label(intercomparison_info: Mapping[str, Any], dataset_id: str) -> str:
    dataset = intercomparison_info["datasets"].get(dataset_id, {})
    return dataset.get("dataset_label") or dataset.get("system_label") or dataset_id


def _native_group_arrays(
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    X: Dict[str, np.ndarray] = {}
    Y: Dict[str, np.ndarray] = {}
    YE: Dict[str, Optional[np.ndarray]] = {}

    for dataset_id, dataset in group.get("datasets", {}).items():
        signal = dataset.get("signal")
        vertical = dataset.get("metadata", {}).get(vertical_scale)
        if not isinstance(signal, xr.DataArray) or not isinstance(vertical, xr.DataArray):
            continue

        X[dataset_id] = _to_plot_units(vertical, vertical_scale)
        Y[dataset_id] = np.asarray(signal.values, dtype=float)
        YE[dataset_id] = _array_values(dataset.get("error"))

    return X, Y, YE


def _harmonized_group_arrays(
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    signals = group.get("signals")
    if not isinstance(signals, xr.DataArray) or "dataset" not in signals.dims:
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

    dataset_ids = [str(v) for v in signals["dataset"].values]
    for dataset_id in dataset_ids:
        X[dataset_id] = x_common.copy()
        Y[dataset_id] = np.asarray(signals.sel(dataset=dataset_id).values, dtype=float)
        if isinstance(errors, xr.DataArray) and dataset_id in [str(v) for v in errors["dataset"].values]:
            YE[dataset_id] = np.asarray(errors.sel(dataset=dataset_id).values, dtype=float)
        else:
            YE[dataset_id] = None

    return X, Y, YE


def _reference_molecular(
    group: Mapping[str, Any],
    *,
    reference_dataset: str,
    vertical_scale: str,
    plot_native_scale: bool,
) -> Optional[Dict[str, np.ndarray]]:
    if not group.get("plot_molecular", False):
        return None

    if plot_native_scale:
        dataset = group.get("datasets", {}).get(reference_dataset, {})
        profile = dataset.get("molecular", {}).get("profile")
        vertical = dataset.get("metadata", {}).get(vertical_scale)
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
    dataset_values = [str(v) for v in molecular["dataset"].values]
    if reference_dataset not in dataset_values:
        return None

    return {
        "x": _to_plot_units(vertical, vertical_scale),
        "y": np.asarray(molecular.sel(dataset=reference_dataset).values, dtype=float),
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


def _auto_y_lims(Y, molecular, *, use_log_y_scale: bool) -> list:
    values = list(Y.values())
    if molecular is not None:
        values.append(molecular["y"])
    vals = _finite_concat(values)
    if vals.size == 0:
        return [1.0e-6, 1.0] if use_log_y_scale else [0.0, 1.0]

    if use_log_y_scale:
        positive = vals[vals > 0]
        if positive.size == 0:
            print("      Warning: no positive values available for logarithmic y scale; using fallback limits")
            return [1.0e-6, 1.0]
        low = float(np.nanmin(positive)) / 2.0
        high = float(np.nanmax(positive)) * 2.0
        if high <= low:
            high = low * 10.0
        return [low, high]

    low = float(np.nanmin(vals))
    high = float(np.nanmax(vals))
    if high <= low:
        spread = max(abs(low) * 0.1, 1.0)
        return [low - spread, high + spread]
    pad = 0.05 * (high - low)
    return [low - pad, high + pad]


def _relative_to_reference(
    Y: Mapping[str, np.ndarray],
    YE: Mapping[str, Optional[np.ndarray]],
    *,
    reference_dataset: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, Optional[np.ndarray]]]:
    if reference_dataset not in Y:
        raise ValueError(f"Reference dataset {reference_dataset!r} is missing from the plotted group")

    y_ref = np.asarray(Y[reference_dataset], dtype=float)
    e_ref = YE.get(reference_dataset)
    e_ref_arr = None if e_ref is None else np.asarray(e_ref, dtype=float)

    relative: Dict[str, np.ndarray] = {}
    relative_error: Dict[str, Optional[np.ndarray]] = {}

    for dataset_id, y in Y.items():
        if dataset_id == reference_dataset:
            continue
        y = np.asarray(y, dtype=float)
        if y.shape != y_ref.shape:
            raise ValueError(
                f"Cannot calculate differences for {dataset_id!r}: its shape {y.shape} "
                f"does not match the reference shape {y_ref.shape}."
            )

        rel = np.full_like(y, np.nan, dtype=float)
        valid = np.isfinite(y) & np.isfinite(y_ref) & (y_ref != 0)
        rel[valid] = (y[valid] - y_ref[valid]) / y_ref[valid]
        relative[dataset_id] = rel

        err = YE.get(dataset_id)
        if err is None and e_ref_arr is None:
            relative_error[dataset_id] = None
            continue

        err_arr = np.zeros_like(y, dtype=float) if err is None else np.asarray(err, dtype=float)
        ref_err_arr = np.zeros_like(y_ref, dtype=float) if e_ref_arr is None else e_ref_arr
        sigma = np.full_like(y, np.nan, dtype=float)
        valid_err = valid & np.isfinite(err_arr) & np.isfinite(ref_err_arr)
        sigma[valid_err] = np.sqrt(
            (err_arr[valid_err] / y_ref[valid_err]) ** 2
            + (
                y[valid_err] * ref_err_arr[valid_err]
                / (y_ref[valid_err] ** 2)
            ) ** 2
        )
        relative_error[dataset_id] = sigma

    return relative, relative_error


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
    reference_dataset = intercomparison_info["reference_dataset"]

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
        reference_dataset=reference_dataset,
        vertical_scale=vertical_scale,
        plot_native_scale=plot_native_scale,
    )
    molecular = _smooth_molecular(
        molecular, group_settings, apply_smoothing=apply_smoothing
    )

    use_log_y_scale = bool(group_settings.get("use_log_y_scale", False))

    x_lims = list(group_settings.get("x_lims") or _auto_x_lims(X))
    y_lims = list(
        group_settings.get("y_lims")
        or _auto_y_lims(Y, molecular, use_log_y_scale=use_log_y_scale)
    )

    if plot_native_scale:
        relative = {}
        relative_error = {}
    else:
        relative, relative_error = _relative_to_reference(
            Y,
            YE,
            reference_dataset=reference_dataset,
        )

    lib = IntercomparisonLibraries(
        intercomparison_info=dict(intercomparison_info),
        group_id=group_id,
        group_kind=group_kind,
        group=dict(group_settings),
    )
    text_generator = GenerateIntercomparisonText(lib)

    dataset_labels = {
        dataset_id: _dataset_label(intercomparison_info, dataset_id)
        for dataset_id in group.get("datasets", {})
    }

    normalisation_region = None
    if vertical_scale in PHYSICAL_VERTICAL_SCALES and group_settings.get("normalisation"):
        normalisation_region = group_settings.get("normalisation_region")

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
        "relative_difference_lims": list(group_settings["relative_difference_lims"]),
        "use_log_y_scale": use_log_y_scale,
        "normalisation_region": normalisation_region,
        "dataset_labels": dataset_labels,
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
    if group_kind == "channel_group":
        print(f"        y scale: {'log' if use_log_y_scale else 'linear'}")
    else:
        print("        y scale: linear")
    if plot_native_scale:
        print("        right panel: empty (native scales are not aligned)")
    else:
        print(f"        right panel: relative differences to {reference_dataset!r}")

    return plot_intercomparison.generate_plot(
        X=X,
        Y=Y,
        YE=YE,
        relative=relative,
        relative_error=relative_error,
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
