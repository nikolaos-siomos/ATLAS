#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Vertical harmonization for ATLAS intercomparison bundles.

``general.vertical_method`` selects either interpolation onto the coarsest
native dataset grid within each group or conservative overlap-weighted vertical binning.
``general.vertical_scale`` selects the coordinate used by either method.

For vertical binning, bin width is resolved group-wise:
1. group ``vertical_bin_width`` override;
2. ``[general] vertical_bin_width``;
3. if the requested width is smaller than the coarsest nominal native step,
   warn and use the coarsest native step instead.

For physical scales (range, height_agl, height_asl), first_bin_left_edge,
vertical_min, vertical_max and vertical_bin_width are specified in kilometres
and converted to metres internally. For ``bins`` they are interpreted directly
in bin units.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Tuple

import numpy as np
import xarray as xr


VERTICAL_SCALES = ("bins", "range", "height_agl", "height_asl")
VERTICAL_METHODS = ("interpolation", "vertical_binning")
PHYSICAL_VERTICAL_SCALES = ("range", "height_agl", "height_asl")
VERTICAL_METADATA_KEYS = ("bins", "range", "height_agl", "height_asl")


def _copy_nested(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _copy_nested(item) for key, item in value.items()}
    return value


def _to_internal_units(value: Optional[float], vertical_scale: str) -> Optional[float]:
    if value is None:
        return None
    value = float(value)
    if vertical_scale in PHYSICAL_VERTICAL_SCALES:
        return 1000.0 * value
    return value


def _from_internal_units(value: float, vertical_scale: str) -> float:
    if vertical_scale in PHYSICAL_VERTICAL_SCALES:
        return float(value) / 1000.0
    return float(value)


def _units_label(vertical_scale: str) -> str:
    return "km" if vertical_scale in PHYSICAL_VERTICAL_SCALES else "bins"


def _as_1d_vertical(value: Any, *, location: str) -> xr.DataArray:
    if not isinstance(value, xr.DataArray):
        raise ValueError(f"{location}: vertical coordinate is missing")
    if "bins" not in value.dims:
        raise ValueError(f"{location}: vertical coordinate has no 'bins' dimension")

    extra_dims = [dim for dim in value.dims if dim != "bins"]
    for dim in extra_dims:
        if value.sizes[dim] != 1:
            raise ValueError(
                f"{location}: vertical coordinate must be one-dimensional after "
                f"product selection; dimension {dim!r} has size {value.sizes[dim]}"
            )
        value = value.squeeze(dim, drop=True)
    return value


def _valid_vertical_values(vertical: xr.DataArray) -> np.ndarray:
    values = np.asarray(vertical.values, dtype=float)
    return values[np.isfinite(values)]


def _dataset_bounds(
    dataset: Mapping[str, Any],
    *,
    vertical_scale: str,
    location: str,
) -> Tuple[float, float]:
    vertical = _as_1d_vertical(
        dataset.get("metadata", {}).get(vertical_scale),
        location=location,
    )
    values = _valid_vertical_values(vertical)
    if values.size == 0:
        raise ValueError(f"{location}: vertical coordinate contains no finite values")
    return float(values.min()), float(values.max())


def _common_bounds(
    group_kind: str,
    group_id: str,
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
    vertical_min: Optional[float],
    vertical_max: Optional[float],
) -> Tuple[float, float]:
    bounds = [
        _dataset_bounds(
            dataset,
            vertical_scale=vertical_scale,
            location=(
                f"[{group_kind}:{group_id}] entry {entry_id!r} "
                f"{vertical_scale!r}"
            ),
        )
        for entry_id, dataset in group.get("entries", {}).items()
    ]
    if not bounds:
        raise ValueError(f"[{group_kind}:{group_id}] contains no participating entrys")

    low = max(bound[0] for bound in bounds)
    high = min(bound[1] for bound in bounds)

    if vertical_min is not None:
        low = max(low, vertical_min)
    if vertical_max is not None:
        high = min(high, vertical_max)

    if not low < high:
        raise ValueError(
            f"[{group_kind}:{group_id}] entries have no common valid overlap "
            f"on vertical_scale={vertical_scale!r} after applying vertical limits"
        )
    return low, high


def _prepare_for_interp(
    value: xr.DataArray,
    vertical: xr.DataArray,
    *,
    location: str,
) -> xr.DataArray:
    if "bins" not in value.dims:
        raise ValueError(f"{location}: value has no 'bins' dimension")

    vertical = _as_1d_vertical(vertical, location=f"{location} vertical coordinate")
    if value.sizes["bins"] != vertical.sizes["bins"]:
        raise ValueError(
            f"{location}: value has {value.sizes['bins']} bins but its vertical "
            f"coordinate has {vertical.sizes['bins']}"
        )

    z = np.asarray(vertical.values, dtype=float)
    finite = np.isfinite(z)
    if finite.sum() < 2:
        raise ValueError(f"{location}: fewer than two finite vertical points are available")

    positions = np.flatnonzero(finite)
    z_finite = z[positions]
    order = np.argsort(z_finite)
    positions = positions[order]
    z_sorted = z_finite[order]

    # Remove duplicate vertical coordinates; xarray interpolation requires a
    # unique monotonic index. Keep the first occurrence of each value.
    _, unique_pos = np.unique(z_sorted, return_index=True)
    unique_pos = np.sort(unique_pos)
    positions = positions[unique_pos]
    z_sorted = z_sorted[unique_pos]

    out = value.isel(bins=positions)
    out = out.assign_coords(_vertical=("bins", z_sorted))
    return out.swap_dims({"bins": "_vertical"})


def _interp_value(
    value: Any,
    vertical: xr.DataArray,
    target_vertical: np.ndarray,
    *,
    location: str,
) -> Any:
    if not isinstance(value, xr.DataArray) or "bins" not in value.dims:
        return value

    prepared = _prepare_for_interp(value, vertical, location=location)
    target = xr.DataArray(target_vertical, dims=("_vertical",), coords={"_vertical": target_vertical})
    out = prepared.interp(_vertical=target)

    # ``swap_dims`` keeps the original ``bins`` coordinate as an auxiliary
    # coordinate. Drop it before renaming the temporary interpolation
    # dimension back to ``bins`` to avoid a name conflict in xarray.
    if "bins" in out.coords:
        out = out.drop_vars("bins")

    out = out.rename({"_vertical": "bins"})
    out = out.assign_coords(bins=np.arange(out.sizes["bins"], dtype=float))
    return out


def _interp_error(
    value: Any,
    vertical: xr.DataArray,
    target_vertical: np.ndarray,
    *,
    location: str,
) -> Any:
    """Linearly interpolate a value while propagating independent errors.

    If a target point lies between native points ``z0`` and ``z1``, the same
    linear interpolation weights used for the signal are applied to the two
    uncertainties in quadrature::

        sigma = sqrt((w0 * sigma0)**2 + (w1 * sigma1)**2)

    Thus an exact native point keeps its original uncertainty, while points
    between native bins have a reduced uncertainty. For approximately equal
    adjacent errors, a midpoint has ``sigma_out ~= sigma_native / sqrt(2)``.
    """
    if not isinstance(value, xr.DataArray) or "bins" not in value.dims:
        return value

    prepared = _prepare_for_interp(value, vertical, location=location)
    source_z = np.asarray(prepared["_vertical"].values, dtype=float)
    target_z = np.asarray(target_vertical, dtype=float)

    if source_z.size < 2:
        raise ValueError(f"{location}: fewer than two source points are available")

    # Locate the bracketing source points for every target coordinate.
    right = np.searchsorted(source_z, target_z, side="left")
    exact = (right < source_z.size) & np.isclose(
        source_z[np.clip(right, 0, source_z.size - 1)],
        target_z,
        rtol=0.0,
        atol=np.finfo(float).eps * 16.0,
    )

    outside = (target_z < source_z[0]) | (target_z > source_z[-1])
    right = np.clip(right, 1, source_z.size - 1)
    left = right - 1

    # Exact source coordinates should not be smoothed: use that native point
    # with unit weight and zero weight on the neighbour.
    exact_index = np.searchsorted(source_z, target_z, side="left")
    exact_index = np.clip(exact_index, 0, source_z.size - 1)
    left = np.where(exact, exact_index, left)
    right = np.where(exact, exact_index, right)

    z0 = source_z[left]
    z1 = source_z[right]
    span = z1 - z0
    w1 = np.zeros_like(target_z, dtype=float)
    nonzero = span != 0.0
    w1[nonzero] = (target_z[nonzero] - z0[nonzero]) / span[nonzero]
    w1 = np.where(exact, 0.0, w1)
    w0 = 1.0 - w1

    target_dim = xr.DataArray(
        np.arange(target_z.size, dtype=int),
        dims=("_target_vertical",),
    )
    left_indexer = xr.DataArray(left, dims=("_target_vertical",))
    right_indexer = xr.DataArray(right, dims=("_target_vertical",))

    source = prepared.rename({"_vertical": "_source_vertical"})
    if "bins" in source.coords:
        source = source.drop_vars("bins")

    sigma0 = source.isel(_source_vertical=left_indexer)
    sigma1 = source.isel(_source_vertical=right_indexer)
    sigma0 = sigma0.assign_coords(_target_vertical=target_dim)
    sigma1 = sigma1.assign_coords(_target_vertical=target_dim)

    weight0 = xr.DataArray(w0, dims=("_target_vertical",))
    weight1 = xr.DataArray(w1, dims=("_target_vertical",))

    out = np.sqrt((weight0 * sigma0) ** 2 + (weight1 * sigma1) ** 2)
    out = out.where(~xr.DataArray(outside, dims=("_target_vertical",)))
    out = out.rename({"_target_vertical": "bins"})
    out = out.assign_coords(bins=np.arange(out.sizes["bins"], dtype=float))
    return out


def _interpolate_dataset(
    dataset: dict[str, Any],
    *,
    source_vertical: xr.DataArray,
    target_vertical: np.ndarray,
    vertical_scale: str,
    location: str,
) -> None:
    dataset["signal"] = _interp_value(
        dataset.get("signal"), source_vertical, target_vertical,
        location=f"{location} signal",
    )
    dataset["error"] = _interp_error(
        dataset.get("error"), source_vertical, target_vertical,
        location=f"{location} error",
    )

    molecular = dataset.get("molecular", {})
    molecular["profile"] = _interp_value(
        molecular.get("profile"), source_vertical, target_vertical,
        location=f"{location} molecular profile",
    )

    metadata = dataset.get("metadata", {})
    old_metadata = dict(metadata)
    for key in VERTICAL_METADATA_KEYS:
        value = old_metadata.get(key)
        if key == vertical_scale:
            metadata[key] = xr.DataArray(
                target_vertical,
                dims=("bins",),
                coords={"bins": np.arange(target_vertical.size, dtype=float)},
                name=getattr(value, "name", None),
            )
        else:
            metadata[key] = _interp_value(
                value,
                source_vertical,
                target_vertical,
                location=f"{location} metadata {key!r}",
            )


def _finite_monotonic_centers(vertical: xr.DataArray, *, location: str) -> np.ndarray:
    centers = np.asarray(vertical.values, dtype=float)
    if centers.ndim != 1:
        raise ValueError(f"{location}: vertical coordinate must be one-dimensional")

    centers = centers[np.isfinite(centers)]
    if centers.size < 2:
        raise ValueError(f"{location}: at least two finite vertical points are required")

    diffs = np.diff(centers)
    if np.all(diffs > 0):
        return centers
    if np.all(diffs < 0):
        return centers[::-1]
    raise ValueError(f"{location}: vertical coordinate must be strictly monotonic")


def _source_bin_edges(centers: np.ndarray, *, location: str) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    if centers.size < 2:
        raise ValueError(f"{location}: at least two vertical points are required")

    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def _nominal_native_step(vertical: xr.DataArray, *, location: str) -> float:
    """Return a robust nominal native step for one selected vertical scale."""
    centers = _finite_monotonic_centers(vertical, location=location)
    steps = np.diff(centers)
    steps = steps[np.isfinite(steps) & (steps > 0)]
    if steps.size == 0:
        raise ValueError(f"{location}: no positive native vertical steps are available")
    return float(np.median(steps))


def _automatic_group_bin_width(
    group_kind: str,
    group_id: str,
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
) -> tuple[float, dict[str, float]]:
    steps: dict[str, float] = {}
    for entry_id, dataset in group.get("entries", {}).items():
        vertical = _as_1d_vertical(
            dataset.get("metadata", {}).get(vertical_scale),
            location=f"[{group_kind}:{group_id}] entry {entry_id!r} {vertical_scale!r}",
        )
        steps[entry_id] = _nominal_native_step(
            vertical,
            location=f"[{group_kind}:{group_id}] entry {entry_id!r} {vertical_scale!r}",
        )

    if not steps:
        raise ValueError(f"[{group_kind}:{group_id}] contains no participating entrys")
    return max(steps.values()), steps


def _group_source_extent(
    group_kind: str,
    group_id: str,
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
) -> tuple[float, float]:
    lows = []
    highs = []
    for entry_id, dataset in group.get("entries", {}).items():
        vertical = _as_1d_vertical(
            dataset.get("metadata", {}).get(vertical_scale),
            location=f"[{group_kind}:{group_id}] entry {entry_id!r} {vertical_scale!r}",
        )
        centers = _finite_monotonic_centers(
            vertical,
            location=f"[{group_kind}:{group_id}] entry {entry_id!r} {vertical_scale!r}",
        )
        edges = _source_bin_edges(
            centers,
            location=f"[{group_kind}:{group_id}] entry {entry_id!r} {vertical_scale!r}",
        )
        lows.append(float(edges[0]))
        highs.append(float(edges[-1]))

    if not lows:
        raise ValueError(f"[{group_kind}:{group_id}] contains no participating entrys")
    return min(lows), max(highs)


def _target_bin_edges(
    *,
    grid_start: float,
    width: float,
    source_high: float,
    vertical_min: Optional[float],
    vertical_max: Optional[float],
) -> np.ndarray:
    if width <= 0:
        raise ValueError("vertical_bin_width must be greater than zero")

    low_requested = grid_start if vertical_min is None else max(grid_start, vertical_min)
    start_index = int(np.ceil((low_requested - grid_start) / width - 1.0e-12))
    low_edge = grid_start + start_index * width

    if vertical_max is None:
        high_requested = source_high
        n_bins = int(np.ceil((high_requested - low_edge) / width - 1.0e-12))
    else:
        high_requested = vertical_max
        n_bins = int(np.floor((high_requested - low_edge) / width + 1.0e-12))

    if n_bins < 1:
        raise ValueError(
            "The requested vertical limits contain fewer than one complete target bin"
        )

    return low_edge + np.arange(n_bins + 1, dtype=float) * width


def _conservative_overlap_weights(
    vertical: xr.DataArray,
    target_edges: np.ndarray,
    *,
    location: str,
) -> tuple[xr.DataArray, bool]:
    """Build source-to-target physical overlap lengths once per dataset."""
    vertical = _as_1d_vertical(vertical, location=location)
    raw_centers = np.asarray(vertical.values, dtype=float)
    if not np.all(np.isfinite(raw_centers)):
        raise ValueError(f"{location}: vertical coordinate contains non-finite values")

    diffs = np.diff(raw_centers)
    descending = bool(np.all(diffs < 0))
    if not (np.all(diffs > 0) or descending):
        raise ValueError(f"{location}: vertical coordinate must be strictly monotonic")

    centers = raw_centers[::-1] if descending else raw_centers
    source_edges = _source_bin_edges(centers, location=location)

    source_left = source_edges[:-1][None, :]
    source_right = source_edges[1:][None, :]
    target_left = np.asarray(target_edges[:-1], dtype=float)[:, None]
    target_right = np.asarray(target_edges[1:], dtype=float)[:, None]

    overlap = np.minimum(source_right, target_right) - np.maximum(source_left, target_left)
    overlap = np.clip(overlap, 0.0, None)

    weights = xr.DataArray(
        overlap,
        dims=("target_bins", "source_bins"),
        coords={
            "target_bins": np.arange(overlap.shape[0], dtype=int),
            "source_bins": np.arange(overlap.shape[1], dtype=int),
        },
    )
    return weights, descending


def _prepare_value_for_binning(
    value: xr.DataArray,
    *,
    descending: bool,
    location: str,
) -> xr.DataArray:
    if "bins" not in value.dims:
        raise ValueError(f"{location}: value has no 'bins' dimension")
    if descending:
        value = value.isel(bins=slice(None, None, -1))
    if "bins" in value.coords:
        value = value.drop_vars("bins")
    value = value.rename({"bins": "source_bins"})
    return value.assign_coords(source_bins=np.arange(value.sizes["source_bins"], dtype=int))


def _bin_value_conservative(
    value: Any,
    weights: xr.DataArray,
    *,
    descending: bool,
    error: bool,
    location: str,
) -> Any:
    if not isinstance(value, xr.DataArray) or "bins" not in value.dims:
        return value

    prepared = _prepare_value_for_binning(
        value,
        descending=descending,
        location=location,
    )
    if prepared.sizes["source_bins"] != weights.sizes["source_bins"]:
        raise ValueError(
            f"{location}: value has {prepared.sizes['source_bins']} bins but the "
            f"vertical overlap map has {weights.sizes['source_bins']}"
        )

    valid = prepared.notnull()
    effective_weights = weights.where(valid, 0.0)
    denominator = effective_weights.sum("source_bins")

    if error:
        numerator = np.sqrt(((prepared.fillna(0.0) * weights) ** 2).sum("source_bins"))
    else:
        numerator = (prepared.fillna(0.0) * weights).sum("source_bins")

    out = (numerator / denominator).where(denominator > 0)
    out = out.rename({"target_bins": "bins"})
    return out.assign_coords(bins=np.arange(out.sizes["bins"], dtype=float))


def _bin_dataset(
    dataset: dict[str, Any],
    *,
    source_vertical: xr.DataArray,
    edges: np.ndarray,
    vertical_scale: str,
    location: str,
) -> None:
    centers = 0.5 * (edges[:-1] + edges[1:])
    weights, descending = _conservative_overlap_weights(
        source_vertical,
        edges,
        location=f"{location} {vertical_scale!r}",
    )

    dataset["signal"] = _bin_value_conservative(
        dataset.get("signal"), weights,
        descending=descending, error=False, location=f"{location} signal",
    )
    dataset["error"] = _bin_value_conservative(
        dataset.get("error"), weights,
        descending=descending, error=True, location=f"{location} error",
    )

    molecular = dataset.get("molecular", {})
    molecular["profile"] = _bin_value_conservative(
        molecular.get("profile"), weights,
        descending=descending, error=False, location=f"{location} molecular profile",
    )

    metadata = dataset.get("metadata", {})
    old_metadata = dict(metadata)
    for key in VERTICAL_METADATA_KEYS:
        value = old_metadata.get(key)
        if key == vertical_scale:
            metadata[key] = xr.DataArray(
                centers,
                dims=("bins",),
                coords={"bins": np.arange(centers.size, dtype=float)},
                name=getattr(value, "name", None),
            )
        else:
            metadata[key] = _bin_value_conservative(
                value, weights,
                descending=descending,
                error=False,
                location=f"{location} metadata {key!r}",
            )


def _combine_group_arrays(group_kind: str, group_id: str, group: dict[str, Any]) -> None:
    """Create aligned group arrays with dimensions ``(entry, bins)``."""
    entries = group.get("entries", {})
    if not entries:
        return

    vertical_scale = group.get("vertical_scale")
    vertical_grid = group.get("vertical_grid")
    if not isinstance(vertical_grid, xr.DataArray):
        raise ValueError(f"[{group_kind}:{group_id}] vertical_grid is missing after binning")

    def _combine(key: str) -> Optional[xr.DataArray]:
        arrays = []
        labels = []
        for entry_id, dataset in entries.items():
            value = dataset.get(key)
            if not isinstance(value, xr.DataArray) or "bins" not in value.dims:
                continue
            arrays.append(value)
            labels.append(entry_id)
        if not arrays:
            return None

        combined = xr.concat(arrays, dim=xr.IndexVariable("entry", labels), join="exact")
        z = np.asarray(vertical_grid.values)
        if vertical_scale == "bins":
            combined = combined.assign_coords(bins=z)
        else:
            combined = combined.assign_coords({vertical_scale: ("bins", z)})
        return combined

    signals = _combine("signal")
    errors = _combine("error")
    if signals is not None:
        group["signals"] = signals
    if errors is not None:
        group["errors"] = errors

    molecular_arrays = []
    molecular_labels = []
    for entry_id, dataset in entries.items():
        value = dataset.get("molecular", {}).get("profile")
        if not isinstance(value, xr.DataArray) or "bins" not in value.dims:
            continue
        molecular_arrays.append(value)
        molecular_labels.append(entry_id)

    if molecular_arrays:
        molecular = xr.concat(
            molecular_arrays,
            dim=xr.IndexVariable("entry", molecular_labels),
            join="exact",
        )
        z = np.asarray(vertical_grid.values)
        if vertical_scale == "bins":
            molecular = molecular.assign_coords(bins=z)
        else:
            molecular = molecular.assign_coords({vertical_scale: ("bins", z)})
        group["molecular_profiles"] = molecular


def _harmonize_group_interpolation(
    group_kind: str,
    group_id: str,
    group: dict[str, Any],
    *,
    vertical_scale: str,
    vertical_min: Optional[float],
    vertical_max: Optional[float],
) -> None:
    """Interpolate onto the coarsest native grid in the group.

    The target grid is taken from the participating entry with the largest
    nominal native step on ``vertical_scale``. This deliberately avoids
    upscaling any coarser dataset onto a finer grid. The common vertical
    overlap and optional ``vertical_min`` / ``vertical_max`` limits are still
    applied before interpolation.
    """
    entries = group.get("entries", {})
    if not entries:
        raise ValueError(f"[{group_kind}:{group_id}] contains no participating entrys")

    _, native_steps = _automatic_group_bin_width(
        group_kind,
        group_id,
        group,
        vertical_scale=vertical_scale,
    )
    worst_step = max(native_steps.values())

    # Keep entry order deterministic when two entries have the same
    # coarsest nominal step.
    tolerance = max(1.0e-12, abs(worst_step) * 1.0e-9)
    target_entry_id = next(
        entry_id
        for entry_id in entries
        if abs(native_steps[entry_id] - worst_step) <= tolerance
    )
    target_dataset = entries[target_entry_id]

    low, high = _common_bounds(
        group_kind, group_id, group,
        vertical_scale=vertical_scale,
        vertical_min=vertical_min,
        vertical_max=vertical_max,
    )

    target_native_vertical = _as_1d_vertical(
        target_dataset.get("metadata", {}).get(vertical_scale),
        location=(
            f"[{group_kind}:{group_id}] interpolation target entry "
            f"{target_entry_id!r} {vertical_scale!r}"
        ),
    )
    z_target = np.asarray(target_native_vertical.values, dtype=float)
    target_vertical = z_target[
        np.isfinite(z_target) & (z_target >= low) & (z_target <= high)
    ]
    target_vertical = np.unique(target_vertical)
    if target_vertical.size < 2:
        raise ValueError(
            f"[{group_kind}:{group_id}] coarsest native grid from entry "
            f"{target_entry_id!r} contains fewer than two points inside the "
            "common vertical overlap"
        )

    for entry_id, dataset in entries.items():
        source_vertical = _as_1d_vertical(
            dataset.get("metadata", {}).get(vertical_scale),
            location=(
                f"[{group_kind}:{group_id}] entry {entry_id!r} "
                f"{vertical_scale!r}"
            ),
        )
        _interpolate_dataset(
            dataset,
            source_vertical=source_vertical,
            target_vertical=target_vertical,
            vertical_scale=vertical_scale,
            location=f"[{group_kind}:{group_id}] entry {entry_id!r}",
        )

    group["vertical_grid"] = xr.DataArray(
        target_vertical,
        dims=("bins",),
        coords={"bins": np.arange(target_vertical.size, dtype=float)},
        name=vertical_scale,
    )
    group["vertical_scale"] = vertical_scale
    group["vertical_method"] = "interpolation"
    group["interpolation_target_entry"] = target_entry_id
    group["interpolation_target_step"] = _from_internal_units(
        worst_step, vertical_scale
    )
    group["native_vertical_steps"] = {
        entry_id: _from_internal_units(step, vertical_scale)
        for entry_id, step in native_steps.items()
    }



def _harmonize_group_binning(
    group_kind: str,
    group_id: str,
    group: dict[str, Any],
    *,
    general: Mapping[str, Any],
    vertical_scale: str,
) -> None:
    requested_width = group.get("vertical_bin_width")
    if requested_width is not None:
        width_source = "group"
    else:
        requested_width = general.get("vertical_bin_width")
        width_source = "general"

    # Always determine the native steps. They define the minimum sensible
    # conservative-bin width for this group.
    worst_step_internal, native_steps = _automatic_group_bin_width(
        group_kind,
        group_id,
        group,
        vertical_scale=vertical_scale,
    )

    if requested_width is None:
        width_internal = worst_step_internal
        width_source = "automatic coarsest native step"
    else:
        width_internal = _to_internal_units(requested_width, vertical_scale)
        if width_internal is None or width_internal <= 0:
            raise ValueError(f"[{group_kind}:{group_id}] vertical_bin_width must be > 0")

        # Do not create output bins finer than the coarsest native sampling in
        # the group. Such a grid would imply vertical resolution that at least
        # one participating entry does not actually provide.
        tolerance = max(1.0e-12, abs(worst_step_internal) * 1.0e-9)
        if width_internal + tolerance < worst_step_internal:
            requested_display = _from_internal_units(width_internal, vertical_scale)
            worst_display = _from_internal_units(worst_step_internal, vertical_scale)
            print(
                f"      Warning: [{group_kind}:{group_id}] requested vertical_bin_width "
                f"{requested_display:.6g} {_units_label(vertical_scale)} is smaller than "
                f"the coarsest native step {worst_display:.6g} "
                f"{_units_label(vertical_scale)}. Using the coarsest native step instead."
            )
            width_internal = worst_step_internal
            width_source = "coarsest native step fallback"

    grid_start = _to_internal_units(general.get("first_bin_left_edge", 0.0), vertical_scale)
    vertical_min = _to_internal_units(general.get("vertical_min"), vertical_scale)
    vertical_max = _to_internal_units(general.get("vertical_max"), vertical_scale)
    if grid_start is None:
        grid_start = 0.0

    _, source_high = _group_source_extent(
        group_kind,
        group_id,
        group,
        vertical_scale=vertical_scale,
    )
    edges = _target_bin_edges(
        grid_start=grid_start,
        width=width_internal,
        source_high=source_high,
        vertical_min=vertical_min,
        vertical_max=vertical_max,
    )

    before = {}
    for entry_id, dataset in group.get("entries", {}).items():
        signal = dataset.get("signal")
        if isinstance(signal, xr.DataArray):
            before[entry_id] = signal.sizes.get("bins")

        source_vertical = _as_1d_vertical(
            dataset.get("metadata", {}).get(vertical_scale),
            location=f"[{group_kind}:{group_id}] entry {entry_id!r} {vertical_scale!r}",
        )
        _bin_dataset(
            dataset,
            source_vertical=source_vertical,
            edges=edges,
            vertical_scale=vertical_scale,
            location=f"[{group_kind}:{group_id}] entry {entry_id!r}",
        )

    # The selected vertical scale is represented by the geometric center
    # of every conservative destination bin.
    centers = 0.5 * (edges[:-1] + edges[1:])
    group["vertical_grid"] = xr.DataArray(
        centers,
        dims=("bins",),
        coords={"bins": np.arange(centers.size, dtype=float)},
        name=vertical_scale,
    )
    group["vertical_scale"] = vertical_scale
    group["vertical_bin_width"] = _from_internal_units(width_internal, vertical_scale)
    group["vertical_bin_width_source"] = width_source
    if native_steps:
        group["native_vertical_steps"] = {
            entry_id: _from_internal_units(step, vertical_scale)
            for entry_id, step in native_steps.items()
        }

    _combine_group_arrays(group_kind, group_id, group)

    print(
        f"    - {group_id}: width={group['vertical_bin_width']:.6g} "
        f"{_units_label(vertical_scale)} ({width_source}), "
        f"grid={centers.size} bins"
    )
    if native_steps:
        for entry_id, step in group["native_vertical_steps"].items():
            print(
                f"        {entry_id}: native step={step:.6g} "
                f"{_units_label(vertical_scale)}; {before.get(entry_id, '?')} -> {centers.size} bins"
            )
    else:
        for entry_id in group.get("entries", {}):
            print(f"        {entry_id}: {before.get(entry_id, '?')} -> {centers.size} bins")


def harmonize_intercomparison_vertical(
    intercomparison_info: Mapping[str, Any],
    intercomparison_bundles: Mapping[str, Any],
) -> dict[str, Any]:
    """Harmonize active groups using interpolation or conservative binning."""
    out = _copy_nested(intercomparison_bundles)
    general = intercomparison_info["general"]
    vertical_scale = general["vertical_scale"]
    vertical_method = general["vertical_method"]

    if vertical_scale not in VERTICAL_SCALES:
        raise ValueError(
            f"Unsupported vertical_scale {vertical_scale!r}; expected one of {VERTICAL_SCALES}"
        )
    if vertical_method not in VERTICAL_METHODS:
        raise ValueError(
            f"Unsupported vertical_method {vertical_method!r}; expected one of {VERTICAL_METHODS}"
        )

    vertical_min = _to_internal_units(general.get("vertical_min"), vertical_scale)
    vertical_max = _to_internal_units(general.get("vertical_max"), vertical_scale)

    print()
    print("-----------------------------------------------")
    print("Vertical harmonization")
    print("-----------------------------------------------")
    print(f"Vertical scale: {vertical_scale}")
    print(f"Method: {vertical_method}")

    if vertical_method == "vertical_binning":
        print(
            f"First bin left edge: {general.get('first_bin_left_edge', 0.0)} "
            f"{_units_label(vertical_scale)}"
        )
        if general.get("vertical_bin_width") is None:
            print("General bin width: automatic coarsest native step unless overridden by a group")
        else:
            print(
                f"General bin width: {general['vertical_bin_width']} "
                f"{_units_label(vertical_scale)}"
            )

    if general.get("vertical_min") is not None or general.get("vertical_max") is not None:
        print(
            f"Retained limits: {general.get('vertical_min')} to "
            f"{general.get('vertical_max')} {_units_label(vertical_scale)}"
        )

    for group_kind, collection_key in (
        ("channel_group", "channel_groups"),
        ("pair_group", "pair_groups"),
    ):
        groups = out.get(collection_key, {})
        title = "Channel groups" if group_kind == "channel_group" else "Pair groups"
        if not groups:
            print(f"{title}: none")
            continue

        print(f"{title}:")
        for group_id, group in groups.items():
            if vertical_method == "interpolation":
                _harmonize_group_interpolation(
                    group_kind,
                    group_id,
                    group,
                    vertical_scale=vertical_scale,
                    vertical_min=vertical_min,
                    vertical_max=vertical_max,
                )
                _combine_group_arrays(group_kind, group_id, group)
                print(
                    f"    - {group_id}: interpolated to coarsest native grid "
                    f"from {group['interpolation_target_entry']!r} "
                    f"(step={group['interpolation_target_step']:.6g} "
                    f"{_units_label(vertical_scale)}, "
                    f"{group['vertical_grid'].sizes['bins']} bins)"
                )
                for entry_id, step in group.get("native_vertical_steps", {}).items():
                    print(
                        f"        {entry_id}: native step={step:.6g} "
                        f"{_units_label(vertical_scale)}"
                    )
            else:
                _harmonize_group_binning(
                    group_kind,
                    group_id,
                    group,
                    general=general,
                    vertical_scale=vertical_scale,
                )
                group["vertical_method"] = "vertical_binning"

    print("Vertical harmonization complete.")
    return out
