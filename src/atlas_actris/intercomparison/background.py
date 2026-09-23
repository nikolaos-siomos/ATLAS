#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Background correction for ATLAS intercomparison bundles.

This module is deliberately separate from temporal filtering and normalization.
For channel and pair groups it follows the existing ATLAS concept: estimate one
mean background over a configured far-range region and subtract that scalar from
the profile/ratio. The signal uncertainty is left unchanged, matching the
existing ATLAS background-correction step. Pair range metadata are prepared from
the pair's reference channel during bundle construction.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import xarray as xr


def _copy_nested(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _copy_nested(item) for key, item in value.items()}
    return value


def _as_scalar(value: Any, *, location: str) -> float:
    if isinstance(value, xr.DataArray):
        value = value.squeeze(drop=True).values
    arr = np.asarray(value)
    if arr.size != 1:
        raise ValueError(f"{location}: expected one scalar value, got shape {arr.shape}")
    return float(arr.reshape(-1)[0])


def _reference_background_region_from_metadata(
    group_id: str,
    group: Mapping[str, Any],
) -> list[float] | None:
    """Derive a physical range [km] from reference channel background bins."""

    reference_entry = group.get("reference_entry")
    reference = group.get("entries", {}).get(reference_entry)
    if reference is None:
        return None

    metadata = reference.get("metadata", {})
    channel_info = metadata.get("channel_info")
    ranges = metadata.get("range")
    bins = metadata.get("bins")
    if channel_info is None or ranges is None or bins is None:
        return None

    try:
        low_bin = _as_scalar(
            channel_info.sel(parameters="background_low_bin"),
            location=f"[channel_group:{group_id}] reference background_low_bin",
        )
        high_bin = _as_scalar(
            channel_info.sel(parameters="background_high_bin"),
            location=f"[channel_group:{group_id}] reference background_high_bin",
        )
    except (KeyError, ValueError):
        return None

    try:
        mask = (bins >= low_bin) & (bins <= high_bin)
        selected = ranges.where(mask, drop=True)
        if selected.sizes.get("bins", 0) == 0:
            return None
        return [
            float(selected.min(skipna=True).values) / 1000.0,
            float(selected.max(skipna=True).values) / 1000.0,
        ]
    except Exception:
        return None


def _resolve_background_region(
    group_kind: str,
    group_id: str,
    group: Mapping[str, Any],
) -> list[float]:
    """Return the already resolved group background interval in km."""

    region = group.get("background_region")
    if region in (None, []):
        raise ValueError(
            f"[{group_kind}:{group_id}] background_correction=True but "
            "background_region is empty after resolving [general] defaults."
        )
    return [float(region[0]), float(region[1])]


def _background_correct_dataset(
    dataset: dict[str, Any],
    *,
    region_km: list[float],
    vertical_scale: str,
    location: str,
) -> None:
    signal = dataset.get("signal")
    vertical = dataset.get("metadata", {}).get(vertical_scale)

    if not isinstance(signal, xr.DataArray):
        raise ValueError(f"{location}: signal is not an xarray.DataArray")
    if "bins" not in signal.dims:
        raise ValueError(f"{location}: signal does not contain a 'bins' dimension")
    if not isinstance(vertical, xr.DataArray):
        raise ValueError(
            f"{location}: {vertical_scale!r} vertical metadata are required "
            "for background correction"
        )

    low_m = 1000.0 * float(region_km[0])
    high_m = 1000.0 * float(region_km[1])
    mask = (vertical >= low_m) & (vertical <= high_m)

    background_signal = signal.where(mask)
    n_bins = background_signal.notnull().sum("bins")
    background = background_signal.mean("bins", skipna=True)

    if bool(np.asarray((n_bins == 0).all().values)):
        raise ValueError(
            f"{location}: background region {region_km} km contains no valid bins"
        )

    dataset["background"] = background
    dataset["signal"] = signal - background
    dataset["background_region"] = list(region_km)


def apply_intercomparison_background_correction(
    intercomparison_info: Mapping[str, Any],
    intercomparison_bundles: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply configured background correction as a standalone workflow step."""

    out = _copy_nested(intercomparison_bundles)
    vertical_scale = intercomparison_info["general"]["vertical_scale"]

    print()
    print("-----------------------------------------------")
    print("Background correction")
    print("-----------------------------------------------")
    print(f"Vertical scale: {vertical_scale}")

    channel_groups = out.get("channel_groups", {})
    if not channel_groups:
        print("Channel groups: none")
    else:
        print("Channel groups:")
        for group_id, group in channel_groups.items():
            if not group.get("background_correction", False):
                print(f"    - {group_id}: disabled -> unchanged")
                continue

            region = _resolve_background_region("channel_group", group_id, group)
            print(f"    - {group_id}: enabled, region = {region} km")
            for dataset_id, dataset in group.get("entries", {}).items():
                _background_correct_dataset(
                    dataset,
                    region_km=region,
                    vertical_scale=vertical_scale,
                    location=f"[channel_group:{group_id}] entry {dataset_id!r}",
                )
                print(f"        {dataset_id}: background estimated and subtracted")

    pair_groups = out.get("pair_groups", {})
    if pair_groups:
        print("Pair groups:")
        for group_id, group in pair_groups.items():
            if not group.get("background_correction", False):
                print(f"    - {group_id}: disabled -> unchanged")
                continue

            region = _resolve_background_region("pair_group", group_id, group)
            print(f"    - {group_id}: enabled, region = {region} km")
            for dataset_id, dataset in group.get("entries", {}).items():
                _background_correct_dataset(
                    dataset,
                    region_km=region,
                    vertical_scale=vertical_scale,
                    location=f"[pair_group:{group_id}] entry {dataset_id!r}",
                )
                print(f"        {dataset_id}: background estimated and subtracted")
    else:
        print("Pair groups: none")

    print("Background correction complete.")
    print()
    return out
