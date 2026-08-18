#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Normalization for ATLAS intercomparison bundles.

Every active group is normalized to an absolute reference target over the
configured vertical interval; profiles are not normalized to one.

Channel groups:
- normalise_to_molecular=False -> reference measured channel profile
- normalise_to_molecular=True  -> reference molecular channel profile

Pair groups:
- always normalize to the reference measured pair ratio
- molecular_ratio is retained only as an optional plotting reference
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import xarray as xr


def _copy_nested(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _copy_nested(item) for key, item in value.items()}
    return value


def _collapse_time(value: Any) -> Any:
    if not isinstance(value, xr.DataArray) or "time" not in value.dims:
        return value
    if value.sizes["time"] == 1:
        return value.squeeze("time", drop=True)
    return value.mean("time", skipna=True)


def _region_mean(
    value: xr.DataArray,
    vertical: xr.DataArray,
    region_km: list[float],
    *,
    location: str,
) -> xr.DataArray:
    if "bins" not in value.dims:
        raise ValueError(f"{location}: value does not contain a 'bins' dimension")

    low_m = 1000.0 * float(region_km[0])
    high_m = 1000.0 * float(region_km[1])
    mask = (vertical >= low_m) & (vertical <= high_m)
    selected = value.where(mask)
    n_bins = selected.notnull().sum("bins")
    if bool(np.asarray((n_bins == 0).all().values)):
        raise ValueError(
            f"{location}: normalization region {region_km} km contains no valid bins"
        )
    return selected.mean("bins", skipna=True)


def _resolved_region(group_kind: str, group_id: str, group: Mapping[str, Any]) -> list[float]:
    region = group.get("normalisation_region")
    if region in (None, []):
        raise ValueError(
            f"[{group_kind}:{group_id}] normalisation=True but normalisation_region "
            "is empty after resolving [general] defaults."
        )
    return [float(region[0]), float(region[1])]


def _reference_dataset(group_kind: str, group_id: str, group: Mapping[str, Any]) -> dict[str, Any]:
    reference_id = group.get("reference_dataset")
    reference = group.get("datasets", {}).get(reference_id)
    if reference is None:
        raise ValueError(
            f"[{group_kind}:{group_id}] reference dataset {reference_id!r} "
            "does not participate in this active group."
        )
    return reference


def _measured_reference_target(
    group_kind: str,
    group_id: str,
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
    region: list[float],
) -> tuple[xr.DataArray, str]:
    reference = _reference_dataset(group_kind, group_id, group)
    signal = reference.get("signal")
    vertical = reference.get("metadata", {}).get(vertical_scale)

    if not isinstance(signal, xr.DataArray):
        raise ValueError(f"[{group_kind}:{group_id}] reference signal is missing")
    if not isinstance(vertical, xr.DataArray):
        raise ValueError(
            f"[{group_kind}:{group_id}] reference dataset has no "
            f"{vertical_scale!r} vertical metadata"
        )

    target = _region_mean(
        signal,
        vertical,
        region,
        location=f"[{group_kind}:{group_id}] reference measured signal",
    )
    return target, "reference measured signal"


def _channel_target(
    group_id: str,
    group: Mapping[str, Any],
    *,
    vertical_scale: str,
    region: list[float],
) -> tuple[xr.DataArray, str]:
    if not bool(group.get("normalise_to_molecular", False)):
        return _measured_reference_target(
            "channel_group", group_id, group,
            vertical_scale=vertical_scale,
            region=region,
        )

    reference = _reference_dataset("channel_group", group_id, group)
    molecular = _collapse_time(reference.get("molecular", {}).get("profile"))
    vertical = reference.get("metadata", {}).get(vertical_scale)

    if not isinstance(molecular, xr.DataArray):
        source = reference.get("molecular", {}).get("profile_source", "molecular")
        raise ValueError(
            f"[channel_group:{group_id}] normalise_to_molecular=True but the "
            f"reference dataset has no molecular profile from source {source!r}."
        )
    if not isinstance(vertical, xr.DataArray):
        raise ValueError(
            f"[channel_group:{group_id}] reference dataset has no "
            f"{vertical_scale!r} vertical metadata"
        )

    target = _region_mean(
        molecular,
        vertical,
        region,
        location=f"[channel_group:{group_id}] reference molecular profile",
    )
    return target, "reference molecular profile"


def _apply_target_to_group(
    group_kind: str,
    group_id: str,
    group: dict[str, Any],
    *,
    vertical_scale: str,
    region: list[float],
    target_mean: xr.DataArray,
    target_description: str,
) -> None:
    group["normalisation_region"] = region
    group["normalisation_target"] = target_description

    for dataset_id, dataset in group.get("datasets", {}).items():
        signal = dataset.get("signal")
        error = dataset.get("error")
        vertical = dataset.get("metadata", {}).get(vertical_scale)

        if not isinstance(signal, xr.DataArray):
            raise ValueError(f"[{group_kind}:{group_id}] dataset {dataset_id!r}: signal missing")
        if not isinstance(vertical, xr.DataArray):
            raise ValueError(
                f"[{group_kind}:{group_id}] dataset {dataset_id!r}: "
                f"{vertical_scale!r} vertical metadata missing"
            )

        source_mean = _region_mean(
            signal,
            vertical,
            region,
            location=f"[{group_kind}:{group_id}] dataset {dataset_id!r} signal",
        )
        factor = target_mean / source_mean
        factor = factor.where(np.isfinite(factor) & (source_mean != 0))

        dataset["normalisation_factor"] = factor
        dataset["signal"] = signal * factor
        if isinstance(error, xr.DataArray):
            dataset["error"] = error * abs(factor)


def apply_intercomparison_normalization(
    intercomparison_info: Mapping[str, Any],
    intercomparison_bundles: Mapping[str, Any],
) -> dict[str, Any]:
    """Normalize channel and pair groups to their configured reference target."""

    out = _copy_nested(intercomparison_bundles)
    vertical_scale = intercomparison_info["general"]["vertical_scale"]

    print()
    print("-----------------------------------------------")
    print("Normalization")
    print("-----------------------------------------------")
    print(f"Vertical scale: {vertical_scale}")

    channel_groups = out.get("channel_groups", {})
    if not channel_groups:
        print("Channel groups: none")
    else:
        print("Channel groups:")
        for group_id, group in channel_groups.items():
            if not group.get("normalisation", False):
                print(f"    - {group_id}: disabled -> unchanged")
                continue

            region = _resolved_region("channel_group", group_id, group)
            target_mean, target_description = _channel_target(
                group_id,
                group,
                vertical_scale=vertical_scale,
                region=region,
            )
            _apply_target_to_group(
                "channel_group",
                group_id,
                group,
                vertical_scale=vertical_scale,
                region=region,
                target_mean=target_mean,
                target_description=target_description,
            )
            print(
                f"    - {group_id}: normalized over {region} km "
                f"to {target_description} (absolute reference level, not 1)"
            )
            for dataset_id in group.get("datasets", {}):
                print(f"        {dataset_id}: normalization factor calculated and applied")

    pair_groups = out.get("pair_groups", {})
    if not pair_groups:
        print("Pair groups: none")
    else:
        print("Pair groups:")
        for group_id, group in pair_groups.items():
            if not group.get("normalisation", False):
                print(f"    - {group_id}: disabled -> unchanged")
                continue

            region = _resolved_region("pair_group", group_id, group)
            target_mean, target_description = _measured_reference_target(
                "pair_group",
                group_id,
                group,
                vertical_scale=vertical_scale,
                region=region,
            )
            _apply_target_to_group(
                "pair_group",
                group_id,
                group,
                vertical_scale=vertical_scale,
                region=region,
                target_mean=target_mean,
                target_description=target_description,
            )
            print(
                f"    - {group_id}: normalized over {region} km "
                "to reference measured pair ratio (absolute reference level, not 1)"
            )
            for dataset_id in group.get("datasets", {}):
                print(f"        {dataset_id}: normalization factor calculated and applied")

    print("Normalization complete.")
    print()
    return out
