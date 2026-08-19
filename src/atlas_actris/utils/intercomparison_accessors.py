#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Small group-wise accessors for ATLAS intercomparison bundles.

Each public function extracts one specific type of information from one
channel or pair group and returns a dictionary keyed by dataset ID.  The
helpers intentionally do not rebuild the complete bundle structure; they are
meant to make analysis and plotting code concise while keeping data grouped by
comparison group.
"""

from __future__ import annotations

from typing import Any, Mapping


VERTICAL_SCALES = ("bins", "range", "height_agl", "height_asl")


def _get_group(
    intercomparison_bundles: Mapping[str, Any],
    *,
    collection: str,
    group_id: str,
) -> Mapping[str, Any]:
    groups = intercomparison_bundles.get(collection, {})
    if group_id not in groups:
        available = ", ".join(groups.keys()) or "none"
        kind = "channel" if collection == "channel_groups" else "pair"
        raise KeyError(
            f"Unknown {kind} group {group_id!r}. Available groups: {available}"
        )
    return groups[group_id]


def _collect_dataset_entry(
    group: Mapping[str, Any],
    key: str,
) -> dict[str, Any]:
    return {
        dataset_id: dataset.get(key)
        for dataset_id, dataset in group.get("datasets", {}).items()
    }


def _collect_dataset_metadata(
    group: Mapping[str, Any],
    key: str,
) -> dict[str, Any]:
    return {
        dataset_id: dataset.get("metadata", {}).get(key)
        for dataset_id, dataset in group.get("datasets", {}).items()
    }


def _collect_molecular_profiles(group: Mapping[str, Any]) -> dict[str, Any]:
    return {
        dataset_id: dataset.get("molecular", {}).get("profile")
        for dataset_id, dataset in group.get("datasets", {}).items()
    }


def _validate_vertical_scale(vertical_scale: str) -> str:
    vertical_scale = str(vertical_scale).strip().lower()
    if vertical_scale not in VERTICAL_SCALES:
        raise ValueError(
            f"vertical_scale must be one of {VERTICAL_SCALES}, "
            f"got {vertical_scale!r}"
        )
    return vertical_scale


# -----------------------------------------------------------------------------
# Channel groups
# -----------------------------------------------------------------------------

def collect_channel_group_signals(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
) -> dict[str, Any]:
    """Return measured/processed signals for one channel group.

    Returns
    -------
    dict
        ``{dataset_id: signal}``
    """
    group = _get_group(
        intercomparison_bundles,
        collection="channel_groups",
        group_id=group_id,
    )
    return _collect_dataset_entry(group, "signal")


def collect_channel_group_errors(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
) -> dict[str, Any]:
    """Return signal uncertainties for one channel group."""
    group = _get_group(
        intercomparison_bundles,
        collection="channel_groups",
        group_id=group_id,
    )
    return _collect_dataset_entry(group, "error")


def collect_channel_group_vertical(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
    vertical_scale: str,
) -> dict[str, Any]:
    """Return one vertical coordinate for every dataset in a channel group.

    ``vertical_scale`` may be ``bins``, ``range``, ``height_agl``, or
    ``height_asl``.
    """
    vertical_scale = _validate_vertical_scale(vertical_scale)
    group = _get_group(
        intercomparison_bundles,
        collection="channel_groups",
        group_id=group_id,
    )
    return _collect_dataset_metadata(group, vertical_scale)


def collect_channel_group_molecular(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
) -> dict[str, Any]:
    """Return molecular attenuated-backscatter profiles for a channel group.

    The bundle loader has already selected ``opto_parameters='atten_bsc'``.
    """
    group = _get_group(
        intercomparison_bundles,
        collection="channel_groups",
        group_id=group_id,
    )
    return _collect_molecular_profiles(group)


# -----------------------------------------------------------------------------
# Pair groups
# -----------------------------------------------------------------------------

def collect_pair_group_signals(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
) -> dict[str, Any]:
    """Return measured/processed pair ratios (e.g. VLDR) for one pair group."""
    group = _get_group(
        intercomparison_bundles,
        collection="pair_groups",
        group_id=group_id,
    )
    return _collect_dataset_entry(group, "signal")


def collect_pair_group_errors(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
) -> dict[str, Any]:
    """Return pair-ratio uncertainties for one pair group."""
    group = _get_group(
        intercomparison_bundles,
        collection="pair_groups",
        group_id=group_id,
    )
    return _collect_dataset_entry(group, "error")


def collect_pair_group_vertical(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
    vertical_scale: str,
) -> dict[str, Any]:
    """Return one vertical coordinate for every dataset in a pair group."""
    vertical_scale = _validate_vertical_scale(vertical_scale)
    group = _get_group(
        intercomparison_bundles,
        collection="pair_groups",
        group_id=group_id,
    )
    return _collect_dataset_metadata(group, vertical_scale)


def collect_pair_group_molecular(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
) -> dict[str, Any]:
    """Return molecular-ratio profiles for one pair group."""
    group = _get_group(
        intercomparison_bundles,
        collection="pair_groups",
        group_id=group_id,
    )
    return _collect_molecular_profiles(group)


def collect_pair_group_molecular_info(
    intercomparison_bundles: Mapping[str, Any],
    group_id: str,
) -> dict[str, Any]:
    """Return molecular metadata (``molecular_info``) for one pair group."""
    group = _get_group(
        intercomparison_bundles,
        collection="pair_groups",
        group_id=group_id,
    )
    return {
        dataset_id: dataset.get("molecular", {}).get("metadata")
        for dataset_id, dataset in group.get("datasets", {}).items()
    }
