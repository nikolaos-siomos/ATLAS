#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare comparison-oriented channel and pair groups from dataset stages."""

from __future__ import annotations

from typing import Any, Mapping, Protocol


class StageStoreProtocol(Protocol):
    def get_entry(self, dataset_id: str, qa_test: str, parameter: str) -> Any:
        ...

    def get_optional_entry(self, dataset_id: str, qa_test: str, parameter: str) -> Any:
        ...


DATASET_METADATA_PARAMETERS = (
    "system_info",
    "radiosonde_info",
    "time_info",
)

CHANNEL_METADATA_PARAMETERS = (
    "channel_info",
    "bins",
    "range",
    "height_agl",
    "height_asl",
    "shots",
)

PAIR_METADATA_PARAMETERS = (
    "pol_cal_info",
)

PAIR_VERTICAL_PARAMETERS = (
    "bins",
    "range",
    "height_agl",
    "height_asl",
)

CHANNEL_MOLECULAR_PROFILE_SOURCE = "molecular"
PAIR_MOLECULAR_PROFILE_SOURCE = "molecular_ratio"
PAIR_MOLECULAR_METADATA_SOURCE = "molecular_info"


def _load_entry(
    stage_store: StageStoreProtocol,
    *,
    dataset_id: str,
    qa_test: str,
    parameter: str,
    config_field: str,
    group_kind: str,
    group_id: str,
) -> Any:
    """Load one source and translate importer failures into INI-oriented errors."""

    try:
        return stage_store.get_entry(dataset_id, qa_test, parameter)
    except (ValueError, KeyError, FileNotFoundError) as exc:
        detail = exc.args[0] if getattr(exc, "args", None) else str(exc)
        raise ValueError(
            f"Intercomparison configuration error while preparing "
            f"[{group_kind}:{group_id}].\n"
            f"Dataset: {dataset_id!r}\n"
            f"QA test: {qa_test!r}\n"
            f"Requested source: {parameter!r}\n"
            f"Configuration field: [{ 'dataset:' + dataset_id }] {config_field}\n\n"
            f"The requested source could not be loaded from this QA test. "
            f"Check that {config_field} points to a parameter exported under "
            f"qa_test={qa_test!r}.\n"
            f"Importer detail: {detail}"
        ) from None


def _available_ids(value: Any, coordinate: str, *, limit: int = 20) -> str:
    """Return a compact list of available product IDs for an error message."""

    try:
        raw = value.coords[coordinate].values
        values = [str(item) for item in raw.tolist()]
    except Exception:
        return ""

    if not values:
        return ""
    shown = values[:limit]
    suffix = " ..." if len(values) > limit else ""
    return ", ".join(shown) + suffix


def _select_product(value: Any, product_id: str, *, kind: str, location: str) -> Any:
    """Select one ATLAS channel or pair using the fixed coordinate name."""

    product_id = product_id.lower()

    if kind == "channel":
        coordinate = "channel"
    elif kind == "pair":
        coordinate = "pair"
    else:
        raise ValueError(f"Unsupported product kind: {kind!r}")

    try:
        return value.sel({coordinate: product_id})
    except Exception:
        available = _available_ids(value, coordinate)
        extra = (
            f" Available {kind} IDs: {available}."
            if available
            else ""
        )
        raise ValueError(
            f"{location}: {kind} ID {product_id!r} was not found in "
            f"the exported source along coordinate {coordinate!r}.{extra} "
            f"Check the corresponding atlas_{kind}_id value in the INI file."
        ) from None


def _load_optional_metadata(
    stage_store: StageStoreProtocol,
    *,
    dataset_id: str,
    qa_test: str,
    parameter: str,
) -> Any:
    """Load optional metadata under the dataset QA test; absent entries become None."""

    return stage_store.get_optional_entry(dataset_id, qa_test, parameter)


def _select_optional_product_metadata(
    value: Any,
    product_id: str,
    *,
    kind: str,
    location: str,
) -> Any:
    """Select one channel/pair from optional metadata, preserving None."""

    if value is None:
        return None
    return _select_product(value, product_id, kind=kind, location=location)



def _derive_molecular_pair_id(pair_id: str, *, location: str) -> str:
    """Derive the ATLAS molecular pair ID from an 8-character pair ID.

    Molecular pair IDs use the same identifier except that character 6
    (Python index 5) is always ``m``.
    """

    normalized = pair_id.strip().lower()
    if len(normalized) != 8:
        raise ValueError(
            f"{location}: atlas_pair_id {normalized!r} must contain exactly 8 "
            f"characters so the molecular pair ID can be derived."
        )
    return normalized[:5] + "m" + normalized[6:]


def _select_optional_molecular_pair_metadata(
    value: Any,
    molecular_pair_id: str,
    *,
    atlas_pair_id: str,
    parameter: str,
    group_id: str,
    dataset_id: str,
) -> Any:
    """Select a derived molecular pair ID from optional pair metadata."""

    if value is None:
        return None

    try:
        return value.sel(pair=molecular_pair_id)
    except Exception:
        available = _available_ids(value, "pair")
        extra = f" Available pair IDs: {available}." if available else ""
        raise ValueError(
            f"pair group {group_id!r}, dataset {dataset_id!r}, metadata {parameter!r}: "
            f"derived molecular pair ID {molecular_pair_id!r} was not found along "
            f"coordinate 'pair'. It was derived from atlas_pair_id {atlas_pair_id!r} "
            f"by replacing character 6 with 'm'.{extra}"
        ) from None



def _pair_reference_channel_id(
    pol_cal_info: Any,
    *,
    group_id: str,
    dataset_id: str,
) -> str | None:
    """Return the pair's reflected/reference channel ID from ``pol_cal_info``."""

    if pol_cal_info is None:
        return None
    try:
        value = pol_cal_info.sel(parameters="ch_r")
        raw = value.squeeze(drop=True).values
        channel_id = str(raw.item() if hasattr(raw, "item") else raw).strip().lower()
    except Exception:
        return None
    return channel_id or None


def _load_pair_vertical_metadata(
    stage_store: StageStoreProtocol,
    *,
    dataset_id: str,
    qa_test: str,
    channel_id: str | None,
    group_id: str,
) -> dict[str, Any]:
    """Load pair vertical coordinates from the pair's ``ch_r`` channel."""

    if not channel_id:
        return {parameter: None for parameter in PAIR_VERTICAL_PARAMETERS}

    out: dict[str, Any] = {}
    for parameter in PAIR_VERTICAL_PARAMETERS:
        entry = _load_optional_metadata(
            stage_store,
            dataset_id=dataset_id,
            qa_test=qa_test,
            parameter=parameter,
        )
        out[parameter] = _select_optional_product_metadata(
            entry,
            channel_id,
            kind="channel",
            location=(
                f"pair group {group_id!r}, dataset {dataset_id!r}, "
                f"vertical metadata {parameter!r} via ch_r={channel_id!r}"
            ),
        )
    return out

def _prepare_dataset_metadata(
    intercomparison_info: Mapping[str, Any],
    stage_store: StageStoreProtocol,
) -> dict[str, dict[str, Any]]:
    """Load dataset-level metadata once per dataset under its configured QA test."""

    out: dict[str, dict[str, Any]] = {}
    for dataset_id, dataset_cfg in intercomparison_info["datasets"].items():
        qa_test = dataset_cfg["qa_test"]
        metadata = {
            parameter: _load_optional_metadata(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=parameter,
            )
            for parameter in DATASET_METADATA_PARAMETERS
        }
        out[dataset_id] = {
            "system_label": dataset_cfg.get("system_label"),
            "dataset_label": dataset_cfg.get("dataset_label"),
            "reference": bool(dataset_cfg.get("reference", False)),
            "qa_test": qa_test,
            "metadata": metadata,
        }
    return out


def _prepare_channel_groups(
    intercomparison_info: Mapping[str, Any],
    stage_store: StageStoreProtocol,
) -> dict[str, dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    reference_dataset = intercomparison_info["reference_dataset"]
    dataset_configs = intercomparison_info["datasets"]

    for group_id, group in intercomparison_info.get("channel_groups", {}).items():
        datasets: dict[str, dict[str, Any]] = {}

        for dataset_id, group_dataset in group["datasets"].items():
            if not group_dataset.get("atlas_channel_id"):
                continue
            dataset_cfg = dataset_configs[dataset_id]
            qa_test = dataset_cfg["qa_test"]
            signal_source = dataset_cfg["signal_source"]
            error_source = dataset_cfg["signal_error_source"]
            channel_id = group_dataset["atlas_channel_id"]

            signal_entry = _load_entry(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=signal_source,
                config_field="signal_source",
                group_kind="channel_group",
                group_id=group_id,
            )
            error_entry = _load_entry(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=error_source,
                config_field="signal_error_source",
                group_kind="channel_group",
                group_id=group_id,
            )

            metadata = {}
            for parameter in CHANNEL_METADATA_PARAMETERS:
                metadata_entry = _load_optional_metadata(
                    stage_store,
                    dataset_id=dataset_id,
                    qa_test=qa_test,
                    parameter=parameter,
                )
                metadata[parameter] = _select_optional_product_metadata(
                    metadata_entry,
                    channel_id,
                    kind="channel",
                    location=(
                        f"channel group {group_id!r}, dataset {dataset_id!r}, "
                        f"metadata {parameter!r}"
                    ),
                )

            molecular_entry = _load_optional_metadata(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=CHANNEL_MOLECULAR_PROFILE_SOURCE,
            )
            molecular_profile = _select_optional_product_metadata(
                molecular_entry,
                channel_id,
                kind="channel",
                location=(
                    f"channel group {group_id!r}, dataset {dataset_id!r}, "
                    f"molecular source {CHANNEL_MOLECULAR_PROFILE_SOURCE!r}"
                ),
            )

            datasets[dataset_id] = {
                "system_label": dataset_cfg.get("system_label"),
                "dataset_label": dataset_cfg.get("dataset_label"),
                "qa_test": qa_test,
                "signal_source": signal_source,
                "signal_error_source": error_source,
                "atlas_channel_id": channel_id,
                "metadata": metadata,
                "molecular": {
                    "profile_source": CHANNEL_MOLECULAR_PROFILE_SOURCE,
                    "profile": molecular_profile,
                    "metadata_source": None,
                    "metadata": None,
                },
                "signal": _select_product(
                    signal_entry,
                    channel_id,
                    kind="channel",
                    location=(
                        f"channel group {group_id!r}, dataset {dataset_id!r}, "
                        f"source {signal_source!r}"
                    ),
                ),
                "error": _select_product(
                    error_entry,
                    channel_id,
                    kind="channel",
                    location=(
                        f"channel group {group_id!r}, dataset {dataset_id!r}, "
                        f"source {error_source!r}"
                    ),
                ),
            }

        if not datasets:
            continue

        groups[group_id] = {
            "label": group["label"],
            "reference_dataset": reference_dataset,
            "reference_atlas_channel_id": group["reference_atlas_channel_id"],
            "background_correction": group["background_correction"],
            "background_region": group["background_region"],
            "normalisation": group["normalisation"],
            "normalisation_region": group["normalisation_region"],
            "normalise_to_molecular": group["normalise_to_molecular"],
            "plot_molecular": group["plot_molecular"],
            "datasets": datasets,
        }

    return groups


def _prepare_pair_groups(
    intercomparison_info: Mapping[str, Any],
    stage_store: StageStoreProtocol,
) -> dict[str, dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    reference_dataset = intercomparison_info["reference_dataset"]
    dataset_configs = intercomparison_info["datasets"]

    for group_id, group in intercomparison_info.get("pair_groups", {}).items():
        datasets: dict[str, dict[str, Any]] = {}

        for dataset_id, group_dataset in group["datasets"].items():
            if not group_dataset.get("atlas_pair_id"):
                continue
            dataset_cfg = dataset_configs[dataset_id]
            qa_test = dataset_cfg["qa_test"]
            pair_source = dataset_cfg["pair_source"]
            error_source = dataset_cfg["pair_error_source"]
            pair_id = group_dataset["atlas_pair_id"]
            molecular_pair_id = _derive_molecular_pair_id(
                pair_id,
                location=f"pair group {group_id!r}, dataset {dataset_id!r}",
            )

            pair_entry = _load_entry(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=pair_source,
                config_field="pair_source",
                group_kind="pair_group",
                group_id=group_id,
            )
            error_entry = _load_entry(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=error_source,
                config_field="pair_error_source",
                group_kind="pair_group",
                group_id=group_id,
            )

            metadata = {}
            for parameter in PAIR_METADATA_PARAMETERS:
                metadata_entry = _load_optional_metadata(
                    stage_store,
                    dataset_id=dataset_id,
                    qa_test=qa_test,
                    parameter=parameter,
                )
                metadata[parameter] = _select_optional_product_metadata(
                    metadata_entry,
                    pair_id,
                    kind="pair",
                    location=(
                        f"pair group {group_id!r}, dataset {dataset_id!r}, "
                        f"metadata {parameter!r}"
                    ),
                )

            molecular_profile_entry = _load_optional_metadata(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=PAIR_MOLECULAR_PROFILE_SOURCE,
            )
            molecular_profile = _select_optional_molecular_pair_metadata(
                molecular_profile_entry,
                molecular_pair_id,
                atlas_pair_id=pair_id,
                parameter=PAIR_MOLECULAR_PROFILE_SOURCE,
                group_id=group_id,
                dataset_id=dataset_id,
            )

            molecular_info_entry = _load_optional_metadata(
                stage_store,
                dataset_id=dataset_id,
                qa_test=qa_test,
                parameter=PAIR_MOLECULAR_METADATA_SOURCE,
            )
            molecular_info = _select_optional_molecular_pair_metadata(
                molecular_info_entry,
                molecular_pair_id,
                atlas_pair_id=pair_id,
                parameter=PAIR_MOLECULAR_METADATA_SOURCE,
                group_id=group_id,
                dataset_id=dataset_id,
            )

            pair_channel_id = _pair_reference_channel_id(
                metadata.get("pol_cal_info"),
                group_id=group_id,
                dataset_id=dataset_id,
            )
            metadata.update(
                _load_pair_vertical_metadata(
                    stage_store,
                    dataset_id=dataset_id,
                    qa_test=qa_test,
                    channel_id=pair_channel_id,
                    group_id=group_id,
                )
            )

            datasets[dataset_id] = {
                "system_label": dataset_cfg.get("system_label"),
                "dataset_label": dataset_cfg.get("dataset_label"),
                "qa_test": qa_test,
                "pair_source": pair_source,
                "pair_error_source": error_source,
                "atlas_pair_id": pair_id,
                "molecular_pair_id": molecular_pair_id,
                "vertical_channel_id": pair_channel_id,
                "metadata": metadata,
                "molecular": {
                    "profile_source": PAIR_MOLECULAR_PROFILE_SOURCE,
                    "profile": molecular_profile,
                    "metadata_source": PAIR_MOLECULAR_METADATA_SOURCE,
                    "metadata": molecular_info,
                },
                "signal": _select_product(
                    pair_entry,
                    pair_id,
                    kind="pair",
                    location=(
                        f"pair group {group_id!r}, dataset {dataset_id!r}, "
                        f"source {pair_source!r}"
                    ),
                ),
                "error": _select_product(
                    error_entry,
                    pair_id,
                    kind="pair",
                    location=(
                        f"pair group {group_id!r}, dataset {dataset_id!r}, "
                        f"source {error_source!r}"
                    ),
                ),
            }

        if not datasets:
            continue

        groups[group_id] = {
            "label": group["label"],
            "reference_dataset": reference_dataset,
            "reference_atlas_pair_id": group["reference_atlas_pair_id"],
            "background_correction": group["background_correction"],
            "background_region": group["background_region"],
            "normalisation": group["normalisation"],
            "normalisation_region": group["normalisation_region"],
            "plot_molecular": group["plot_molecular"],
            "datasets": datasets,
        }

    return groups


def prepare_intercomparison_bundles(
    intercomparison_info: Mapping[str, Any],
    stage_store: StageStoreProtocol,
) -> dict[str, dict[str, dict[str, Any]]]:
    """Build comparison-oriented groups without modifying source data."""

    return {
        "datasets": _prepare_dataset_metadata(intercomparison_info, stage_store),
        "channel_groups": _prepare_channel_groups(intercomparison_info, stage_store),
        "pair_groups": _prepare_pair_groups(intercomparison_info, stage_store),
    }
