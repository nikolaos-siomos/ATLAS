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



def _select_channel_molecular_profile(value: Any, *, location: str) -> Any:
    """Select the attenuated molecular backscatter profile used for channels.

    The exported ``molecular`` source contains several optical parameters along
    ``opto_parameters``.  Intercomparison processing uses only ``atten_bsc``.
    """

    if value is None:
        return None

    try:
        return value.sel(opto_parameters="atten_bsc")
    except Exception:
        available = _available_ids(value, "opto_parameters")
        extra = (
            f" Available opto_parameters: {available}."
            if available
            else ""
        )
        raise ValueError(
            f"{location}: molecular source does not contain "
            f"opto_parameters='atten_bsc'.{extra}"
        ) from None

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



def _diagnostic_scalar(value: Any) -> Any:
    """Resolve a lazy scalar for diagnostic printing only."""
    try:
        if hasattr(value, "compute"):
            value = value.compute()
        if hasattr(value, "item"):
            value = value.item()
    except Exception:
        pass
    return value


def _diagnostic_finite_count(value: Any) -> int | None:
    """Count finite/non-null values without modifying the source object."""
    try:
        import numpy as np

        finite = np.isfinite(value)
        count = finite.sum()
        return int(_diagnostic_scalar(count))
    except Exception:
        try:
            count = value.notnull().sum()
            return int(_diagnostic_scalar(count))
        except Exception:
            return None


def _diagnostic_size(value: Any) -> int | None:
    try:
        return int(value.size)
    except Exception:
        return None


def _diagnose_channel_selection(
    *,
    group_id: str,
    entry_id: str,
    dataset_id: str,
    qa_test: str,
    signal_source: str,
    error_source: str,
    channel_id: str,
    signal_entry: Any,
    selected_signal: Any,
    selected_error: Any,
) -> None:
    """Print read-only diagnostics for one channel entry."""

    print()
    print("    BUNDLE DIAGNOSTIC")
    print(f"        group: [channel_group:{group_id}]")
    print(f"        entry: {entry_id!r}")
    print(f"        dataset: {dataset_id!r}")
    print(f"        qa_test: {qa_test!r}")
    print(f"        signal_source: {signal_source!r}")
    print(f"        error_source: {error_source!r}")
    print(f"        requested channel: {channel_id!r}")

    try:
        print(
            f"        raw signal dims / shape / dtype: "
            f"{signal_entry.dims} / {signal_entry.shape} / {signal_entry.dtype}"
        )
    except Exception:
        print(f"        raw signal type: {type(signal_entry)!r}")

    try:
        channel_values = [str(v) for v in signal_entry.coords["channel"].values.tolist()]
        print(f"        raw signal channel count: {len(channel_values)}")
        print(f"        requested channel present: {channel_id.lower() in [v.lower() for v in channel_values]}")
        print(f"        available channels: {', '.join(channel_values[:20])}" + (" ..." if len(channel_values) > 20 else ""))
    except Exception as exc:
        print(f"        could not inspect channel coordinate: {exc}")

    signal_count = _diagnostic_finite_count(selected_signal)
    signal_size = _diagnostic_size(selected_signal)
    error_count = _diagnostic_finite_count(selected_error)
    error_size = _diagnostic_size(selected_error)

    try:
        print(
            f"        selected signal dims / shape / dtype: "
            f"{selected_signal.dims} / {selected_signal.shape} / {selected_signal.dtype}"
        )
    except Exception:
        pass
    print(f"        selected signal finite values: {signal_count}/{signal_size}")

    try:
        print(
            f"        selected error dims / shape / dtype: "
            f"{selected_error.dims} / {selected_error.shape} / {selected_error.dtype}"
        )
    except Exception:
        pass
    print(f"        selected error finite values: {error_count}/{error_size}")

    try:
        print(f"        selected signal attrs: {dict(selected_signal.attrs)}")
    except Exception:
        pass


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
    entry_id: str | None = None,
) -> Any:
    """Select a derived molecular pair ID from optional pair metadata."""

    if value is None:
        return None

    try:
        return value.sel(pair=molecular_pair_id)
    except Exception:
        available = _available_ids(value, "pair")
        extra = f" Available pair IDs: {available}." if available else ""
        entry_part = f", entry {entry_id!r}" if entry_id is not None else ""
        raise ValueError(
            f"pair group {group_id!r}{entry_part}, dataset {dataset_id!r}, metadata {parameter!r}: "
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
    dataset_configs = intercomparison_info["datasets"]

    for group_id, group in intercomparison_info.get("channel_groups", {}).items():
        entries: dict[str, dict[str, Any]] = {}

        for entry_id, group_entry in group["entries"].items():
            dataset_id = group_entry["dataset_id"]
            dataset_cfg = dataset_configs[dataset_id]
            qa_test = dataset_cfg["qa_test"]
            signal_source = dataset_cfg["signal_source"]
            error_source = dataset_cfg["signal_error_source"]
            channel_id = group_entry["atlas_channel_id"]

            signal_entry = _load_entry(
                stage_store, dataset_id=dataset_id, qa_test=qa_test,
                parameter=signal_source, config_field="signal_source",
                group_kind="channel_group", group_id=group_id,
            )
            error_entry = _load_entry(
                stage_store, dataset_id=dataset_id, qa_test=qa_test,
                parameter=error_source, config_field="signal_error_source",
                group_kind="channel_group", group_id=group_id,
            )

            metadata = {}
            for parameter in CHANNEL_METADATA_PARAMETERS:
                metadata_entry = _load_optional_metadata(
                    stage_store, dataset_id=dataset_id, qa_test=qa_test, parameter=parameter,
                )
                metadata[parameter] = _select_optional_product_metadata(
                    metadata_entry, channel_id, kind="channel",
                    location=(
                        f"channel group {group_id!r}, entry {entry_id!r}, "
                        f"dataset {dataset_id!r}, metadata {parameter!r}"
                    ),
                )

            molecular_entry = _load_optional_metadata(
                stage_store, dataset_id=dataset_id, qa_test=qa_test,
                parameter=CHANNEL_MOLECULAR_PROFILE_SOURCE,
            )
            molecular_profile = _select_optional_product_metadata(
                molecular_entry, channel_id, kind="channel",
                location=(
                    f"channel group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}, "
                    f"molecular source {CHANNEL_MOLECULAR_PROFILE_SOURCE!r}"
                ),
            )
            molecular_profile = _select_channel_molecular_profile(
                molecular_profile,
                location=(
                    f"channel group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}, "
                    f"molecular source {CHANNEL_MOLECULAR_PROFILE_SOURCE!r}"
                ),
            )

            selected_signal = _select_product(
                signal_entry, channel_id, kind="channel",
                location=(
                    f"channel group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}, "
                    f"source {signal_source!r}"
                ),
            )
            selected_error = _select_product(
                error_entry, channel_id, kind="channel",
                location=(
                    f"channel group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}, "
                    f"source {error_source!r}"
                ),
            )

            _diagnose_channel_selection(
                group_id=group_id,
                entry_id=entry_id,
                dataset_id=dataset_id,
                qa_test=qa_test,
                signal_source=signal_source,
                error_source=error_source,
                channel_id=channel_id,
                signal_entry=signal_entry,
                selected_signal=selected_signal,
                selected_error=selected_error,
            )

            entries[entry_id] = {
                "dataset_id": dataset_id,
                "entry_label": group_entry.get("label"),
                "reference": bool(group_entry.get("reference", False)),
                "system_label": dataset_cfg.get("system_label"),
                "dataset_label": dataset_cfg.get("dataset_label"),
                "qa_test": qa_test,
                "signal_source": signal_source,
                "signal_error_source": error_source,
                "atlas_channel_id": channel_id,
                "metadata": metadata,
                "molecular": {
                    "profile_source": CHANNEL_MOLECULAR_PROFILE_SOURCE,
                    "profile_parameter": "atten_bsc",
                    "profile": molecular_profile,
                    "metadata_source": None,
                    "metadata": None,
                },
                "signal": selected_signal,
                "error": selected_error,
            }

        if not entries:
            continue

        reference_entry = group["reference_entry"]
        groups[group_id] = {
            "label": group["label"],
            "reference_entry": reference_entry,
            "reference_dataset": entries[reference_entry]["dataset_id"],
            "reference_atlas_channel_id": group["reference_atlas_channel_id"],
            "background_correction": group["background_correction"],
            "background_region": group["background_region"],
            "normalisation": group["normalisation"],
            "normalisation_region": group["normalisation_region"],
            "normalise_to_molecular": group["normalise_to_molecular"],
            "plot_molecular": group["plot_molecular"],
            "entries": entries,
        }

    return groups


def _prepare_pair_groups(
    intercomparison_info: Mapping[str, Any],
    stage_store: StageStoreProtocol,
) -> dict[str, dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    dataset_configs = intercomparison_info["datasets"]

    for group_id, group in intercomparison_info.get("pair_groups", {}).items():
        entries: dict[str, dict[str, Any]] = {}

        for entry_id, group_entry in group["entries"].items():
            dataset_id = group_entry["dataset_id"]
            dataset_cfg = dataset_configs[dataset_id]
            qa_test = dataset_cfg["qa_test"]
            pair_source = dataset_cfg["pair_source"]
            error_source = dataset_cfg["pair_error_source"]
            pair_id = group_entry["atlas_pair_id"]
            molecular_pair_id = _derive_molecular_pair_id(
                pair_id,
                location=f"pair group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}",
            )

            pair_entry = _load_entry(
                stage_store, dataset_id=dataset_id, qa_test=qa_test, parameter=pair_source,
                config_field="pair_source", group_kind="pair_group", group_id=group_id,
            )
            error_entry = _load_entry(
                stage_store, dataset_id=dataset_id, qa_test=qa_test, parameter=error_source,
                config_field="pair_error_source", group_kind="pair_group", group_id=group_id,
            )

            metadata = {}
            for parameter in PAIR_METADATA_PARAMETERS:
                metadata_entry = _load_optional_metadata(
                    stage_store, dataset_id=dataset_id, qa_test=qa_test, parameter=parameter,
                )
                metadata[parameter] = _select_optional_product_metadata(
                    metadata_entry, pair_id, kind="pair",
                    location=(
                        f"pair group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}, "
                        f"metadata {parameter!r}"
                    ),
                )

            molecular_profile_entry = _load_optional_metadata(
                stage_store, dataset_id=dataset_id, qa_test=qa_test,
                parameter=PAIR_MOLECULAR_PROFILE_SOURCE,
            )
            molecular_profile = _select_optional_molecular_pair_metadata(
                molecular_profile_entry, molecular_pair_id, atlas_pair_id=pair_id,
                parameter=PAIR_MOLECULAR_PROFILE_SOURCE, group_id=group_id, dataset_id=dataset_id,
                entry_id=entry_id,
            )

            molecular_info_entry = _load_optional_metadata(
                stage_store, dataset_id=dataset_id, qa_test=qa_test,
                parameter=PAIR_MOLECULAR_METADATA_SOURCE,
            )
            molecular_info = _select_optional_molecular_pair_metadata(
                molecular_info_entry, molecular_pair_id, atlas_pair_id=pair_id,
                parameter=PAIR_MOLECULAR_METADATA_SOURCE, group_id=group_id, dataset_id=dataset_id,
                entry_id=entry_id,
            )

            pair_channel_id = _pair_reference_channel_id(
                metadata.get("pol_cal_info"), group_id=group_id, dataset_id=dataset_id,
            )
            metadata.update(
                _load_pair_vertical_metadata(
                    stage_store, dataset_id=dataset_id, qa_test=qa_test,
                    channel_id=pair_channel_id, group_id=group_id,
                )
            )

            entries[entry_id] = {
                "dataset_id": dataset_id,
                "entry_label": group_entry.get("label"),
                "reference": bool(group_entry.get("reference", False)),
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
                    pair_entry, pair_id, kind="pair",
                    location=(
                        f"pair group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}, "
                        f"source {pair_source!r}"
                    ),
                ),
                "error": _select_product(
                    error_entry, pair_id, kind="pair",
                    location=(
                        f"pair group {group_id!r}, entry {entry_id!r}, dataset {dataset_id!r}, "
                        f"source {error_source!r}"
                    ),
                ),
            }

        if not entries:
            continue

        reference_entry = group["reference_entry"]
        groups[group_id] = {
            "label": group["label"],
            "reference_entry": reference_entry,
            "reference_dataset": entries[reference_entry]["dataset_id"],
            "reference_atlas_pair_id": group["reference_atlas_pair_id"],
            "background_correction": group["background_correction"],
            "background_region": group["background_region"],
            "normalisation": group["normalisation"],
            "normalisation_region": group["normalisation_region"],
            "plot_molecular": group["plot_molecular"],
            "entries": entries,
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
