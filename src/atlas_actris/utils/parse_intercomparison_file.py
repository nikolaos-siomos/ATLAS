#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse and validate the ATLAS intercomparison initialization file.

The parser validates the INI structure, resolves general and dataset defaults,
expands channel/pair group IDs across configured datasets, and returns a nested
dictionary. Metadata-dependent defaults remain unresolved until exported-stage
metadata are loaded later in the workflow.
"""

from __future__ import annotations

import configparser
import re
from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping


TIME_FORMATS = (
    "%H%M",
    "%Y%m%d",
    "%Y%m%d_%H",
    "%Y%m%d_%H%M",
    "%Y%m%d_%H%M%S",
)

QA_TESTS = (
    "ray",
    "drk",
    "pcb",
    "tlc",
    "tlc_rin",
    "ray_pcb",
    "pcb_aux",
    "trg",
    "dtm",
    "cam",
)

VERTICAL_SCALES = ("bins", "range", "height_agl", "height_asl")
VERTICAL_METHODS = ("interpolation", "vertical_binning")

GENERAL_SCHEMA: dict[str, dict[str, Any]] = {
    "output_folder": {"dtype": Path, "default": Path("./analysis"), "category": "optional"},
    "overwrite_output": {"dtype": bool, "default": False, "category": "optional"},
    "default_qa_test": {"dtype": str, "default": "ray", "allowed": QA_TESTS, "category": "optional"},
    "default_signal_source": {"dtype": str, "default": "profile", "category": "optional"},
    "default_signal_error_source": {"dtype": str, "default": "profile_error", "category": "optional"},
    "default_pair_source": {"dtype": str, "default": "pol_cal_ratio_mean", "category": "optional"},
    "default_pair_error_source": {"dtype": str, "default": "pol_cal_ratio_error_mean", "category": "optional"},
    "vertical_scale": {"dtype": str, "default": "height_asl", "allowed": VERTICAL_SCALES, "category": "optional"},
    "vertical_method": {"dtype": str, "default": "interpolation", "allowed": VERTICAL_METHODS, "category": "optional"},
    "vertical_bin_width": {"dtype": float, "default": 0.1, "min": 0.0, "category": "optional"},
    "first_bin_left_edge": {"dtype": float, "default": 0.0, "category": "optional"},
    "vertical_min": {"dtype": float, "default": None, "category": "optional"},
    "vertical_max": {"dtype": float, "default": None, "category": "optional"},
    "plot_native_scale": {"dtype": bool, "default": False, "category": "optional"},
    "slice_measurement": {"dtype": str, "default": [], "is_list": True, "category": "optional"},
    "exclude_measurement": {"dtype": str, "default": [], "is_list": True, "category": "optional"},
    "default_channel_background_correction": {"dtype": bool, "default": False, "category": "optional"},
    "default_channel_background_region": {"dtype": float, "default": [18.0, 22.0], "is_list": True, "size": 2, "category": "optional"},
    "default_channel_normalisation": {"dtype": bool, "default": True, "category": "optional"},
    "default_channel_normalisation_region": {"dtype": float, "default": [4.0, 6.0], "is_list": True, "size": 2, "category": "optional"},
    "default_channel_normalise_to_molecular": {"dtype": bool, "default": True, "category": "optional"},
    "default_channel_plot_molecular": {"dtype": bool, "default": True, "category": "optional"},
    "default_pair_background_correction": {"dtype": bool, "default": False, "category": "optional"},
    "default_pair_background_region": {"dtype": float, "default": [18.0, 22.0], "is_list": True, "size": 2, "category": "optional"},
    "default_pair_normalisation": {"dtype": bool, "default": False, "category": "optional"},
    "default_pair_normalisation_region": {"dtype": float, "default": [7.5, 9.0], "is_list": True, "size": 2, "category": "optional"},
    "default_pair_plot_molecular": {"dtype": bool, "default": True, "category": "optional"},
    "dpi": {"dtype": int, "default": 150, "min": 1, "category": "optional"},
    "color_reduction": {"dtype": bool, "default": False, "category": "optional"},
}

PLOTTING_SCHEMA: dict[str, dict[str, Any]] = {
    "x_lims": {"dtype": float, "default": [], "is_list": True, "size": 2, "category": "optional"},
    "x_tick": {"dtype": float, "default": 2.0, "min": 0.0, "category": "optional"},
    "relative_difference_lims": {"dtype": float, "default": [-0.4, 0.4], "is_list": True, "size": 2, "category": "optional"},
    "channel_y_lims": {"dtype": float, "default": [], "is_list": True, "size": 2, "category": "optional"},
    "channel_smooth": {"dtype": bool, "default": True, "category": "optional"},
    "channel_smoothing_range": {"dtype": float, "default": [0.05, 35.0], "is_list": True, "size": 2, "category": "optional"},
    "channel_smoothing_window": {"dtype": float, "default": 0.5, "min": 0.0, "category": "optional"},
    "pair_y_lims": {"dtype": float, "default": [], "is_list": True, "size": 2, "category": "optional"},
    "pair_smooth": {"dtype": bool, "default": False, "category": "optional"},
    "pair_smoothing_range": {"dtype": float, "default": [0.05, 15.0], "is_list": True, "size": 2, "category": "optional"},
    "pair_smoothing_window": {"dtype": float, "default": 0.5, "min": 0.0, "category": "optional"},
}

DATASET_SCHEMA: dict[str, dict[str, Any]] = {
    "stage_path": {"dtype": Path, "default": None, "category": "mandatory"},
    "reference": {"dtype": bool, "default": False, "category": "optional"},
    "system_label": {"dtype": str, "default": None, "category": "optional"},
    "dataset_label": {"dtype": str, "default": None, "category": "optional"},
    "qa_test": {"dtype": str, "default": None, "allowed": QA_TESTS, "category": "optional"},
    "signal_source": {"dtype": str, "default": None, "category": "optional"},
    "signal_error_source": {"dtype": str, "default": None, "category": "optional"},
    "pair_source": {"dtype": str, "default": None, "category": "optional"},
    "pair_error_source": {"dtype": str, "default": None, "category": "optional"},
}

CHANNEL_GROUP_SCHEMA: dict[str, dict[str, Any]] = {
    "label": {"dtype": str, "default": None, "category": "optional"},
    "background_correction": {"dtype": bool, "default": None, "category": "optional"},
    "background_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "normalisation": {"dtype": bool, "default": None, "category": "optional"},
    "normalisation_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "normalise_to_molecular": {"dtype": bool, "default": None, "category": "optional"},
    "plot_molecular": {"dtype": bool, "default": None, "category": "optional"},
    "vertical_bin_width": {"dtype": float, "default": None, "min": 0.0, "category": "optional"},
    "smooth": {"dtype": bool, "default": None, "category": "optional"},
    "smoothing_range": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "smoothing_window": {"dtype": float, "default": None, "min": 0.0, "category": "optional"},
    "x_lims": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "x_tick": {"dtype": float, "default": None, "min": 0.0, "category": "optional"},
    "y_lims": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "relative_difference_lims": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "use_log_y_scale": {"dtype": bool, "default": True, "category": "optional"},
}

PAIR_GROUP_SCHEMA: dict[str, dict[str, Any]] = {
    "label": {"dtype": str, "default": None, "category": "optional"},
    "background_correction": {"dtype": bool, "default": None, "category": "optional"},
    "background_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "normalisation": {"dtype": bool, "default": None, "category": "optional"},
    "normalisation_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "plot_molecular": {"dtype": bool, "default": None, "category": "optional"},
    "vertical_bin_width": {"dtype": float, "default": None, "min": 0.0, "category": "optional"},
    "smooth": {"dtype": bool, "default": None, "category": "optional"},
    "smoothing_range": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "smoothing_window": {"dtype": float, "default": None, "min": 0.0, "category": "optional"},
    "x_lims": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "x_tick": {"dtype": float, "default": None, "min": 0.0, "category": "optional"},
    "y_lims": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "relative_difference_lims": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "use_log_y_scale": {"dtype": bool, "default": False, "category": "optional"},
}

SECTION_SCHEMAS = {
    "general": GENERAL_SCHEMA,
    "plotting": PLOTTING_SCHEMA,
    "dataset": DATASET_SCHEMA,
    "channel_group": CHANNEL_GROUP_SCHEMA,
    "pair_group": PAIR_GROUP_SCHEMA,
}

_SECTION_PATTERN = re.compile(r"^(general|plotting|dataset|channel_group|pair_group)(?::(.+))?$")


def _split_list(value: str) -> list[str]:
    return [part.strip() for part in re.split(r"[,;]", value) if part.strip()]


def _parse_bool(value: str, *, location: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "yes", "1", "on"}:
        return True
    if normalized in {"false", "no", "0", "off"}:
        return False
    raise ValueError(f"{location}: expected a boolean, got {value!r}")


def _parse_scalar(value: str, dtype: type, *, location: str) -> Any:
    if dtype is str:
        return value.strip() or None
    if dtype is bool:
        return _parse_bool(value, location=location)
    if dtype is int:
        try:
            return int(value)
        except ValueError as exc:
            raise ValueError(f"{location}: expected an integer, got {value!r}") from exc
    if dtype is float:
        try:
            return float(value)
        except ValueError as exc:
            raise ValueError(f"{location}: expected a number, got {value!r}") from exc
    if dtype is Path:
        return Path(value.strip())
    raise TypeError(f"{location}: unsupported schema dtype {dtype!r}")


def _parse_value(raw: str, meta: Mapping[str, Any], *, location: str) -> Any:
    if not raw.strip():
        return deepcopy(meta.get("default"))

    dtype = meta["dtype"]
    if meta.get("is_list"):
        parts = _split_list(raw)
        value = [_parse_scalar(part, dtype, location=location) for part in parts]
    else:
        value = _parse_scalar(raw, dtype, location=location)

    allowed = meta.get("allowed")
    if allowed is not None and value not in allowed:
        raise ValueError(f"{location}: expected one of {tuple(allowed)}, got {value!r}")

    if value is not None and not isinstance(value, list):
        if "min" in meta and value <= meta["min"]:
            raise ValueError(f"{location}: value must be greater than {meta['min']}")
        if "max" in meta and value > meta["max"]:
            raise ValueError(f"{location}: value must be at most {meta['max']}")

    size = meta.get("size")
    if size is not None and value not in (None, []):
        allowed_sizes = size if isinstance(size, list) else [size]
        if len(value) not in allowed_sizes:
            raise ValueError(f"{location}: expected list size {allowed_sizes}, got {len(value)}")

    return value


def _resolve_path(path: Path | None, ini_dir: Path) -> Path | None:
    if path is None:
        return None
    if not path.is_absolute():
        path = ini_dir / path
    return path.expanduser().resolve()


def _validate_time_value(value: str, *, location: str) -> None:
    from datetime import datetime

    for fmt in TIME_FORMATS:
        try:
            datetime.strptime(value, fmt)
            return
        except ValueError:
            continue
    raise ValueError(
        f"{location}: unsupported time {value!r}. Accepted formats: "
        "HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, yyyymmdd_HHMMSS"
    )


def _validate_time_ranges(values: list[str], *, location: str) -> list[tuple[str, str]]:
    if len(values) % 2 != 0:
        raise ValueError(f"{location}: values must be provided as repeating start, stop pairs")
    ranges: list[tuple[str, str]] = []
    for index in range(0, len(values), 2):
        start, stop = values[index], values[index + 1]
        _validate_time_value(start, location=f"{location}[{index}]")
        _validate_time_value(stop, location=f"{location}[{index + 1}]")
        ranges.append((start, stop))
    return ranges


def _parse_fixed_section(
    parser: configparser.ConfigParser,
    section: str,
    schema: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    raw = parser[section]
    unknown = set(raw) - set(schema)
    if unknown:
        raise ValueError(f"[{section}]: unknown parameters: {sorted(unknown)}")

    result: dict[str, Any] = {}
    for key, meta in schema.items():
        value = _parse_value(raw.get(key, ""), meta, location=f"[{section}] {key}")
        if meta.get("category") == "mandatory" and value is None:
            raise ValueError(f"[{section}] {key}: mandatory value is missing")
        result[key] = value
    return result


def _parse_group_section(
    parser: configparser.ConfigParser,
    section: str,
    schema: Mapping[str, Mapping[str, Any]],
    *,
    dataset_ids: set[str],
    id_key: str,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    """Parse fixed group settings plus ``<dataset_id>.<id_key>`` mappings."""

    raw = parser[section]
    fixed_raw: dict[str, str] = {}
    dataset_raw: dict[str, dict[str, str]] = {dataset_id: {} for dataset_id in dataset_ids}

    for key, value in raw.items():
        if key in schema:
            fixed_raw[key] = value
            continue

        if "." not in key:
            raise ValueError(
                f"[{section}]: unknown parameter {key!r}. Dataset-specific mappings "
                f"must use '<dataset_id>.{id_key}'"
            )

        dataset_id, parameter = key.split(".", 1)
        if dataset_id not in dataset_ids:
            raise ValueError(f"[{section}] {key}: unknown dataset ID {dataset_id!r}")
        if parameter != id_key:
            raise ValueError(
                f"[{section}] {key}: unsupported dataset parameter {parameter!r}; "
                f"expected {id_key!r}"
            )
        dataset_raw[dataset_id][id_key] = value

    values: dict[str, Any] = {}
    for key, meta in schema.items():
        value = _parse_value(fixed_raw.get(key, ""), meta, location=f"[{section}] {key}")
        if meta.get("category") == "mandatory" and value is None:
            raise ValueError(f"[{section}] {key}: mandatory value is missing")
        values[key] = value

    parsed_datasets: dict[str, dict[str, Any]] = {}
    for dataset_id, mapping in dataset_raw.items():
        parsed: dict[str, Any] = {}
        if id_key in mapping:
            token = mapping[id_key].strip().lower()
            if token == "off":
                raise ValueError(
                    f"[{section}] {dataset_id}.{id_key}: 'off' is no longer used; "
                    "leave the value empty or omit the parameter to exclude this dataset"
                )
            parsed[id_key] = token or None
        parsed_datasets[dataset_id] = parsed

    return values, parsed_datasets


def _resolve_dataset_defaults(
    dataset: dict[str, Any],
    general: Mapping[str, Any],
) -> dict[str, Any]:
    dataset["qa_test"] = dataset["qa_test"] or general["default_qa_test"]
    dataset["signal_source"] = dataset["signal_source"] or general["default_signal_source"]
    dataset["signal_error_source"] = (
        dataset["signal_error_source"] or general["default_signal_error_source"]
    )
    dataset["pair_source"] = dataset["pair_source"] or general["default_pair_source"]
    dataset["pair_error_source"] = (
        dataset["pair_error_source"] or general["default_pair_error_source"]
    )
    return dataset




def _resolve_group_defaults(
    values: dict[str, Any],
    general: Mapping[str, Any],
    plotting: Mapping[str, Any],
    *,
    kind: str,
) -> dict[str, Any]:
    """Resolve empty group processing and plotting settings."""

    processing_prefix = f"default_{kind}_"
    processing_keys = [
        "background_correction",
        "background_region",
        "normalisation",
        "normalisation_region",
        "plot_molecular",
    ]
    if kind == "channel":
        processing_keys.append("normalise_to_molecular")

    for key in processing_keys:
        if values.get(key) is None:
            values[key] = deepcopy(general[processing_prefix + key])

    plotting_map = {
        "x_lims": "x_lims",
        "x_tick": "x_tick",
        "y_lims": f"{kind}_y_lims",
        "relative_difference_lims": "relative_difference_lims",
        "smooth": f"{kind}_smooth",
        "smoothing_range": f"{kind}_smoothing_range",
        "smoothing_window": f"{kind}_smoothing_window",
    }
    for key, plotting_key in plotting_map.items():
        if values.get(key) is None:
            values[key] = deepcopy(plotting[plotting_key])

    # use_log_y_scale intentionally comes only from the group schema:
    # True by default for channel groups, False by default for pair groups.
    # vertical_bin_width remains unresolved here. The vertical processor resolves
    # group override -> [general] -> coarsest native-step fallback.
    return values

def parse_intercomparison_ini(filepath: str | Path) -> dict[str, Any]:
    """Load and resolve an ATLAS intercomparison initialization file."""

    ini_path = Path(filepath).expanduser().resolve()
    if not ini_path.exists():
        raise FileNotFoundError(f"Intercomparison initialization file not found: {ini_path}")

    parser = configparser.ConfigParser(interpolation=None, strict=True)
    parser.optionxform = str.lower
    with ini_path.open("r", encoding="utf-8") as stream:
        parser.read_file(stream)

    if "general" not in parser:
        raise ValueError("The intercomparison INI must contain a [general] section")

    section_index: dict[str, list[tuple[str, str | None]]] = {
        "dataset": [],
        "channel_group": [],
        "pair_group": [],
    }
    for section in parser.sections():
        match = _SECTION_PATTERN.match(section)
        if not match:
            raise ValueError(
                f"Unknown section [{section}]. Allowed sections are [general], [plotting], "
                "[dataset:<id>], [channel_group:<id>], and [pair_group:<id>]"
            )
        kind, section_id = match.groups()
        if kind in {"general", "plotting"}:
            if section != kind:
                raise ValueError(f"The {kind} section must be named exactly [{kind}]")
            continue
        if not section_id or not section_id.strip():
            raise ValueError(f"[{section}]: section ID is missing")
        section_index[kind].append((section, section_id.strip()))

    if not section_index["dataset"]:
        raise ValueError("At least one [dataset:<id>] section is required")
    if not section_index["channel_group"] and not section_index["pair_group"]:
        raise ValueError(
            "At least one [channel_group:<id>] or [pair_group:<id>] section is required"
        )

    ini_dir = ini_path.parent
    general = _parse_fixed_section(parser, "general", GENERAL_SCHEMA)
    if "plotting" in parser:
        plotting = _parse_fixed_section(parser, "plotting", PLOTTING_SCHEMA)
    else:
        plotting = {key: deepcopy(meta.get("default")) for key, meta in PLOTTING_SCHEMA.items()}

    general["output_folder"] = _resolve_path(general["output_folder"], ini_dir)
    general["slice_measurement"] = _validate_time_ranges(
        general["slice_measurement"], location="[general] slice_measurement"
    )
    general["exclude_measurement"] = _validate_time_ranges(
        general["exclude_measurement"], location="[general] exclude_measurement"
    )

    if (
        general["vertical_min"] is not None
        and general["vertical_max"] is not None
        and general["vertical_min"] >= general["vertical_max"]
    ):
        raise ValueError("[general] vertical_min must be smaller than vertical_max")

    datasets: dict[str, dict[str, Any]] = {}
    for section, dataset_id in section_index["dataset"]:
        if dataset_id in datasets:
            raise ValueError(f"Duplicate dataset ID: {dataset_id!r}")
        dataset = _parse_fixed_section(parser, section, DATASET_SCHEMA)
        dataset["stage_path"] = _resolve_path(dataset["stage_path"], ini_dir)
        datasets[dataset_id] = _resolve_dataset_defaults(dataset, general)

    references = [
        dataset_id
        for dataset_id, dataset in datasets.items()
        if dataset["reference"]
    ]
    if len(references) != 1:
        raise ValueError(f"Exactly one dataset must set reference=True; found {references}")
    reference_dataset = references[0]
    dataset_ids = set(datasets)

    channel_groups: dict[str, dict[str, Any]] = {}
    for section, group_id in section_index["channel_group"]:
        if group_id in channel_groups:
            raise ValueError(f"Duplicate channel-group ID: {group_id!r}")

        values, dataset_values = _parse_group_section(
            parser,
            section,
            CHANNEL_GROUP_SCHEMA,
            dataset_ids=dataset_ids,
            id_key="atlas_channel_id",
        )
        values = _resolve_group_defaults(values, general, plotting, kind="channel")

        resolved_datasets: dict[str, dict[str, Any]] = {}
        for dataset_id in datasets:
            channel_id = dataset_values[dataset_id].get("atlas_channel_id")
            if channel_id is None:
                continue
            resolved_datasets[dataset_id] = {"atlas_channel_id": channel_id}

        # An empty group is treated as if it was not defined. This allows users
        # to keep placeholder sections in the INI without triggering imports.
        if not resolved_datasets:
            continue

        reference_id = resolved_datasets.get(reference_dataset, {}).get("atlas_channel_id")
        if reference_id is None:
            raise ValueError(
                f"[{section}] {reference_dataset}.atlas_channel_id: "
                "the reference dataset must provide a channel ID for an active group"
            )

        channel_groups[group_id] = {
            "label": values["label"] or reference_id,
            "reference_atlas_channel_id": reference_id,
            "datasets": resolved_datasets,
            "background_correction": values["background_correction"],
            "background_region": values["background_region"],
            "normalisation": values["normalisation"],
            "normalisation_region": values["normalisation_region"],
            "normalise_to_molecular": values["normalise_to_molecular"],
            "plot_molecular": values["plot_molecular"],
            "vertical_bin_width": values["vertical_bin_width"],
            "smooth": values["smooth"],
            "smoothing_range": values["smoothing_range"],
            "smoothing_window": values["smoothing_window"],
            "x_lims": values["x_lims"],
            "x_tick": values["x_tick"],
            "y_lims": values["y_lims"],
            "relative_difference_lims": values["relative_difference_lims"],
            "use_log_y_scale": values["use_log_y_scale"],
        }

    pair_groups: dict[str, dict[str, Any]] = {}
    for section, group_id in section_index["pair_group"]:
        if group_id in pair_groups:
            raise ValueError(f"Duplicate pair-group ID: {group_id!r}")

        values, dataset_values = _parse_group_section(
            parser,
            section,
            PAIR_GROUP_SCHEMA,
            dataset_ids=dataset_ids,
            id_key="atlas_pair_id",
        )
        values = _resolve_group_defaults(values, general, plotting, kind="pair")

        resolved_datasets: dict[str, dict[str, Any]] = {}
        for dataset_id in datasets:
            pair_id = dataset_values[dataset_id].get("atlas_pair_id")
            if pair_id is None:
                continue
            resolved_datasets[dataset_id] = {"atlas_pair_id": pair_id}

        # An empty group is treated as if it was not defined.
        if not resolved_datasets:
            continue

        reference_id = resolved_datasets.get(reference_dataset, {}).get("atlas_pair_id")
        if reference_id is None:
            raise ValueError(
                f"[{section}] {reference_dataset}.atlas_pair_id: "
                "the reference dataset must provide a pair ID for an active group"
            )

        pair_groups[group_id] = {
            "label": values["label"] or reference_id,
            "reference_atlas_pair_id": reference_id,
            "datasets": resolved_datasets,
            "background_correction": values["background_correction"],
            "background_region": values["background_region"],
            "normalisation": values["normalisation"],
            "normalisation_region": values["normalisation_region"],
            "plot_molecular": values["plot_molecular"],
            "vertical_bin_width": values["vertical_bin_width"],
            "smooth": values["smooth"],
            "smoothing_range": values["smoothing_range"],
            "smoothing_window": values["smoothing_window"],
            "x_lims": values["x_lims"],
            "x_tick": values["x_tick"],
            "y_lims": values["y_lims"],
            "relative_difference_lims": values["relative_difference_lims"],
            "use_log_y_scale": values["use_log_y_scale"],
        }

    return {
        "ini_file": ini_path,
        "general": general,
        "plotting": plotting,
        "reference_dataset": reference_dataset,
        "datasets": datasets,
        "channel_groups": channel_groups,
        "pair_groups": pair_groups,
    }


parse_intercomparison_file = parse_intercomparison_ini
