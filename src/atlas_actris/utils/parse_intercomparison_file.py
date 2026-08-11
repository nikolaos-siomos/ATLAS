#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse and validate the ATLAS intercomparison initialization file.

The parser is intentionally limited to the first module phase: it validates the
INI structure, resolves ordinary defaults, expands channel/pair IDs to all
configured systems, resolves per-system signal/error overrides, and returns a
nested dictionary. Metadata-dependent defaults (for example lidar labels and
reference-channel normalization/background regions) remain ``None`` and are
resolved later, after exported stages are loaded.
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

VERTICAL_SCALES = ("height_asl", "height_agl", "range")
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
    "vertical_binning": {"dtype": float, "default": None, "min": 0.0, "category": "optional"},
    "vertical_min": {"dtype": float, "default": None, "category": "optional"},
    "vertical_max": {"dtype": float, "default": None, "category": "optional"},
    "slice_measurement": {"dtype": str, "default": [], "is_list": True, "category": "optional"},
    "exclude_measurement": {"dtype": str, "default": [], "is_list": True, "category": "optional"},
    "background_correction": {"dtype": bool, "default": False, "category": "optional"},
    "background_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "normalisation": {"dtype": bool, "default": True, "category": "optional"},
    "normalisation_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "normalise_to_molecular": {"dtype": bool, "default": True, "category": "optional"},
    "plot_molecular": {"dtype": bool, "default": True, "category": "optional"},
    "dpi": {"dtype": int, "default": 150, "min": 1, "category": "optional"},
    "color_reduction": {"dtype": bool, "default": False, "category": "optional"},
}

SYSTEM_SCHEMA: dict[str, dict[str, Any]] = {
    "stage_path": {"dtype": Path, "default": None, "category": "mandatory"},
    "reference": {"dtype": bool, "default": False, "category": "optional"},
    "label": {"dtype": str, "default": None, "category": "optional"},
}

CHANNEL_SCHEMA: dict[str, dict[str, Any]] = {
    "label": {"dtype": str, "default": None, "category": "optional"},
    "force_qa_test": {"dtype": str, "default": None, "allowed": QA_TESTS, "category": "optional"},
    "background_correction": {"dtype": bool, "default": None, "category": "optional"},
    "background_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "normalisation": {"dtype": bool, "default": None, "category": "optional"},
    "normalisation_region": {"dtype": float, "default": None, "is_list": True, "size": 2, "category": "optional"},
    "normalise_to_molecular": {"dtype": bool, "default": None, "category": "optional"},
    "plot_molecular": {"dtype": bool, "default": None, "category": "optional"},
}

PAIR_SCHEMA: dict[str, dict[str, Any]] = {
    "label": {"dtype": str, "default": None, "category": "optional"},
    "force_qa_test": {"dtype": str, "default": None, "allowed": QA_TESTS, "category": "optional"},
}

SECTION_SCHEMAS = {
    "general": GENERAL_SCHEMA,
    "system": SYSTEM_SCHEMA,
    "channel": CHANNEL_SCHEMA,
    "pair": PAIR_SCHEMA,
}

_SECTION_PATTERN = re.compile(r"^(general|system|channel|pair)(?::(.+))?$")


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


def _parse_dynamic_data_section(
    parser: configparser.ConfigParser,
    section: str,
    schema: Mapping[str, Mapping[str, Any]],
    *,
    system_ids: set[str],
    product_prefix: str,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    raw = parser[section]
    fixed_raw: dict[str, str] = {}
    system_raw: dict[str, dict[str, str]] = {system_id: {} for system_id in system_ids}

    if product_prefix == "channel":
        allowed_override_names = {"signal_source", "signal_error_source"}
    elif product_prefix == "pair":
        allowed_override_names = {"pair_source", "pair_error_source"}
    else:
        raise ValueError(f"Unsupported product prefix: {product_prefix!r}")

    for key, value in raw.items():
        if key in schema:
            fixed_raw[key] = value
            continue

        if "." in key:
            system_id, override = key.split(".", 1)
            if system_id not in system_ids:
                raise ValueError(f"[{section}] {key}: unknown system ID {system_id!r}")
            if override not in allowed_override_names:
                raise ValueError(f"[{section}] {key}: unsupported system override {override!r}")
            system_raw[system_id][override] = value
            continue

        if key in system_ids:
            system_raw[key]["atlas_id"] = value
            continue

        raise ValueError(f"[{section}]: unknown parameter {key!r}")

    values: dict[str, Any] = {}
    for key, meta in schema.items():
        value = _parse_value(fixed_raw.get(key, ""), meta, location=f"[{section}] {key}")
        if meta.get("category") == "mandatory" and value is None:
            raise ValueError(f"[{section}] {key}: mandatory value is missing")
        values[key] = value

    parsed_system_raw: dict[str, dict[str, Any]] = {}
    for system_id, overrides in system_raw.items():
        parsed: dict[str, Any] = {}
        if "atlas_id" in overrides:
            token = overrides["atlas_id"].strip()
            parsed["atlas_id"] = None if token.lower() == "off" else (token or None)
            parsed["disabled"] = token.lower() == "off"
        for override_name in allowed_override_names:
            if override_name in overrides:
                parsed[override_name] = overrides[override_name].strip() or None
        parsed_system_raw[system_id] = parsed

    return values, parsed_system_raw


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
        "system": [],
        "channel": [],
        "pair": [],
    }
    for section in parser.sections():
        match = _SECTION_PATTERN.match(section)
        if not match:
            raise ValueError(
                f"Unknown section [{section}]. Allowed sections are [general], "
                "[system:<id>], [channel:<id>], and [pair:<id>]"
            )
        kind, section_id = match.groups()
        if kind == "general":
            if section != "general":
                raise ValueError("The general section must be named exactly [general]")
            continue
        if not section_id or not section_id.strip():
            raise ValueError(f"[{section}]: section ID is missing")
        section_index[kind].append((section, section_id.strip()))

    if not section_index["system"]:
        raise ValueError("At least one [system:<id>] section is required")
    if not section_index["channel"] and not section_index["pair"]:
        raise ValueError("At least one [channel:<id>] or [pair:<id>] section is required")

    ini_dir = ini_path.parent
    general = _parse_fixed_section(parser, "general", GENERAL_SCHEMA)
    general["output_folder"] = _resolve_path(general["output_folder"], ini_dir)
    general["slice_measurement"] = _validate_time_ranges(
        general["slice_measurement"], location="[general] slice_measurement"
    )
    general["exclude_measurement"] = _validate_time_ranges(
        general["exclude_measurement"], location="[general] exclude_measurement"
    )

    if general["vertical_method"] == "vertical_binning" and general["vertical_binning"] is None:
        raise ValueError("[general] vertical_binning is mandatory when vertical_method=vertical_binning")
    if general["vertical_method"] == "interpolation" and general["vertical_binning"] is not None:
        raise ValueError("[general] vertical_binning must be empty when vertical_method=interpolation")
    if (
        general["vertical_min"] is not None
        and general["vertical_max"] is not None
        and general["vertical_min"] >= general["vertical_max"]
    ):
        raise ValueError("[general] vertical_min must be smaller than vertical_max")

    systems: dict[str, dict[str, Any]] = {}
    for section, system_id in section_index["system"]:
        if system_id in systems:
            raise ValueError(f"Duplicate system ID: {system_id!r}")
        system = _parse_fixed_section(parser, section, SYSTEM_SCHEMA)
        system["stage_path"] = _resolve_path(system["stage_path"], ini_dir)
        systems[system_id] = system

    references = [system_id for system_id, system in systems.items() if system["reference"]]
    if len(references) != 1:
        raise ValueError(f"Exactly one system must set reference=True; found {references}")
    reference_system = references[0]

    system_ids = set(systems)

    channels: dict[str, dict[str, Any]] = {}
    for section, comparison_id in section_index["channel"]:
        if comparison_id in channels:
            raise ValueError(f"Duplicate channel comparison ID: {comparison_id!r}")
        values, system_values = _parse_dynamic_data_section(
            parser,
            section,
            CHANNEL_SCHEMA,
            system_ids=system_ids,
            product_prefix="channel",
        )

        reference_values = system_values[reference_system]
        if reference_values.get("disabled"):
            raise ValueError(
                f"[{section}] {reference_system}: the reference system cannot be disabled"
            )

        reference_id = reference_values.get("atlas_id")
        if reference_id is None:
            raise ValueError(
                f"[{section}] {reference_system}: the reference system channel ID is mandatory"
            )

        reference_signal = (
            reference_values.get("signal_source")
            or general["default_signal_source"]
        )
        reference_error = (
            reference_values.get("signal_error_source")
            or general["default_signal_error_source"]
        )

        resolved_systems: dict[str, dict[str, Any]] = {}
        for system_id in systems:
            raw_system = system_values[system_id]
            if raw_system.get("disabled"):
                continue

            atlas_id = raw_system.get("atlas_id") or reference_id
            signal = raw_system.get("signal_source") or reference_signal
            error = raw_system.get("signal_error_source") or reference_error
            resolved_systems[system_id] = {
                "atlas_channel_id": atlas_id,
                "signal_source": signal,
                "signal_error_source": error,
            }

        channels[comparison_id] = {
            "label": values["label"] or reference_id,
            "reference_channel_id": reference_id,
            "qa_test": values["force_qa_test"] or general["default_qa_test"],
            "systems": resolved_systems,
            "background_correction": (
                general["background_correction"]
                if values["background_correction"] is None
                else values["background_correction"]
            ),
            "background_region": (
                values["background_region"]
                if values["background_region"] not in (None, [])
                else general["background_region"]
            ),
            "normalisation": (
                general["normalisation"]
                if values["normalisation"] is None
                else values["normalisation"]
            ),
            "normalisation_region": (
                values["normalisation_region"]
                if values["normalisation_region"] not in (None, [])
                else general["normalisation_region"]
            ),
            "normalise_to_molecular": (
                general["normalise_to_molecular"]
                if values["normalise_to_molecular"] is None
                else values["normalise_to_molecular"]
            ),
            "plot_molecular": (
                general["plot_molecular"]
                if values["plot_molecular"] is None
                else values["plot_molecular"]
            ),
        }

    pairs: dict[str, dict[str, Any]] = {}
    for section, comparison_id in section_index["pair"]:
        if comparison_id in pairs:
            raise ValueError(f"Duplicate pair comparison ID: {comparison_id!r}")
        values, system_values = _parse_dynamic_data_section(
            parser,
            section,
            PAIR_SCHEMA,
            system_ids=system_ids,
            product_prefix="pair",
        )

        reference_values = system_values[reference_system]
        if reference_values.get("disabled"):
            raise ValueError(
                f"[{section}] {reference_system}: the reference system cannot be disabled"
            )

        reference_id = reference_values.get("atlas_id")
        if reference_id is None:
            raise ValueError(
                f"[{section}] {reference_system}: the reference system pair ID is mandatory"
            )

        reference_signal = (
            reference_values.get("pair_source")
            or general["default_pair_source"]
        )
        reference_error = (
            reference_values.get("pair_error_source")
            or general["default_pair_error_source"]
        )

        resolved_systems: dict[str, dict[str, Any]] = {}
        for system_id in systems:
            raw_system = system_values[system_id]
            if raw_system.get("disabled"):
                continue

            atlas_id = raw_system.get("atlas_id") or reference_id
            signal = raw_system.get("pair_source") or reference_signal
            error = raw_system.get("pair_error_source") or reference_error
            resolved_systems[system_id] = {
                "atlas_pair_id": atlas_id,
                "pair_source": signal,
                "pair_error_source": error,
            }

        pairs[comparison_id] = {
            "label": values["label"] or reference_id,
            "reference_pair_id": reference_id,
            "qa_test": values["force_qa_test"] or general["default_qa_test"],
            "systems": resolved_systems,
        }

    return {
        "ini_file": ini_path,
        "general": general,
        "reference_system": reference_system,
        "systems": systems,
        "channels": channels,
        "pairs": pairs,
    }


# Backward-friendly alias for the future command-line module.
parse_intercomparison_file = parse_intercomparison_ini
