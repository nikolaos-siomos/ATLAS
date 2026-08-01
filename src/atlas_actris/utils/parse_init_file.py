#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:16:34 2025

@author: nikos
"""

from __future__ import annotations

import re
import os
import numpy as np
import configparser
from pathlib import Path
from pprint import pprint
from datetime import datetime
from utils.printouts import print_header, endpoint
from utils.error_classes import ConfigError, CustomWarning
from typing import Any, Dict, List, Optional, Union, Mapping, Iterable, Sequence
from utils.cookbook import collect_stages

Number = Union[int, float]

# -------------------------------------------------------------------
# SCHEMA
#   dtype:        expected Python type (str, int, float)
#   default:      scalar default (always scalar; expanded to lists as needed)
#   is_list:      True if value is a list in INI (comma or semicolon-separated)
#   allowed:      optional list of allowed values (checked only on non-empty values)
#   category:     "mandatory" | "recommended" | "optional"
#   min/max:      optional numeric bounds
#   check_path:   Indicates whether a string must be checked as file or folder for existence
#   size:         The list ust have a specific size if not empty
# -------------------------------------------------------------------

qa_tests = [
    "ray",
    "pcb",
    "tlc",
    "tlc_rin",
    "drk",
    "trg",
    "dtm",
    "ray_pcb",
    "pcb_aux",
    "cam",
]

quicklooks = [
    "ray",
    "pcb",
    "tlc",
    "tlc_rin",
    "drk",
    "ray_pcb",
    "pcb_aux",
    "vldr",
]

qa_measurement_folders = [
    "ray",
    "pcb_p45",
    "pcb_m45",
    "tlc_north",
    "tlc_east",
    "tlc_south",
    "tlc_west",
    "tlc_inner",
    "tlc_outer",
    "trg",
    "dtm",
    "ray_pcb",
    "pcb_aux_p45",
    "pcb_aux_m45",
    "drk",
    "drk_ray",
    "drk_pcb",
    "drk_tlc",
    "drk_tlc_rin",
    "drk_trg",
    "drk_dtm",
    "drk_ray_pcb",
    "drk_pcb_aux",
]

# User-facing bundle aliases accepted only by slice_measurement and
# exclude_measurement. These aliases are expanded to concrete measurement keys
# before the final caller_info dictionary is returned.
#
# Keep qa_tests, quicklooks, and qa_measurement_folders unchanged: external code
# may import or rely on those names and values.
slice_exclude_bundles = {
    "tlc": ["tlc_north", "tlc_east", "tlc_south", "tlc_west"],
    "tlc_rin": ["tlc_inner", "tlc_outer"],
    "pcb": ["pcb_p45", "pcb_m45"],
    "pcb_aux": ["pcb_aux_p45", "pcb_aux_m45"],
}

version_warning_nrm = (
    "\nSince update 1.0.0 the ray suffix has replaced the nrm suffix "
    "for Rayleigh measurements. Any folders with the nrm suffix "
    "will be renamed using the ray suffix\n"
)

version_warning_pcb = (
    "\nSince update 1.0.0 the p45 and m45 folders have replaced the "
    "+45 and -45 folders for pol. cal measurements. Any folder named "
    "+45 (-45) will be automatically renamed to p45 (m45), respectively\n"
)

slice_exclude_allowed_keys = sorted(
    set(qa_measurement_folders) | set(slice_exclude_bundles.keys())
)

height_units = ["m_asl", "m_agl", "km_asl", "km_agl"]
pressure_units = ["Pa", "hPa", "atm"]
temperature_units = ["K", "C", "Cx10"]
humidity_units = ["percent", "fraction"]

default_export_stage = ['pol_cal_complete']
default_mean_signal_stages = [
    'common_preprocessing_complete', 
    'preprocessing_complete', 
    'dark_preprocessing_complete'
    ]
default_signal_stages = ['screening_complete']

allowed_stages = collect_stages()
        
allowed_vertical_scales = ["range", "height_agl", "height_asl"]

SCHEMA: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "scc_compatible_format":     {"dtype": bool, "default": False, "is_list": False, "category": "optional"},
    "export_hoi_cfg":            {"dtype": str,  "default": "0",   "is_list": False, "category": "optional", "allowed": ["0", "1", "2"]},
    "scc_configuration_id":      {"dtype": str,  "default": None,  "is_list": False, "category": "optional"},

    "parent_folder":             {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "dir"},
    "atlas_configuration_file":  {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "file"},
    "atlas_settings_file":       {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "file"},
    "radiosonde_folder":         {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "dir"},
    "radiosonde_file":           {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "file"},

    "process":                 {"dtype": str,  "default": qa_tests,                   "is_list": True,  "category": "optional", "allowed": qa_tests + ["off"]},
    "process_qck":             {"dtype": str,  "default": quicklooks,                 "is_list": True,  "category": "optional", "allowed": quicklooks + ["off"]},
    "vertical_scale":          {"dtype": str,  "default": "range",                    "is_list": False, "category": "optional", "allowed": allowed_vertical_scales},
    "view_mean_signal_stages": {"dtype": str,  "default": default_mean_signal_stages, "is_list": True,  "category": "optional", "allowed": allowed_stages},
    "view_signal_stages":      {"dtype": str,  "default": default_signal_stages,      "is_list": True,  "category": "optional", "allowed": allowed_stages},
    "dpi":                     {"dtype": int,  "default": 150,                        "is_list": False, "category": "optional"},
    "color_reduction":         {"dtype": bool, "default": False,                      "is_list": False, "category": "optional"},
    "output_folder":           {"dtype": str,  "default": None,                       "is_list": False, "category": "optional"},
    "overwrite_output":        {"dtype": bool, "default": False,                      "is_list": False, "category": "optional"},
    "expert_analyst":          {"dtype": str,  "default": None,                       "is_list": False, "category": "optional"},
    "export_stages":           {"dtype": str,  "default": default_export_stage,       "is_list": True,  "category": "optional", "allowed": allowed_stages},
    "export_all":              {"dtype": bool, "default": False,                      "is_list": False, "category": "optional"},

    "select_channels":          {"dtype": str, "default": [], "is_list": True, "category": "optional"},
    "exclude_wavelength":       {"dtype": str, "default": [], "is_list": True, "category": "optional"},
    "exclude_telescope_type":   {"dtype": str, "default": [], "is_list": True, "category": "optional", "allowed": ["n", "f", "x", "m", "g", "y", "l", "h", "z"]},
    "exclude_channel_type":     {"dtype": str, "default": [], "is_list": True, "category": "optional", "allowed": ["p", "c", "t", "v", "r", "a", "f"]},
    "exclude_acquisition_mode": {"dtype": str, "default": [], "is_list": True, "category": "optional", "allowed": ["a", "p", "g"]},
    "exclude_channel_subtype":  {"dtype": str, "default": [], "is_list": True, "category": "optional", "allowed": ["r", "t", "n", "o", "w", "c", "h", "l", "a", "m", "b", "s", "x"]},

    "max_height_agl":               {"dtype": float, "default": 40.,     "is_list": False, "category": "optional"},
    "low_shot_threshold":           {"dtype": float, "default": 0.9,     "is_list": False, "category": "optional", "min": 0., "max": 0.999},
    "trim_overflows":               {"dtype": int,   "default": 0,       "is_list": False, "category": "optional", "allowed": [0, 1, 2, 3]},
    "low_res_averaging_rate":       {"dtype": str,   "default": '1H',    "is_list": False, "category": "optional"},
    "low_res_averaging_threshold":  {"dtype": float, "default": 0.5,     "is_list": False, "category": "optional", "min": 0., "max": 1.},
    "high_res_averaging_rate":      {"dtype": str,   "default": '10min', "is_list": False, "category": "optional"},
    "high_res_averaging_threshold": {"dtype": float, "default": 0.5,     "is_list": False, "category": "optional", "min": 0., "max": 1.},
    "max_adjacent_overflows":       {"dtype": int,   "default": 5,       "is_list": False, "category": "optional", "min": 0},
    "slice_measurement":            {"dtype": str,   "default": [],      "is_list": True,  "category": "optional"},
    "exclude_measurement":          {"dtype": str,   "default": [],      "is_list": True,  "category": "optional"},

    "ray":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "pcb":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "tlc":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "tlc_rin":      {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "drk":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "trg":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "dtm":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "ray_pcb":      {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "pcb_aux":      {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "cam":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "relative"},

    "files_per_quadrant":       {"dtype": int,   "default": None,     "is_list": False, "category": "optional", "min": 0},
    "files_per_ring":           {"dtype": int,   "default": None,     "is_list": False, "category": "optional", "min": 0},
    "rsonde_skip_header":       {"dtype": int,   "default": 1,        "is_list": False, "category": "optional", "min": 0},
    "rsonde_skip_footer":       {"dtype": int,   "default": 0,        "is_list": False, "category": "optional", "min": 0},
    "rsonde_delimiter":         {"dtype": str,   "default": "S",      "is_list": False, "category": "optional", "allowed": ["S", "C"]},
    "rsonde_column_index":      {"dtype": int,   "default": [2, 1, 3, 5], "is_list": True, "category": "optional", "min": 0, "size": [3, 4]},
    "rsonde_column_units":      {"dtype": str,   "default": ["m_asl", "hPa", "C", "percent"], "is_list": True, "category": "optional", "size": [3, 4], "allowed": height_units + pressure_units + temperature_units + humidity_units},
    "rsonde_station_altitude":  {"dtype": float, "default": None,     "is_list": False, "category": "optional", "min": 0,    "max": 5000},
    "cloudnet_station_name":    {"dtype": str,   "default": None,     "is_list": False, "category": "optional"},
    "rsonde_station_name":      {"dtype": str,   "default": None,     "is_list": False, "category": "optional"},
    "rsonde_station_wmo_id":    {"dtype": int,   "default": None,     "is_list": False, "category": "optional", "min": 0},
}


# -------------------------------------------------------------------
# INI SECTION LAYOUT
#
# The initialization parser stores all values in a flat caller_info dict.
# Sections are used only to validate user-facing INI files and to keep the
# templates/documentation organized.  Section names are intentionally strict
# lowercase strings.  Valid keys are still defined exclusively by SCHEMA.
# -------------------------------------------------------------------

INIT_FILE_SECTIONS: Dict[str, List[str]] = {
    "configuration": [
        "scc_compatible_format",
        "export_hoi_cfg",
        "scc_configuration_id",
    ],
    "explicit_paths": [
        "parent_folder",
        "atlas_configuration_file",
        "atlas_settings_file",
        "radiosonde_folder",
        "radiosonde_file",
        "output_folder",
    ],
    "general_options": [
        "process",
        "process_qck",
        "vertical_scale",
        "view_mean_signal_stages",
        "view_signal_stages",
        "dpi",
        "color_reduction",
        "overwrite_output",
        "expert_analyst",
        "export_stages",
        "export_all",
    ],
    "filter_channels": [
        "select_channels",
        "exclude_wavelength",
        "exclude_telescope_type",
        "exclude_channel_type",
        "exclude_acquisition_mode",
        "exclude_channel_subtype",
    ],
    "trimming_options": [
        "max_height_agl",
        "low_shot_threshold",
        "trim_overflows",
        "low_res_averaging_rate",
        "low_res_averaging_threshold",
        "high_res_averaging_rate",
        "high_res_averaging_threshold",
        "max_adjacent_overflows",
        "slice_measurement",
        "exclude_measurement",
    ],
    "explicit_folders": [
        "ray",
        "pcb",
        "tlc",
        "tlc_rin",
        "drk",
        "trg",
        "dtm",
        "ray_pcb",
        "pcb_aux",
        "cam",
    ],
    "parsing_options": [
        "files_per_quadrant",
        "files_per_ring",
        "rsonde_skip_header",
        "rsonde_skip_footer",
        "rsonde_delimiter",
        "rsonde_column_index",
        "rsonde_column_units",
        "rsonde_station_altitude",
        "cloudnet_station_name",
        "rsonde_station_name",
        "rsonde_station_wmo_id",
    ],
}


def _init_key_to_section() -> Dict[str, str]:
    key_to_section: Dict[str, str] = {}
    duplicates: List[str] = []

    for section, keys in INIT_FILE_SECTIONS.items():
        for key in keys:
            if key in key_to_section:
                duplicates.append(key)
            key_to_section[key] = section

    if duplicates:
        raise ConfigError(
            "Internal initialization section layout contains duplicated key(s): "
            f"{sorted(set(duplicates))}"
        )

    schema_keys = set(SCHEMA.keys())
    section_keys = set(key_to_section.keys())

    missing = schema_keys - section_keys
    unknown = section_keys - schema_keys

    if missing:
        raise ConfigError(
            "Internal initialization section layout is missing schema key(s): "
            f"{sorted(missing)}"
        )

    if unknown:
        raise ConfigError(
            "Internal initialization section layout contains key(s) not present "
            f"in SCHEMA: {sorted(unknown)}"
        )

    return key_to_section


INIT_KEY_TO_SECTION = _init_key_to_section()

explicit_path_keys = {
    "parent_folder",
    "atlas_configuration_file",
    "atlas_settings_file",
    "radiosonde_folder",
    "radiosonde_file",
    "output_folder",
}



def _raise_init_section_and_parameter_errors(config: configparser.ConfigParser) -> None:
    """Validate sections, unknown parameters, and key placement.

    Valid parameter names are defined only by SCHEMA.  INIT_FILE_SECTIONS only
    defines where valid keys are allowed to appear in user-facing INI files.
    The returned caller_info dictionary remains flat and section-independent.
    """

    allowed_sections = set(INIT_FILE_SECTIONS.keys())
    schema_keys = set(SCHEMA.keys())

    unknown_sections: List[str] = []
    unknown_parameters: List[tuple[str, str]] = []
    misplaced_parameters: List[tuple[str, str, str]] = []
    duplicate_parameters: List[tuple[str, str, str]] = []
    seen_parameters: Dict[str, str] = {}

    for section in config.sections():
        if section not in allowed_sections:
            unknown_sections.append(section)
            continue

        for key in config[section].keys():
            key_str = str(key).strip()

            if key_str not in schema_keys:
                unknown_parameters.append((section, key_str))
                continue

            expected_section = INIT_KEY_TO_SECTION[key_str]

            if section != expected_section:
                misplaced_parameters.append((section, key_str, expected_section))

            previous_section = seen_parameters.get(key_str)
            if previous_section is not None:
                duplicate_parameters.append((key_str, previous_section, section))
            else:
                seen_parameters[key_str] = section

    messages: List[str] = []

    if unknown_sections:
        messages.append(
            "Unknown initialization section(s). Section names are strict "
            "lowercase values and must be one of: "
            f"{sorted(allowed_sections)}\n"
            + "\n".join(f"  - [{section}]" for section in unknown_sections)
        )

    if unknown_parameters:
        messages.append(
            "Initialization file contains parameter(s) that are not defined "
            "in the ATLAS initialization schema. Empty declarations are also "
            "invalid.\n"
            + "\n".join(
                f"  - [{section}] {key}"
                for section, key in unknown_parameters
            )
        )

    if misplaced_parameters:
        messages.append(
            "Initialization file contains valid parameter(s) under the wrong "
            "section. Move them to the expected section shown below.\n"
            + "\n".join(
                f"  - [{section}] {key}  -> expected [{expected}]"
                for section, key, expected in misplaced_parameters
            )
        )

    if duplicate_parameters:
        messages.append(
            "Initialization file contains duplicate parameter(s) across "
            "sections. Each parameter may be declared only once.\n"
            + "\n".join(
                f"  - {key}: [{first}] and [{second}]"
                for key, first, second in duplicate_parameters
            )
        )

    if messages:
        raise ConfigError("\n\n".join(messages))


tlc_subfolders = ["north", "east", "south", "west"]
tlc_rin_subfolders = ["inner", "outer"]
pcb_subfolders = ["p45", "m45", "+45", "-45"]

# Legacy physical folder aliases. These are deliberately not SCHEMA entries and
# are not measurement-type identifiers. They only allow existing folders named
# with old suffixes to be detected and renamed to the current physical names
# before the final path registry is built.
legacy_folder_aliases = {
    "abs_nrm": "nrm",
    "abs_drk_nrm": "drk_nrm",
    "abs_nrm_pcb": "nrm_pcb",
    "abs_drk_nrm_pcb": "drk_nrm_pcb",
}


# -------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------

def _is_empty_scalar(v: Any) -> bool:
    return v is None or (isinstance(v, str) and v.strip() == "")


def _ini_root(filepath: str) -> str:
    """Return the absolute folder containing the call_atlas INI file."""

    return os.path.normpath(os.path.dirname(os.path.abspath(filepath)))


def _path_relative_to_ini_root(path: str, filepath: str) -> str:
    """Resolve a path relative to the call_atlas INI folder if needed."""

    path = os.path.expanduser(str(path))

    if not os.path.isabs(path):
        path = os.path.join(_ini_root(filepath), path)

    return os.path.normpath(path)


def _resolve_explicit_paths_relative_to_ini(
    parser_args: Dict[str, Any],
    filepath: str,
) -> Dict[str, Any]:
    """Resolve explicit path entries relative to the call_atlas INI folder.

    This applies only to the user-facing explicit path entries, including
    ``output_folder``.  Relative QA measurement paths such as ``ray`` and
    ``pcb`` keep their existing behavior: they are resolved relative to
    ``parent_folder`` later.
    """

    for key in explicit_path_keys:
        value = parser_args.get(key)

        if value is None:
            continue

        if isinstance(value, str) and value.strip() == "":
            parser_args[key] = None
            continue

        parser_args[key] = _path_relative_to_ini_root(value, filepath)

    return parser_args


def _base_folder_from_parent_folder(parent_folder: str) -> str:
    """Return the folder one level above parent_folder.

    Example
    -------
    /my_drive/station_id/input_data -> /my_drive/station_id
    """

    parent_folder = os.path.normpath(os.path.abspath(os.path.expanduser(str(parent_folder))))
    return os.path.normpath(os.path.dirname(parent_folder))


def _convert_scalar(raw: Optional[str], expected: type, name: str) -> Optional[Union[str, int, float, bool]]:
    """Convert a single value, allowing int-like floats for int. Empty -> None."""
    if raw is None:
        return None

    s = raw.strip()

    if s == "":
        return None

    if expected == int:
        try:
            f = float(s)
        except ValueError:
            raise ConfigError(f"{name}: expected integer, got '{raw}'")

        if f.is_integer():
            return int(f)

        raise ConfigError(f"{name}: expected integer, got non-integer float '{raw}'")

    if expected == float:
        try:
            return float(s)
        except ValueError:
            raise ConfigError(f"{name}: expected float, got '{raw}'")

    if expected == str:
        return s

    if expected == bool:
        if s not in ["True", "False"]:
            raise ConfigError(f"{name}: expected bool, please use either True or False")

        return s == "True"

    raise ConfigError(f"{name}: unsupported dtype {expected}")


def _split_list(raw: Optional[str]) -> List[str]:
    if raw is None:
        return []

    s = raw.strip()

    if s == "":
        return []

    return [x.strip() for x in s.replace(";", ",").split(",") if x.strip() != ""]


def _convert_list(raw: Optional[str], meta: Dict[str, Any], name: str) -> List[Any]:
    items = _split_list(raw)
    out: List[Any] = items

    for i in range(len(items)):
        items[i] = _convert_scalar(items[i], meta["dtype"], name)

    return out


def _check_size(name: str, value: Any, meta: Dict[str, Any]) -> None:
    if meta.get("is_list") and "size" in meta and len(value) > 0:
        if len(value) not in meta["size"]:
            if isinstance(meta["size"], list):
                raise ConfigError(
                    f"{name} must have one of the following sizes: {meta['size']}\n"
                    f"The following values were provided: {value}"
                )
            else:
                raise ConfigError(
                    f"{name} must have size: {meta['size']}\n"
                    f"The following values were provided: {value}"
                )


def _check_limits(name: str, value: Any, meta: Dict[str, Any]) -> None:
    def check_one(v: Any):
        if v is None or isinstance(v, str):
            return

        if "min" in meta and v < meta["min"]:
            raise ConfigError(f"{name}: value {v} below minimum {meta['min']}")

        if "max" in meta and v > meta["max"]:
            raise ConfigError(f"{name}: value {v} above maximum {meta['max']}")

    if meta.get("is_list"):
        if not isinstance(value, list):
            raise ConfigError(
                f"{name}: value {value} is not a list despite being defined "
                "as a list in the SCHEMA"
            )

        if len(value) > 0:
            for v in value or []:
                check_one(v)

    else:
        if isinstance(value, list):
            raise ConfigError(
                f"{name}: value {value} is a list despite being defined "
                "as a scalar parameter in the SCHEMA"
            )

        check_one(value)


def _check_allowed(name: str, value: Any, meta: Dict[str, Any]) -> None:
    if "allowed" not in meta:
        return

    allowed = set(meta["allowed"])

    def check_one(v: Any):
        if v is None:
            return

        if isinstance(v, str) and v.strip() == "":
            return

        if v not in allowed:
            raise ConfigError(f"{name}: value '{v}' not in allowed {sorted(allowed)}")

    if meta.get("is_list"):
        if len(value) > 0:
            for v in value or []:
                check_one(v)
    else:
        check_one(value)


def _warn_recommended(name: str, msg: str, recommended_missing: bool) -> None:
    if recommended_missing:
        print(msg)

    print(f"  --{name}")


# -------------------------------------------------------------------
# Expansion & computed defaults
# -------------------------------------------------------------------

def _fill_with_defaults(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    """Fill empty entries with default values."""

    for key, meta in SCHEMA.items():
        arg = parser_args.get(key)
        default = meta["default"]

        if arg is None or arg == []:
            parser_args[key] = default

    return parser_args


def _expand_slice_exclude_bundles(values: Any) -> Any:
    """Expand slice/exclude bundle aliases to concrete measurement keys.

    Explicit concrete keys always win over bundle aliases within the same
    parameter. For example:

        tlc_inner, 2300, 0100, tlc_rin, 0100, 0200

    becomes:

        tlc_inner, 2300, 0100, tlc_outer, 0100, 0200

    because tlc_inner was explicitly provided, so tlc_rin expands only to the
    remaining ring member.
    """

    if values is None or len(values) == 0:
        return values

    # Let the existing detailed validation raise the triplet-format error.
    if len(values) % 3 != 0:
        return values

    explicit_keys = {
        str(values[i]).strip()
        for i in range(0, len(values), 3)
        if str(values[i]).strip() not in slice_exclude_bundles
    }

    expanded: List[Any] = []

    for i in range(0, len(values), 3):
        key = str(values[i]).strip()
        start_time = values[i + 1]
        stop_time = values[i + 2]

        if key in slice_exclude_bundles:
            keys = [
                member
                for member in slice_exclude_bundles[key]
                if member not in explicit_keys
            ]
        else:
            keys = [key]

        for expanded_key in keys:
            expanded.extend([expanded_key, start_time, stop_time])

    return expanded


def _compute_emitted_wavelength_if_missing(parser_args: Dict[str, Any]) -> Dict[str, Any]:

    det = parser_args.get("detected_wavelength")
    em = parser_args.get("emitted_wavelength")

    if len(det) > 0:
        if len(em) == 0 or all(v is None for v in em):
            for i in range(len(det)):
                if det[i] is None:
                    parser_args["emitted_wavelength"][i] = None
                elif det[i] < 520:
                    parser_args["emitted_wavelength"][i] = 354.71
                elif det[i] < 1000:
                    parser_args["emitted_wavelength"][i] = 532.07
                else:
                    parser_args["emitted_wavelength"][i] = 1064.14

    return parser_args


def _special_checks(parser_args: Dict[str, Any]) -> None:

    ray_averaging_time = parser_args.get("ray_averaging_time")
    ray_qck_averaging_time = parser_args.get("ray_qck_averaging_time")

    pattern = re.compile(r"^\d{1,2}(?:min|h)$")

    def is_valid_time_format(s: str) -> bool:
        return bool(pattern.match(s))

    if ray_averaging_time is not None:
        if not is_valid_time_format(ray_averaging_time):
            raise ConfigError(
                f"The provided ray_averaging_time format was not understood: "
                f"{ray_averaging_time} Allowed formats: 'xmin' or 'xH' where x "
                "is an up to 2 digit positive integer that corresponds either "
                "the number of integers or the number of hours. For exacmple: "
                "10min or 3H"
            )
    else:
        parser_args["ray_averaging_time"] = "all"

    if ray_averaging_time is not None:
        if not is_valid_time_format(ray_qck_averaging_time):
            raise ConfigError(
                f"The provided ray_qck_averaging_time format was not understood: "
                f"{ray_qck_averaging_time} Allowed formats: 'xmin' or 'xH' where x "
                "is an up to 2 digit positive integer that corresponds either "
                "the number of integers or the number of hours. For exacmple: "
                "10min or 3H"
            )
    else:
        parser_args["ray_qck_averaging_time"] = "raw"

    rsonde_column_units = parser_args.get("rsonde_column_units")

    if rsonde_column_units[0] not in height_units:
        raise ConfigError(
            f"The height units provided in rsonde_column_units are wrong "
            f"{rsonde_column_units} Allowed values: {height_units}"
        )

    if rsonde_column_units[1] not in pressure_units:
        raise ConfigError(
            f"The pressure units provided in rsonde_column_units are wrong "
            f"{rsonde_column_units} Allowed values: {pressure_units}"
        )

    if rsonde_column_units[2] not in temperature_units:
        raise ConfigError(
            f"The temperature units provided in rsonde_column_units are wrong "
            f"{rsonde_column_units} Allowed values: {temperature_units}"
        )

    if rsonde_column_units[3] not in humidity_units:
        raise ConfigError(
            f"The humidity units provided in rsonde_column_units are wrong "
            f"{rsonde_column_units} Allowed values: {humidity_units}"
        )

    # Expand user-facing bundle aliases before validation so downstream code
    # receives only concrete measurement keys.
    parser_args["slice_measurement"] = _expand_slice_exclude_bundles(
        parser_args.get("slice_measurement")
    )
    parser_args["exclude_measurement"] = _expand_slice_exclude_bundles(
        parser_args.get("exclude_measurement")
    )

    slice_measurement = parser_args.get("slice_measurement")
    exclude_measurement = parser_args.get("exclude_measurement")

    accepted_time_formats = (
        "HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, "
        "or yyyymmdd_HHMMSS"
    )

    def wrong_format(name: str, identifier: str, value: str) -> str:
        return (
            f"The format of the provided {name} {identifier} {value!r} is wrong. "
            f"Accepted time formats are: {accepted_time_formats}."
        )

    def _valid_time_parts(hour: str, minute: str = "00", second: str = "00") -> bool:
        try:
            h = int(hour)
            m = int(minute)
            s = int(second)
        except Exception:
            return False

        return 0 <= h <= 23 and 0 <= m <= 59 and 0 <= s <= 59

    def _is_valid_slice_time(value: Any) -> bool:
        """Validate slice/exclude time strings.

        Accepted formats:
        - HHMM
        - yyyymmdd
        - yyyymmdd_HH
        - yyyymmdd_HHMM
        - yyyymmdd_HHMMSS
        """

        if value is None:
            return False

        text = str(value).strip()

        # Legacy HHMM format. The measurement date is assigned later in
        # compute_slice_and_exclude().
        if re.fullmatch(r"\d{4}", text):
            return _valid_time_parts(text[:2], text[2:4])

        # Absolute date only: yyyymmdd -> midnight.
        match = re.fullmatch(r"(\d{8})", text)
        if match:
            try:
                datetime.strptime(match.group(1), "%Y%m%d")
                return True
            except ValueError:
                return False

        # Absolute date/time: yyyymmdd_HH[MM[SS]].
        match = re.fullmatch(r"(\d{8})_(\d{2}|\d{4}|\d{6})", text)
        if not match:
            return False

        date_part, time_part = match.groups()

        try:
            datetime.strptime(date_part, "%Y%m%d")
        except ValueError:
            return False

        hour = time_part[:2]
        minute = time_part[2:4] if len(time_part) >= 4 else "00"
        second = time_part[4:6] if len(time_part) == 6 else "00"

        return _valid_time_parts(hour, minute, second)

    example_text = (
        "Only sets of 3 strings are acceptable: test identifier, start time, "
        "and stop time. Accepted time formats are: "
        f"{accepted_time_formats}. Examples: ray, 2330, 0100 or "
        "ray, 20251203_2330, 20251204_0100."
    )

    def _check_slice_exclude_parameter(name: str, values: Any) -> None:
        if values is None or len(values) == 0:
            return

        if len(values) % 3 != 0:
            raise ConfigError(
                f"The provided {name} parameter is wrong {values}\n{example_text}"
            )

        for i in range(0, len(values), 3):
            qa_key = values[i]
            start_time = values[i + 1]
            stop_time = values[i + 2]

            if qa_key not in qa_measurement_folders:
                raise ConfigError(
                    f"{name}: unrecognized qa_test ID provided: "
                    f"{qa_key}\nRecognised IDs: {qa_measurement_folders}"
                )

            if not _is_valid_slice_time(start_time):
                raise ConfigError(wrong_format(name, "start time", start_time))

            if not _is_valid_slice_time(stop_time):
                raise ConfigError(wrong_format(name, "stop time", stop_time))

    _check_slice_exclude_parameter("slice_measurement", slice_measurement)
    _check_slice_exclude_parameter("exclude_measurement", exclude_measurement)

    CustomWarning(
        "Please note that since update 0.6.0 the more generic slice_measurement "
        "parameter can be provided instead of slice_rayleigh. The slice_rayleigh "
        "parameter will be ignored even if provided"
    )


# -------------------------------------------------------------------
# Presence checks
# -------------------------------------------------------------------

def _enforce_mandatory_and_recommended(parser_args: Dict[str, Any]) -> None:
    recommended_first_time = True

    for name, meta in SCHEMA.items():
        cat = meta.get("category", "optional")
        val = parser_args.get(name)

        if meta.get("is_list", False):
            is_empty = (
                val is None
                or len(val) == 0
                or all(_is_empty_scalar(v) for v in val)
            )
        else:
            is_empty = _is_empty_scalar(val)

        if cat == "mandatory" and is_empty:
            raise ConfigError(f"{name} is mandatory and was not provided.")

        if cat == "recommended" and is_empty:
            _warn_recommended(
                name,
                "--Warning: Recomended configuration parameters not provided:",
                recommended_first_time,
            )
            recommended_first_time = False


def _resolve_default_base_paths(parser_args: Dict[str, Any], filepath: str) -> Dict[str, Any]:
    """
    Resolve default folders.

    Explicit paths may be absolute or relative to the folder containing the
    call_atlas INI file.  The old user-facing ``main_data_folder`` parameter has
    been retired.  Internally, keep ``main_data_folder`` equal to the INI folder
    for compatibility with downstream helpers.

    Defaults
    --------
    parent_folder:
        <ini_folder>/input_data

    radiosonde_folder:
        <base_folder>/radiosondes

    where <base_folder> is the folder one level above parent_folder.  For
    example, if parent_folder is /my_drive/station_id/input_data, then
    radiosonde_folder becomes /my_drive/station_id/radiosondes.
    """

    init_folder = _ini_root(filepath)

    # Retained only as an internal compatibility value for autodetect_paths and
    # downstream code that may still expect this key in caller_info.
    parser_args["main_data_folder"] = init_folder

    if parser_args.get("parent_folder") is None:
        parser_args["parent_folder"] = os.path.join(init_folder, "input_data")
    else:
        parser_args["parent_folder"] = os.path.normpath(parser_args["parent_folder"])

    base_folder = _base_folder_from_parent_folder(parser_args["parent_folder"])

    if parser_args.get("radiosonde_folder") is None and parser_args.get("radiosonde_file") is None:
        parser_args["radiosonde_folder"] = os.path.join(base_folder, "radiosondes")

    return parser_args

def _resolve_default_configuration_file(
    parser_args: Dict[str, Any],
    filepath: str,
    *,
    require_existing_local_config: bool = True,
    resolve_remote_target: bool = True,
) -> Dict[str, Any]:
    """
    Resolve the default ATLAS configuration file.

    The default configuration folder is placed next to parent_folder, not next
    to the call_atlas INI file.  For example:

        parent_folder = /my_drive/station_id/input_data
        cfg_folder    = /my_drive/station_id/configurations

    For export_hoi_cfg == "0", the local default must already exist:
        <cfg_folder>/config_file.ini

    For export_hoi_cfg == "1" or "2", the configuration file is expected to be
    downloaded/exported later, so only the target filename is prepared here:
        <cfg_folder>/config_file_{scc_configuration_id}.ini
    """

    export_hoi_cfg = parser_args.get("export_hoi_cfg")
    atlas_configuration_file = parser_args.get("atlas_configuration_file")

    parent_folder = parser_args.get("parent_folder")
    if parent_folder is None:
        parent_folder = os.path.join(_ini_root(filepath), "input_data")

    base_folder = _base_folder_from_parent_folder(parent_folder)
    cfg_folder = os.path.join(base_folder, "configurations")

    if export_hoi_cfg == "0" and atlas_configuration_file is None:
        default_cfg = os.path.join(cfg_folder, "config_file.ini")

        if require_existing_local_config and not os.path.isfile(default_cfg):
            raise ConfigError(
                "atlas_configuration_file was not provided and export_hoi_cfg "
                "is 0, so the default configuration file was expected at:\n"
                f"{default_cfg}"
            )

        parser_args["atlas_configuration_file"] = default_cfg

    elif (
        resolve_remote_target
        and export_hoi_cfg in ["1", "2"]
        and atlas_configuration_file is None
    ):
        scc_configuration_id = parser_args.get("scc_configuration_id")

        if scc_configuration_id is None or str(scc_configuration_id).strip() == "":
            raise ConfigError(
                f"export_hoi_cfg was set to {export_hoi_cfg}, so "
                "scc_configuration_id must be provided to download/export the "
                "SCC configuration file."
            )

        if os.path.exists(cfg_folder) and not os.path.isdir(cfg_folder):
            raise ConfigError(
                "The configurations path needed for atlas_configuration_file "
                f"exists but is not a directory: {cfg_folder}"
            )

        try:
            os.makedirs(cfg_folder, exist_ok=True)
        except Exception as exc:
            raise ConfigError(
                "The configurations folder needed for atlas_configuration_file "
                f"could not be created: {cfg_folder}"
            ) from exc

        parser_args["atlas_configuration_file"] = os.path.join(
            cfg_folder,
            f"config_file_{scc_configuration_id}.ini",
        )

    return parser_args

def _absolute_paths_exist_check(parser_args: Dict[str, Any]) -> Dict[str, Any]:

    for key, meta in SCHEMA.items():
        path = parser_args.get(key)

        if meta.get("check_path") in ["dir", "file"] and path is not None:
            path = os.path.normpath(os.path.expanduser(str(path)))

            if key == "radiosonde_folder":
                if os.path.exists(path) and not os.path.isdir(path):
                    raise ConfigError(
                        f"{key} is provided but points to an existing file, not a directory: {path}"
                    )

                try:
                    os.makedirs(path, exist_ok=True)
                except Exception as exc:
                    raise ConfigError(
                        f"{key} does not point to an existing path and could not be created: {path}"
                    ) from exc

                parser_args[key] = path
                continue

            if not os.path.exists(path):
                raise ConfigError(
                    f"{key} does not point to an existing path: {path}"
                )

            if meta.get("check_path") == "dir":
                if not os.path.isdir(path):
                    raise ConfigError(
                        f"{key} is provided but does not point to an existing directory: {path}"
                    )

            if meta.get("check_path") == "file":
                if not os.path.isfile(path):
                    raise ConfigError(
                        f"{key} is provided but does not point to an existing file: {path}"
                    )

            parser_args[key] = path

    return parser_args

def _relative_paths_exist_check(parser_args: Dict[str, Any]) -> Dict[str, Any]:

    parent_folder = parser_args.get("parent_folder")

    if parent_folder is None:
        raise ConfigError(
            "parent_folder is empty. It should have been already assigned at this stage"
        )

    for key, meta in SCHEMA.items():
        rel_path = parser_args.get(key)

        if meta.get("check_path") == "relative":
            parser_args[f"abs_{key}"] = None

            if key != "drk":
                parser_args[f"abs_drk_{key}"] = None

            if rel_path is not None:
                path = os.path.normpath(os.path.join(parent_folder, rel_path))
                path_drk = os.path.normpath(os.path.join(parent_folder, f"drk_{rel_path}"))

                if not os.path.exists(path):
                    raise ConfigError(
                        f"{key} is provided but does not point to an existing path: {path}"
                    )

                if not os.path.isdir(path):
                    raise ConfigError(
                        f"{key} is provided but does not point to an existing directory: {path}"
                    )

                parser_args[f"abs_{key}"] = path

                if key == "drk":
                    continue

                if os.path.exists(path_drk):
                    parser_args[f"abs_drk_{key}"] = path_drk

            else:
                path = os.path.normpath(os.path.join(parent_folder, key))
                path_drk = os.path.normpath(os.path.join(parent_folder, f"drk_{key}"))

                if os.path.exists(path):
                    parser_args[f"abs_{key}"] = path

                if key == "drk":
                    continue

                if os.path.exists(path_drk):
                    parser_args[f"abs_drk_{key}"] = path_drk

    return parser_args


def _handle_telecover_path_aliases(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle telecover path aliases without deleting active tlc paths.

    Since the internal telecover quadrants were renamed from ``tlc_qua`` to
    ``tlc``, the active base keys are now ``abs_tlc`` and ``abs_drk_tlc``.
    They must be preserved so that ``_special_path_handling`` can expand them
    into ``abs_tlc_north/east/south/west``.

    The old implementation mapped ``abs_tlc`` to ``abs_tlc_qua`` and then
    removed ``abs_tlc``.  That was correct only while ``tlc_qua`` was the
    internal quadrant key.
    """

    return parser_args


def _detect_legacy_folder_aliases(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    """Detect old physical folder names without making them parser keys.

    The INI schema intentionally exposes only the current logical folder names
    (for example ``ray``).  This helper only looks on disk for old folder names
    such as ``nrm`` so that ``_rename_folder`` can rename them to the current
    names before downstream path handling starts.
    """

    parent_folder = parser_args.get("parent_folder")

    if parent_folder is None:
        return parser_args

    parent = Path(parent_folder)

    for abs_key, folder_name in legacy_folder_aliases.items():
        parser_args.setdefault(abs_key, None)
        candidate = parent / folder_name

        if candidate.exists():
            if not candidate.is_dir():
                raise ConfigError(
                    f"Legacy folder alias {folder_name!r} exists but is not a directory: "
                    f"{candidate}"
                )

            parser_args[abs_key] = os.path.normpath(str(candidate))

    return parser_args


def _rename_folder(
    parser_args: Dict[str, Any],
    old_key: str,
    new_key: str,
    new_name: str,
    version_warning: str = "",
) -> Dict[str, Any]:

    old_value = parser_args.get(old_key)
    new_value = parser_args.get(new_key)

    if old_value is not None and new_value is not None:
        raise ConfigError(
            f"{old_key} and {new_key} folders cannot be present at the same time. "
            f"{version_warning}"
        )

    elif old_value is not None and new_value is None:
        CustomWarning(f"folder with the {old_key} suffix was detected. {version_warning}")

        old = Path(old_value)
        new = old.parent / new_name

        old.rename(new)

        parser_args[new_key] = os.path.normpath(str(new))

    parser_args.pop(old_key, None)

    return parser_args


def _special_path_handling(parser_args: Dict[str, Any]) -> Dict[str, Any]:

    parser_args = _detect_legacy_folder_aliases(parser_args)

    # Rename any existing folders with the nrm suffix to ray.
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_nrm",
        new_key="abs_ray",
        new_name="ray",
        version_warning=version_warning_nrm,
    )
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_drk_nrm",
        new_key="abs_drk_ray",
        new_name="drk_ray",
        version_warning=version_warning_nrm,
    )
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_nrm_pcb",
        new_key="abs_ray_pcb",
        new_name="ray_pcb",
        version_warning=version_warning_nrm,
    )
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_drk_nrm_pcb",
        new_key="abs_drk_ray_pcb",
        new_name="drk_ray_pcb",
        version_warning=version_warning_nrm,
    )

    # ray and ray_pcb measurements will be the same if only one of them is provided.
    path_drk = parser_args.get("abs_drk")
    path_ray = parser_args.get("abs_ray")
    path_ray_pcb = parser_args.get("abs_ray_pcb")
    path_pcb = parser_args.get("abs_pcb")

    if path_ray is None and path_ray_pcb is not None:
        parser_args["abs_ray"] = path_ray_pcb
        parser_args["abs_drk_ray"] = parser_args.get("abs_drk_ray_pcb")

    if path_ray_pcb is None and path_ray is not None:
        parser_args["abs_ray_pcb"] = path_ray
        parser_args["abs_drk_ray_pcb"] = parser_args.get("abs_drk_ray")

    # Use common drk measurement if not provided explicitly for a QA test.
    for key in qa_tests:
        if key == "drk":
            continue

        if (
            parser_args.get(f"abs_drk_{key}") is None
            and path_drk is not None
            and parser_args.get(f"abs_{key}") is not None
        ):
            parser_args[f"abs_drk_{key}"] = path_drk

    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")

    parent_folder = os.path.normpath(parser_args["parent_folder"])
    parent_folder_name = os.path.basename(parent_folder)

    if parser_args.get("output_folder") is None:
        base_folder = parser_args.get("main_data_folder")
        if base_folder is None:
            base_folder = os.path.dirname(parent_folder)

        parser_args["output_folder"] = os.path.join(
            base_folder,
            "analysis",
            parent_folder_name,
        )

    else:
        parser_args["output_folder"] = os.path.join(
            parser_args["output_folder"],
            parent_folder_name,
        )

    os.makedirs(parser_args["output_folder"], exist_ok=True)

    if parser_args.get("overwrite_output"):
        parser_args["plot_folder"] = os.path.join(parser_args["output_folder"], "plots")
        parser_args["ascii_folder"] = os.path.join(parser_args["output_folder"], "ascii")
    else:
        parser_args["plot_folder"] = os.path.join(
            parser_args["output_folder"],
            f"plots_{timestamp}",
        )
        parser_args["ascii_folder"] = os.path.join(
            parser_args["output_folder"],
            f"ascii_{timestamp}",
        )

    os.makedirs(parser_args["plot_folder"], exist_ok=True)
    os.makedirs(parser_args["ascii_folder"], exist_ok=True)

    # Create the tlc, tlc_rin, pcb, and pcb_aux related subfolders.
    # Note: abs_tlc may point to a physical folder named 'tlc'.
    _expand_subfolders(
        parser_args,
        base_key="abs_tlc",
        subfolders=tlc_subfolders,
    )
    _expand_subfolders(
        parser_args,
        base_key="abs_tlc_rin",
        subfolders=tlc_rin_subfolders,
    )
    _expand_subfolders(
        parser_args,
        base_key="abs_pcb",
        subfolders=pcb_subfolders,
    )
    _expand_subfolders(
        parser_args,
        base_key="abs_pcb_aux",
        subfolders=pcb_subfolders,
    )

    # Rename any existing +45 and -45 folders to p45 and m45.
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_pcb_+45",
        new_key="abs_pcb_p45",
        new_name="p45",
        version_warning=version_warning_pcb,
    )
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_pcb_-45",
        new_key="abs_pcb_m45",
        new_name="m45",
        version_warning=version_warning_pcb,
    )
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_pcb_aux_+45",
        new_key="abs_pcb_aux_p45",
        new_name="p45",
        version_warning=version_warning_pcb,
    )
    parser_args = _rename_folder(
        parser_args,
        old_key="abs_pcb_aux_-45",
        new_key="abs_pcb_aux_m45",
        new_name="m45",
        version_warning=version_warning_pcb,
    )

    all_none = all(
        parser_args[key] is None
        for key in parser_args.keys()
        if key.startswith("abs")
    )

    if all_none:
        endpoint(1)

    return parser_args

def _canonicalize_telecover_sector_keys(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert expanded telecover sector/ring path keys to loader-friendly names.

    Final telecover key order:
        abs_tlc_north
        abs_tlc_east
        abs_tlc_south
        abs_tlc_west
        abs_tlc_outer
        abs_tlc_inner

    This prevents mixed data_pack keys such as:
        tlc_south
        tlc_rin_inner
    """

    replacements = [
        # Quadrant telecover keys are already loader-friendly after the
        # tlc_qua -> tlc rename, so they must not be self-renamed/popped.

        # Order intentionally outer before inner.
        ("abs_tlc_rin_outer", "abs_tlc_outer"),
        ("abs_tlc_rin_inner", "abs_tlc_inner"),

        # Order intentionally outer before inner.
        ("abs_drk_tlc_rin_outer", "abs_drk_tlc_outer"),
        ("abs_drk_tlc_rin_inner", "abs_drk_tlc_inner"),
    ]

    for old_key, new_key in replacements:
        if old_key not in parser_args:
            continue

        # Protect against accidental self-renames.  Popping in that case would
        # delete a valid key such as abs_tlc_north.
        if old_key == new_key:
            continue

        old_val = parser_args.get(old_key)

        if parser_args.get(new_key) is None:
            parser_args[new_key] = old_val
        elif old_val is not None and parser_args[new_key] != old_val:
            raise ConfigError(
                f"Conflicting telecover paths detected for {new_key}: "
                f"{parser_args[new_key]} and {old_val}"
            )

        parser_args.pop(old_key, None)

    return parser_args

def _expand_subfolders(
    d: dict,
    base_key: str,
    subfolders: Iterable[str] | Mapping[str, str],
    *,
    must_exist: bool = True,
    normalize_to_path: bool = False,
    delete_base: bool = True,
) -> List[str]:
    """
    Expand a base path key into multiple subfolder keys.
    - If base is None, still create subkeys but set them to None.
    - If must_exist=True, subkeys pointing to non-existent folders are set to None.
    - If must_exist=False, subkeys are created regardless of existence.
    - If delete_base=True, remove the original base_key after expansion.
    Returns the list of created subkeys.
    """

    created: List[str] = []
    base = d.get(base_key)

    items = (
        subfolders.items()
        if isinstance(subfolders, Mapping)
        else ((name, name) for name in subfolders)
    )

    for suffix, foldername in items:
        new_key = f"{base_key}_{suffix}"

        if base is None:
            d[new_key] = None
        else:
            p = Path(base) / foldername

            if must_exist and not p.exists():
                d[new_key] = None
            else:
                d[new_key] = p if normalize_to_path else str(p)

        created.append(new_key)

    if delete_base:
        d.pop(base_key, None)

    return created


class DistributionError(Exception):
    """Raised when the file distribution preconditions are not met."""


def _distribute_files_into_sectors(
    cfg: dict,
    *,
    base_key: str,
    files_per_key: str,
    sectors: Sequence[str],
    pattern: str = "*",
    overwrite: bool = False,
) -> None:
    """
    If cfg[files_per_key] is not None, distribute files from cfg[base_key] into subfolders
    named <sector> (e.g., north/east/south/west) in sequential blocks of size files_per_quadrant.
    Enforces:
      - base folder must exist when files_per_quadrant is not None
      - base folder must not be empty when files_per_quadrant is not None
      - total files must be divisible by files_per_quadrant * len(sectors)
    Moves files in the order of sorted names to keep behavior deterministic.
    """

    files_per_quadrant = cfg.get(files_per_key, None)

    if files_per_quadrant is None:
        return

    base_val = cfg.get(base_key, None)

    if base_val is None:
        raise DistributionError(
            f"{base_key} is None but {files_per_key} is set to {files_per_quadrant}."
        )

    base = Path(base_val)

    if not base.exists() or not base.is_dir():
        raise DistributionError(f"{base_key} points to a non-existent directory: {base}")

    files = sorted([p for p in base.glob(pattern) if p.is_file() and p.parent == base and 'temp' not in p.name])

    if not files:
        CustomWarning(
            f"No files detected in {base} (pattern='{pattern}') while "
            f"{files_per_key}={files_per_quadrant}."
        )
        return

    k = len(sectors)
    group = files_per_quadrant * k

    if len(files) % group != 0:
        raise DistributionError(
            f"The {files_per_key} was provided but the file count {len(files)} "
            f"in {base} is not divisible by {files_per_key} * folders "
            f"({files_per_quadrant} * {k} = {group})."
        )

    dest_dirs = []

    for sect in sectors:
        d = base / sect
        d.mkdir(parents=True, exist_ok=True)
        dest_dirs.append(d)

    for i, f in enumerate(files):
        sector_idx = (i // files_per_quadrant) % k
        dest_dir = dest_dirs[sector_idx]
        target = dest_dir / f.name

        if overwrite:
            f.replace(target)
        else:
            if target.exists():
                raise DistributionError(f"Destination already has file: {target}")

            f.rename(target)

    for sect, d in zip(sectors, dest_dirs):
        cfg[f"{base_key}_{sect}"] = d


def assert_pairs_unique(recorder_channel_id, laser_id):
    if len(recorder_channel_id) != len(laser_id):
        raise ConfigError(
            "Duplicate recorder_channel_id values found {recorder_channel_id}. "
            "laser_id is mandatory in this case"
        )

    seen = set()
    dups = []

    for i, pair in enumerate(zip(recorder_channel_id, laser_id)):
        if pair in seen:
            dups.append((i, pair))
        else:
            seen.add(pair)

    if dups:
        raise ConfigError(
            f"Duplicate recorder_channel_id values found {recorder_channel_id}. "
            f"laser_id was provided {laser_id} but it is the same for at least "
            "one of the duplicate recorder_channel_id values"
        )


def read_ini_file(filepath: str) -> Dict[str, Any]:

    """Read, convert, expand from scalar defaults, compute simple defaults, validate, and return dict."""

    config = configparser.ConfigParser(allow_no_value=True, strict=True)
    config.optionxform = str

    read_files = config.read(filepath, encoding="utf-8")

    if not read_files:
        raise ConfigError(
            f"INI file not found or unreadable: {filepath}\n"
            "Make sure the encoding is utf-8"
        )

    _raise_init_section_and_parameter_errors(config)

    parser_args: Dict[str, Any] = {}

    for key, meta in SCHEMA.items():
        found = False

        for section in config.sections():
            if key in config[section]:
                found = True
                raw = config[section][key]

                if meta["is_list"]:
                    parser_args[key] = _convert_list(raw, meta, key)
                else:
                    parser_args[key] = _convert_scalar(raw, meta["dtype"], key)

        if not found:
            if meta["is_list"]:
                parser_args[key] = []
            else:
                parser_args[key] = meta["default"]

    return parser_args


def _get_mtype(d: Dict[str, Any]) -> Dict[str, Any]:

    to_add = {}

    tlc_parts = ["north", "east", "south", "west"]
    tlc_rin_parts = ["outer", "inner"]

    for key in d.keys():
        name = key.removeprefix("abs_")

        if key.startswith("abs_drk"):
            to_add[f"mtype_{name}"] = "drk"

        elif key.startswith("abs_ray"):
            to_add[f"mtype_{name}"] = "nrm"

        elif key.startswith("abs_tlc_"):
            suffix = name.removeprefix("tlc_")

            if suffix in tlc_parts:
                to_add[f"mtype_{name}"] = "tlc"
            elif suffix in tlc_rin_parts:
                to_add[f"mtype_{name}"] = "tlc_rin"

        elif key.startswith("abs_pcb_p45"):
            to_add[f"mtype_{name}"] = "pcb_p45"

        elif key.startswith("abs_pcb_m45"):
            to_add[f"mtype_{name}"] = "pcb_m45"

        elif key.startswith("abs_trg"):
            to_add[f"mtype_{name}"] = "nrm"

        elif key.startswith("abs_dtm"):
            to_add[f"mtype_{name}"] = "nrm"

        elif key.startswith("abs_cam"):
            to_add[f"mtype_{name}"] = "cam"

    d.update(to_add)

    return d


def measurement_type(meas_key: str) -> str:
    """
    Infer the raw-reader measurement type from a caller_info measurement key.

    This replaces the old flat mtype_* entries in caller_info. The returned
    value is still the value expected by the raw file readers.
    """

    if meas_key.startswith("drk"):
        return "drk"

    if meas_key.startswith("ray"):
        return "nrm"

    if meas_key in ["trg", "dtm"]:
        return "nrm"

    if meas_key in ["tlc_north", "tlc_east", "tlc_south", "tlc_west"]:
        return "tlc"

    if meas_key in ["tlc_inner", "tlc_outer"]:
        return "tlc_rin"

    if meas_key == "pcb_p45":
        return "pcb_p45"

    if meas_key == "pcb_m45":
        return "pcb_m45"

    if meas_key == "pcb_aux_p45":
        return "pcb_p45"

    if meas_key == "pcb_aux_m45":
        return "pcb_m45"

    if meas_key.startswith("cam"):
        return "cam"

    raise ConfigError(f"Could not infer measurement type for {meas_key}")


def _build_path_registry(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert temporary flat abs_* keys into a nested path registry.

    Final caller_info entries:
        paths:
            Dict with keys from qa_measurement_folders and values the physical
            folders that should be read exactly once.

        loading_map:
            Dict with logical measurement keys as keys and the already-loaded
            physical measurement keys as values.

    Examples:
        loading_map["ray_pcb"] = "ray"
        loading_map["drk_ray"] = "drk"
        loading_map["drk_tlc"] = "drk"
    """

    path_by_key: Dict[str, str] = {}

    for meas_key in qa_measurement_folders:
        abs_key = f"abs_{meas_key}"
        path = parser_args.get(abs_key)

        if path is not None:
            path_by_key[meas_key] = os.path.normpath(str(path))

    paths: Dict[str, str] = {}
    loading_map: Dict[str, str] = {}
    owner_by_path: Dict[str, str] = {}

    for meas_key in qa_measurement_folders:
        path = path_by_key.get(meas_key)

        if path is None:
            continue

        if path not in owner_by_path:
            owner_by_path[path] = meas_key
            paths[meas_key] = path
        else:
            loading_map[meas_key] = owner_by_path[path]

    parser_args["paths"] = dict(sorted(paths.items()))
    parser_args["loading_map"] = dict(sorted(loading_map.items()))

    # Remove temporary flat absolute-path keys from caller_info.
    for key in list(parser_args.keys()):
        if key.startswith("abs_"):
            parser_args.pop(key)

    # Remove legacy mtype_* keys if they were added by older code paths.
    for key in list(parser_args.keys()):
        if key.startswith("mtype_"):
            parser_args.pop(key)

    return parser_args

def _export_hoi_cfg_check(parser_args: Dict[str, Any]) -> None:
    export_hoi_cfg = parser_args["export_hoi_cfg"]
    atlas_configuration_file = parser_args["atlas_configuration_file"]

    if export_hoi_cfg in ["1", "2"] and atlas_configuration_file is not None:
        raise ConfigError(
            f"export_hoi_cfg was set to {export_hoi_cfg} but "
            "atlas_configuration_file was provided in the explicit_paths section. "
            "Either change export_hoi_cfg to 0 or leave atlas_configuration_file empty. "
        )


def _scc_configuration_id_check(parser_args: Dict[str, Any]) -> None:
    export_hoi_cfg = parser_args["export_hoi_cfg"]
    scc_configuration_id = parser_args.get("scc_configuration_id")

    if export_hoi_cfg in ["1", "2"] and (
        scc_configuration_id is None or str(scc_configuration_id).strip() == ""
    ):
        raise ConfigError(
            f"export_hoi_cfg was set to {export_hoi_cfg}, so "
            "scc_configuration_id must be provided to download/export the SCC "
            "configuration file."
        )


# -------------------------------------------------------------------
# Public API
# -------------------------------------------------------------------

def parse_call_atlas_ini(filepath: str, debug: bool = False) -> Dict[str, Any]:

    print_header(f"Parsing the ATLAS initialization file\n{filepath}")

    # 1) Parse & convert
    parser_args = read_ini_file(filepath)

    # 2) Enforce mandatory/recommended
    _enforce_mandatory_and_recommended(parser_args)

    # 2b) Resolve explicit absolute/relative path entries.  Relative paths are
    # interpreted with respect to the folder containing this initialization file.
    parser_args = _resolve_explicit_paths_relative_to_ini(parser_args, filepath)

    # 3) Resolve default base folders.
    # The old user-facing main_data_folder parameter is no longer read from
    # the INI, but an internal compatibility value is still populated.
    parser_args = _resolve_default_base_paths(parser_args, filepath)

    # 4) Fill with default values if empty
    parser_args = _fill_with_defaults(parser_args)

    # 5) Keep the existing HOI export conflict check. For export_hoi_cfg 1 or 2,
    # atlas_configuration_file must remain empty if it was provided explicitly.
    _export_hoi_cfg_check(parser_args)

    # 5b) For HOI export/download modes, scc_configuration_id is required.
    _scc_configuration_id_check(parser_args)

    # 5c) If no configuration file was provided and local config export is used,
    # fall back to <parent_folder_parent>/configurations/config_file.ini.
    parser_args = _resolve_default_configuration_file(
        parser_args,
        filepath,
        resolve_remote_target=False,
    )

    # 6) Check if all provided paths exist
    parser_args = _absolute_paths_exist_check(parser_args)

    # 6b) For HOI export/download modes, prepare the target configuration
    # filename after the file-existence checks because the file may not exist yet.
    parser_args = _resolve_default_configuration_file(
        parser_args,
        filepath,
        require_existing_local_config=False,
    )

    # # 7) Fill in the explicit paths if they are not available with the autodetected paths
    # parser_args = autodetect_paths(parser_args)

    # 8) Check if the relative QA test folder paths exist and add corresponding absolute paths
    parser_args = _relative_paths_exist_check(parser_args)

    # 8b) Map physical folder/alias 'tlc' to internal key 'tlc'
    parser_args = _handle_telecover_path_aliases(parser_args)

    # 9) Distribute tlc files in the correct sector/ring if requested
    _distribute_files_into_sectors(
        parser_args,
        base_key="abs_tlc",
        files_per_key="files_per_quadrant",
        sectors=tlc_subfolders,
    )

    _distribute_files_into_sectors(
        parser_args,
        base_key="abs_tlc_rin",
        files_per_key="files_per_ring",
        sectors=tlc_rin_subfolders,
    )

    # 10) Fill absolute paths from provided relative paths and existing folders
    parser_args = _special_path_handling(parser_args)
    
    # 10b) Rename expanded telecover keys to loader-friendly names.
    # This prevents mixed keys like tlc_north and tlc_north.
    parser_args = _canonicalize_telecover_sector_keys(parser_args)

    # 11) Collapse temporary abs_* keys into caller_info["paths"] and
    # caller_info["loading_map"].
    parser_args = _build_path_registry(parser_args)

    # 12) Range & allowed checks, only on non-empty values
    for name, meta in SCHEMA.items():
        _check_size(name, parser_args.get(name), meta)
        _check_limits(name, parser_args.get(name), meta)
        _check_allowed(name, parser_args.get(name), meta)

    # 13) Special handling - check the format of slice_measurement and exclude_measurement
    _special_checks(parser_args)

    # 14) Sort keys alphabetically
    parser_args = dict(sorted(parser_args.items()))

    if debug:
        pprint(parser_args)

    return parser_args

