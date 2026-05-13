#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:16:34 2025

@author: nikos
"""

from __future__ import annotations
import configparser
import warnings
from typing import Any, Dict, List, Optional, Union
import numpy as np
from helper_functions.printouts import print_header
from pprint import pprint
from utils.error_classes import ConfigError, CustomWarning

Number = Union[int, float]

# -------------------------------------------------------------------
# SCHEMA
#   dtype:        expected Python type (str, int, float)
#   default:      scalar default (always scalar; expanded to lists as needed)
#   is_list:      True if value is a list in INI (comma or semicolon-separated)
#   allowed:      optional list of allowed values (checked only on non-empty values)
#   category:     "mandatory" | "recommended" | "optional"
#   min/max:      optional numeric bounds
# -------------------------------------------------------------------
SCHEMA: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "station_id":         {"dtype": str,   "default": None, "is_list": False, "category": "mandatory"},
    "lidar_name":         {"dtype": str,   "default": None, "is_list": False, "category": "optional"},
    "station_name":       {"dtype": str,   "default": None, "is_list": False, "category": "optional"},

    "lidar_id":           {"dtype": int,   "default": 1,    "is_list": False, "category": "optional"},
    "version_name":       {"dtype": str,   "default": "",   "is_list": False, "category": "optional"},
    "version_id":         {"dtype": int,   "default": None, "is_list": False, "category": "optional"},
    "configuration_name": {"dtype": str,   "default": "",   "is_list": False, "category": "optional"},
    "configuration_id":   {"dtype": int,   "default": None, "is_list": False, "category": "optional"},
    "station_altitude":   {"dtype": float, "default": None, "is_list": False, "category": "optional", "min": 0,    "max": 5000},
    "station_latitude":   {"dtype": float, "default": None, "is_list": False, "category": "optional", "min": -90,  "max": 90},
    "station_longitude":          {"dtype": float, "default": None, "is_list": False, "category": "optional", "min": -180, "max": 180},
    "zenith_angle":       {"dtype": float, "default": None, "is_list": False, "category": "optional", "min": 0,    "max": 85},
    "azimuth_angle":      {"dtype": float, "default": None, "is_list": False, "category": "optional", "min": 0,    "max": 360},

    # -------------------- [Channels] --------------------
    # Mandatory per comments
    "recorder_channel_id": {"dtype": str,   "default": None, "is_list": True,  "category": "mandatory"},
    "scc_channel_id":      {"dtype": int,   "default": None, "is_list": True,  "category": "optional"},
    "telescope_type":      {"dtype": str,   "default": None, "is_list": True,  "category": "mandatory", "allowed": ["n", "f", "x", "m", "g", "y", "l", "h", "z"]},
    "channel_type":        {"dtype": str,   "default": None, "is_list": True,  "category": "mandatory", "allowed": ["p", "c", "t", "v", "r", "a", "f"]},
    "channel_subtype":     {"dtype": str,   "default": None, "is_list": True,  "category": "mandatory", "allowed": ["r", "t", "n", "o", "w", "c", "h", "l", "a", "m", "b", "s", "x"]},

    # Partly optional (RECOMMENDED)
    "zero_bin":            {"dtype": int,   "default": 0,    "is_list": True,  "category": "recommended"},
    "dead_time":           {"dtype": float, "default": None, "is_list": True,  "category": "recommended", "min": 1.,   "max": 5.},
    "background_low_bin":  {"dtype": int,   "default": None, "is_list": True,  "category": "recommended", "min": 0,    "max": 32768},
    "background_high_bin": {"dtype": int,   "default": None, "is_list": True,  "category": "recommended", "min": 0,    "max": 32768},
    "channel_bandwidth":   {"dtype": float, "default": 1.0,  "is_list": True,  "category": "recommended", "min": 0.05, "max": 150.},
    "G":                   {"dtype": float, "default": None, "is_list": True,  "category": "recommended", "min": -2,   "max": 2.},
    "H":                   {"dtype": float, "default": None, "is_list": True,  "category": "recommended", "min": -2,   "max": 2.},

    # Optional overrides
    # "laser_id":                {"dtype": int,   "default": None, "is_list": True,  "category": "optional", "min": 1, "max": 2},
    "acquisition_mode":            {"dtype": str,   "default": None, "is_list": True, "category": "optional", "allowed": ["a", "p"]},
    "emitted_wavelength":          {"dtype": float, "default": None, "is_list": True, "category": "optional", "min": 200., "max": 3000.},
    "detected_wavelength":         {"dtype": float, "default": None, "is_list": True, "category": "optional", "min": 200., "max": 3000.},
    "bins":                        {"dtype": int,   "default": None, "is_list": True, "category": "optional", "min": 512,  "max": 32768},
    "data_acquisition_range":      {"dtype": float, "default": None, "is_list": True, "category": "optional", "allowed": [20., 100., 500.]},
    "range_resolution":            {"dtype": float, "default": None, "is_list": True, "category": "optional", "min": 1.,  "max": 30.},
    "laser_repetition_rate":       {"dtype": float, "default": None, "is_list": True, "category": "optional", "min": 10., "max": 200.},
    "analog_noise_per_bin":        {"dtype": float, "default": 0.22, "is_list": True, "category": "optional", "min": 0.},
    "analog_noise_scaling_factor": {"dtype": float, "default": 0.7,  "is_list": True, "category": "optional", "min": 0.},
    # "analog_to_digital_resolution": {"dtype": int, "default": None, "is_list": True, "category": "optional", "allowed": [12, 14, 16]},

    # -------------------- [polarization_calibration] --------------------
    "ch_r":                      {"dtype": str,   "default": None, "is_list": True, "category": "recommended"},
    "ch_t":                      {"dtype": str,   "default": None, "is_list": True, "category": "recommended"},
    "K":                         {"dtype": float, "default": 1.0,  "is_list": True, "category": "recommended", "min": 0.5,  "max": 2.},
    "R_to_T_transmission_ratio": {"dtype": float, "default": 1.0,  "is_list": True, "category": "recommended", "min": 1E-3, "max": 1E3},
}

# Groups to simplify length checks/expansion
SYSTEM_KEYS = {
    "station_id",
    "lidar_name",
    "station_name",
    "lidar_id",
    "version_name",
    "version_id",
    "configuration_name",
    "configuration_id",
    "station_altitude",
    "station_latitude",
    "station_longitude",
    "zenith_angle",
    "azimuth_angle"
    }

CHANNEL_KEYS = {
    "recorder_channel_id",
    "scc_channel_id",
    "telescope_type",
    "channel_type",
    "channel_subtype",
    # "laser_id",
    "zero_bin",
    "dead_time",
    "background_low_bin",
    "background_high_bin",
    "G",
    "H",
    "acquisition_mode",
    "detected_wavelength",
    "emitted_wavelength",
    "channel_bandwidth",
    "bins",
    "data_acquisition_range",
    # "analog_to_digital_resolution",
    "range_resolution",
    "laser_repetition_rate"
}
POL_CAL_KEYS = {"ch_r",
                "ch_t",
                "K",
                "R_to_T_transmission_ratio"}

recognized_sections = ["System", "Channels", "polarization_calibration"]

blank_tokens = ["_"]
# -------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------

def _is_empty_scalar(v: Any) -> bool:
    return v is None or (isinstance(v, str) and v.strip() == "")

def _convert_scalar(raw: Optional[str], expected: type, name: str) -> Optional[Union[str, int, float]]:
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
        if s not in ['True', 'False']:
            raise ConfigError(f"{name}: expected bool, please use either True or False")
        if s == 'True':
            s = True
        else:
            s = False
        return s
    raise ConfigError(f"{name}: unsupported dtype {expected}")
    
def normalize_entry(parser_args, key):
    """
    Normalize a single entry in parser_args according to SCHEMA[key]["dtype"].

    Returns the converted value (does NOT modify parser_args).
    """
    if key not in SCHEMA:
        raise KeyError(f"{key} not found in SCHEMA")

    meta = SCHEMA[key]
    dtype = meta["dtype"]
    is_list = meta.get("is_list", False)

    value = parser_args.get(key)

    if value is None:
        return None

    if is_list:
        normalized = []
        for v in value:
            if v is None:
                normalized.append(None)
            elif isinstance(v, str) and v.strip() == "":
                normalized.append(None)
            elif v in blank_tokens:
                normalized.append(None)
            else:
                normalized.append(_convert_scalar(str(v), dtype, key))
        return normalized

    else:
        if value is None:
            return None
        elif isinstance(value, str) and value.strip() == "":
            return None
        elif value in blank_tokens:
            return None
        else:
            return _convert_scalar(str(value), dtype, key)

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
        if items[i] in blank_tokens:
            items[i] = None
        else:
            items[i] = _convert_scalar(items[i], str, name)
    return out

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
            raise ConfigError(f"{name}: value {value} is not a list despite being defined as a list in the SCHEMA")
        else:
            if len(value) > 0:
                for v in value or []:
                    check_one(v)
    else:
        if isinstance(value, list):
            raise ConfigError(f"{name}: value {value} is a list despite being defined as a scalar parameter in the SCHEMA")
        else:
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
    if recommended_missing == True:
        print(msg)
    print(f"  --{name}")

# -------------------------------------------------------------------
# Expansion & computed defaults
# -------------------------------------------------------------------

def _expand_channel_lists_with_defaults(parser_args: Dict[str, Any]) -> None:
    """Ensure every [Channels] list has length N (recorder_channel_id), filling from schema default scalar."""
    recorder_channel_id = parser_args.get("recorder_channel_id") or []
    n_channel = len(recorder_channel_id)
    
    for key in CHANNEL_KEYS:
        arg = parser_args.get(key)
        meta = SCHEMA[key]
        if meta.get("is_list"):
            if arg is None or len(arg) == 0:
                # expand scalar default (which may be None)
                parser_args[key] = [meta.get("default")] * n_channel
            else:
                if len(arg) != n_channel:
                    raise ConfigError(f"{key} length {len(arg)} must equal recorder_channel_id length {n_channel}")

def _expand_pol_cal_lists_with_defaults(parser_args: Dict[str, Any]) -> None:
    """Ensure every [polarisation_calibration] list has length N (ch_r), filling from schema default scalar."""
    ch_r = parser_args.get("ch_r")
    n_pairs = len(ch_r)
    
    if n_pairs > 0:
        # Expand K and R_to_T_transmission_ratio if omitted
        for key in POL_CAL_KEYS:
            if key in ("ch_r", "ch_t"):
                continue
            arg = parser_args.get(key)
            if arg is None or len(arg) == 0:
                parser_args[key] = [SCHEMA[key]["default"]] * n_pairs
            elif len(arg) != n_pairs:
                raise ConfigError(f"{key} length {len(arg)} must equal {n_pairs} length {n_pairs}")
    
    return parser_args

def _compute_dead_time_if_missing(parser_args: Dict[str, Any]) -> None:

    for i, ch in enumerate(parser_args["recorder_channel_id"]):
        dt = parser_args["dead_time"][i]
        aq_mode = parser_args["acquisition_mode"][i]
        if aq_mode == "p" and dt in [None, np.inf, np.nan]:
            parser_args["dead_time"][i] = "3.7"
            
    return parser_args

def _compute_emitted_wavelength_if_missing(parser_args: Dict[str, Any]) -> None:

    for i, ch in enumerate(parser_args["recorder_channel_id"]):
        dt_wv = parser_args["detected_wavelength"][i]
        em_wv = parser_args["emitted_wavelength"][i]
        if dt_wv in [None, np.inf, np.nan]:
            raise Exception(f"detected_wavelength is empty for channel {ch}: {dt_wv}")
        elif em_wv in [None, np.inf, np.nan]:
            if float(dt_wv) <= 520.:
                parser_args["emitted_wavelength"][i] = 354.71
            elif float(dt_wv) > 520. and float(dt_wv) <= 1040.:
                parser_args["emitted_wavelength"][i] = 532.07
            else:
                parser_args["emitted_wavelength"][i] = 1064.14
        if dt_wv % 1 == 0:
            CustomWarning(f"detected wavelength for channel {ch} is too round: {dt_wv}\nConsider providing it with accuracy of 2 decimal points (central wavelength of the interference filter)\n")
              
    return parser_args

def _compute_background_bins_if_missing(parser_args: Dict[str, Any]) -> None:
    
    for i, ch in enumerate(parser_args["recorder_channel_id"]):
        bhi = parser_args["background_high_bin"][i]
        blo = parser_args["background_low_bin"][i]
        zbi = parser_args["zero_bin"][i]
        bns = parser_args["bins"][i]
                
        if zbi in [None, np.inf, np.nan]:
            raise Exception(f"zero_bin is empty for channel {ch}: {zbi}")
        
        if bns in [None, np.inf, np.nan]:
            raise Exception(f"bins is empty for channel {ch}: {bns}")
            
        if (blo == None and bhi != None) or (bhi == None and blo != None):
            raise ConfigError("background_low_bin {blo} and background_high_bin {bhi} must either be both provided or be both left empty")
    
        if blo != None and bhi != None:              
            if blo > bhi:
                raise ConfigError(f"background_low_bin values {blo} must be smaller than background_high_bin values {bhi}")
        
        if blo == None and bhi == None:
            if int(zbi) < -240:
                parser_args["background_low_bin"][i] = 50
                parser_args["background_high_bin"][i] = abs(int(zbi)) - 40
            else:
                parser_args["background_low_bin"][i] = int(0.9*int(bns))
                parser_args["background_high_bin"][i] = int(bns) - 50
        
    return parser_args

def exceeding_limits_error(var1, var2, label1, label2):
    
    errors = []
    for i, (a, b) in enumerate(zip(var1, var2)):
        if b is not None and a > b:
            errors.append(f"index {i}: {a} < {b}")
    
    if errors:
        raise ConfigError(
            "{label1} values must be smaller than the corresponding {label2}:\n" + \
                "{label1} = {var1}\n" + \
                    "{label2} = {var2}"
        )
    
    
def _check_background_limits(parser_args: Dict[str, Any]) -> None:
    
    low_bin = parser_args['background_low_bin']
    high_bin = parser_args['background_high_bin']
    
    bins = parser_args['bins']
    
    exceeding_limits_error(
        var1 = low_bin, 
        var2 = bins, 
        label1 = 'background_low_bin',
        label2 = 'number of bins',
        )

    exceeding_limits_error(
        var1 = high_bin, 
        var2 = bins, 
        label1 = 'background_high_bin',
        label2 = 'number of bins',
        )    
   
    
def _compute_GH_if_missing(parser_args: Dict[str, Any]) -> None:

    for i, ch in enumerate(parser_args["recorder_channel_id"]):
        G = parser_args["G"][i]
        H = parser_args["H"][i]
        ch_t = parser_args["channel_type"][i]
        
        if G in [None, np.inf, np.nan]:
            parser_args["G"][i] = 1.
        if H in [None, np.inf, np.nan]:   
            if ch_t == 'c':
                parser_args["H"][i] = -1.
            elif ch_t == 'p':
                parser_args["H"][i] = 1.
            else:
                parser_args["H"][i] = 0.
                
    return parser_args

# -------------------------------------------------------------------
# Presence checks
# -------------------------------------------------------------------

def _enforce_mandatory_and_recommended(parser_args: Dict[str, Any]) -> None:
    recommended_first_time = True
    for name, meta in SCHEMA.items():
        cat = meta.get("category", "optional")
        val = parser_args.get(name)
        if meta.get("is_list", False):
            # treat list as empty if all elements are None/""
            is_empty = (val is None) or (len(val) == 0) or all(_is_empty_scalar(v) for v in val)
        else:
            is_empty = _is_empty_scalar(val)
        if cat == "mandatory" and is_empty:
            raise ConfigError(f"{name} is mandatory and was not provided.")
        if cat == "recommended" and is_empty:
            _warn_recommended(name, "--Warning: Recomended configuration parameters not provided:", recommended_first_time)
            recommended_first_time = False

def _recorder_channel_id_check(parser_args: Dict[str, Any]) -> None:
    
    recorder_channel_id = parser_args.get("recorder_channel_id")
    laser_id = parser_args.get("laser_id")
    
    if len(recorder_channel_id) == 0:
        raise ConfigError("recorder_channel_id must list at least one channel.")
    n_channel = len(recorder_channel_id)
    # Validate that user-provided channel lists (if present) already have proper lengths
    for key in CHANNEL_KEYS:
        arg = parser_args.get(key)
        # If user provided a non-empty list, enforce same length now
        if len(arg) != 0 and len(arg) != n_channel:
            raise ConfigError(f"{key} length {len(arg)} must equal recorder_channel_id length {n_channel}.")
            
    if len(recorder_channel_id) != len(set(recorder_channel_id)):
        # assert_pairs_unique(recorder_channel_id, laser_id)
        raise ConfigError(f"{key} contains duplicates")

# def assert_pairs_unique(recorder_channel_id, laser_id):
#     if len(recorder_channel_id) != len(laser_id):
#         raise ConfigError("Duplicate recorder_channel_id values found {recorder_channel_id}. laser_id is mandatory in this case")
#     seen = set()
#     dups = []
#     for i, pair in enumerate(zip(recorder_channel_id, laser_id)):
#         if pair in seen:
#             dups.append((i, pair))
#         else:
#             seen.add(pair)
#     if dups:
#         raise ConfigError(f"Duplicate recorder_channel_id values found {recorder_channel_id}. laser_id was provided {laser_id} but it is the same for at least one of the duplicate recorder_channel_id values")
        
def _polarisation_calibration_check(parser_args: Dict[str, Any]) -> None:
    
    ch_r = parser_args.get("ch_r")
    ch_t = parser_args.get("ch_t")
    n_pairs = len(ch_r)
    # Validate that user-provided channel lists (if present) already have proper lengths
    for key in POL_CAL_KEYS:
        arg = parser_args.get(key)
        # If user provided a non-empty list, enforce same length now
        if len(arg) != n_pairs:
            raise ConfigError(f"{key} length {len(arg)} must equal ch_r length {n_pairs}.")

    if ch_r != []:
        if parser_args["recorder_channel_id"] == []:
            raise ConfigError(f"ch_r was provided but the recorder_channel_id is empty. Parameter ch_r must be a subset of recorder_channel_id")

        for r in ch_r:
            if r not in parser_args["recorder_channel_id"]:
                raise ConfigError(f"ch_r: {r} - Value not in recorder_channel_id: {parser_args['recorder_channel_id']}")

    if ch_t != []:
        if parser_args["recorder_channel_id"] == []:
            raise ConfigError(f"ch_t was provided but the recorder_channel_id is empty. Parameter ch_t must be a subset of recorder_channel_id")

        for t in ch_t:
            if r not in parser_args["recorder_channel_id"]:
                raise ConfigError(f"ch_t: {t} - Value not in recorder_channel_id: {parser_args['recorder_channel_id']}")


def read_ini_file(filepath: str)  -> Dict[str, Any]:
    
    """Read, convert, expand from scalar defaults, compute simple defaults, validate, and return dict."""
    config = configparser.ConfigParser(allow_no_value=True, strict=True)

    read_files = config.read(filepath, encoding="utf-8")

    if not read_files:
        raise ConfigError(f"INI file not found or unreadable: {filepath}\nMake sure the encoding is utf-8")

    for section in config.sections():
        if section not in recognized_sections:
            raise ConfigError(f"{section} section is not recognised. Please revise the settings file. Recognized sections: {[sec for sec in recognized_sections]}")
    
    for section in recognized_sections:
        if section not in config.sections():
            CustomWarning(f"{section} section not found in the configuration file")

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
                    parser_args[key] = _convert_scalar(raw, str, key)
            
                  
        if found == False:
            if meta["is_list"]:
                parser_args[key] = []
            else:
                parser_args[key] = meta["default"]
    
    return parser_args

def _normalize_entries(parser_args: Dict[str, Any]) -> None:
    
    for key in parser_args.keys():
        parser_args[key] = normalize_entry(parser_args, key)
    
# -------------------------------------------------------------------
# Public API
# -------------------------------------------------------------------

def parse_atlas_config_file(filepath: str, debug: bool = False) -> Dict[str, Any]:

    print_header(f"Parsing Configuration file:\n{filepath}")

    # Parse & convert
    parser_args = read_ini_file(filepath)    

    # Normalize entries
    _normalize_entries(parser_args) 
      
    # Enforce mandatory/recommended
    _enforce_mandatory_and_recommended(parser_args)

    # Basic presence and uniquiness checks before expansion
    # recorder_channel_id must exist (mandatory) and not be empty
    _recorder_channel_id_check(parser_args)
    
    # Expand channel lists using scalar defaults (including None) when omitted
    _expand_channel_lists_with_defaults(parser_args)

    # Expand polarisation calibration lists using scalar defaults (including None) when omitted
    _expand_pol_cal_lists_with_defaults(parser_args)
    
    # Length check polarisation calibratio
    # all variables must have the same length (either all 0 if not provided or all the same legnth)
    _polarisation_calibration_check(parser_args)
    
    # Range & allowed checks (only on non-empty values)
    for name, meta in SCHEMA.items():
        _check_limits(name, parser_args.get(name), meta)
        _check_allowed(name, parser_args.get(name), meta)

    # Special checks
    # Check if the background limits exceed the number of bins (if provided)
    _check_background_limits(parser_args)
    
    # Sort keys alphabetically
    parser_args = dict(sorted(parser_args.items()))

    if debug:
        pprint(parser_args)
        
    return parser_args
