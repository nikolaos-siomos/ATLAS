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
import os
from helper_functions.printouts import print_header
from pprint import pprint
from utils.error_classes import ConfigError

Number = Union[int, float]

# -------------------------------------------------------------------
# SCHEMA
#   dtype:        expected Python type (str, int, float)
#   default:      scalar default (always scalar; expanded to lists as needed)
#   is_list:      True if value is a list in INI (comma or semicolon-separated)
#   allowed:      optional list of allowed values (checked only on non-empty values)
#   min/max:      optional numeric bounds
#   size:         The list ust have a specific size if not empty
# -------------------------------------------------------------------

SCHEMA_QCK: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "t_lims":      {"dtype": str,   "default": [],        "is_list": True},
    "t_tick":      {"dtype": int,   "default": None,      "is_list": False, "min": 5},
    "y_lims":      {"dtype": float, "default": [0., 14.], "is_list": True,  "min": 0., "max": 100., "size": 2},
    "y_tick":      {"dtype": float, "default": 1.,        "is_list": False, "min": 0., "max": 10},   
    "z_lims":      {"dtype": float, "default": [],        "is_list": True, "size": 2}, 
    "z_max_zone":  {"dtype": float, "default": [0.1, 2.], "is_list": True, "min": 0., "max": 30., "size": 2},

    "smooth":                           {"dtype": bool,  "default": False,       "is_list": False},
    "smoothing_range":                  {"dtype": float, "default": [0.05, 15.], "is_list": True,  "min": 0.,   "max": 100., "size": 2},
    "smoothing_constant_window":        {"dtype": float, "default": None,        "is_list": False, "min": 0.05, "max": 10.},
    "smoothing_progressive_window":     {"dtype": float, "default": [0.05, 0.5], "is_list": True,  "min": 0.05, "max": 10.,  "size": 2},
    "smoothing_exponential":            {"dtype": bool,  "default": False,       "is_list": False},

    "channels":                 {"dtype": str,"default": [], "is_list": True},
    "exclude_telescope_type":   {"dtype": str,"default": [], "is_list": True,"allowed": ["n", "f", "x", "m", "g", "y", "l", "h", "z"]},
    "exclude_channel_type":     {"dtype": str,"default": [], "is_list": True,"allowed": ["p", "c", "t", "v", "r", "a", "f"]},
    "exclude_acquisition_mode": {"dtype": str,"default": [], "is_list": True,"allowed": ["a", "p", "g"]},
    "exclude_channel_subtype":  {"dtype": str,"default": [], "is_list": True,"allowed": ["r", "t", "n", "o", "w", "c", "h", "l", "a", "m", "b", "s", "x"]},
}

SCHEMA_RAY: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "x_lims":      {"dtype": float,   "default": [0.,31.], "is_list": True,  "min": 0.,  "max": 100.},
    "x_tick":      {"dtype": float,   "default": 2.,       "is_list": False, "min": 0.1, "max": 10.},   
    "y_lims":      {"dtype": float,   "default": [],       "is_list": True},

    "use_lin_scale":               {"dtype": bool,    "default": True,      "is_list": False},
    "cross_check_lower_limit":     {"dtype": float,   "default": None,      "is_list": False},
    "normalization_region":        {"dtype": float,   "default": [6.,8.],   "is_list": True,  "min": 0.,   "max": 50.},
    
    "molecular_mask_region":       {"dtype": float,   "default": [2., 30.], "is_list": True,  "min": 0.,    "max": 50.},
    "molecular_mask_window":       {"dtype": float,   "default": 0.5,       "is_list": False, "min": 0.05,   "max": 5.},
    "molecular_mask_window_step":  {"dtype": float,   "default": 0.2,       "is_list": False, "min": 0.05,  "max": 5.},
    "rsem_threshold":              {"dtype": float,   "default": 0.2,       "is_list": False, "min": 0.01,  "max": 1.},
    "first_derivative_threshold":  {"dtype": float,   "default": 2.,        "is_list": False, "min": 0.5,   "max": 5.},
    "second_derivative_threshold": {"dtype": float,   "default": 2.,        "is_list": False, "min": 0.5,   "max": 5.},
    "shapiro_wilk_threshold":      {"dtype": float,   "default": 0.05,      "is_list": False, "min": 0.,    "max": 1.},
    "cross_criterium_threshold":   {"dtype": float,   "default": 1.,        "is_list": False, "min": 0.5,   "max": 5.},
    
    "smooth":                           {"dtype": bool,  "default": True,        "is_list": False},
    "smoothing_range":                  {"dtype": float, "default": [0.05, 31.], "is_list": True,  "min": 0.,   "max": 100., "size": 2},
    "smoothing_constant_window":        {"dtype": float, "default": 0.5,         "is_list": False, "min": 0.05, "max": 10.},
    "smoothing_progressive_window":     {"dtype": float, "default": [],          "is_list": True,  "min": 0.05, "max": 10.,  "size": 2},
    "smoothing_exponential":            {"dtype": bool,  "default": False,       "is_list": False},

    "channels":                 {"dtype": str,"default": [], "is_list": True},
    "exclude_telescope_type":   {"dtype": str,"default": [], "is_list": True,"allowed": ["n", "f", "x", "m", "g", "y", "l", "h", "z"]},
    "exclude_channel_type":     {"dtype": str,"default": ["a","f"], "is_list": True,"allowed": ["p", "c", "t", "v", "r", "a", "f"]},
    "exclude_acquisition_mode": {"dtype": str,"default": [], "is_list": True,"allowed": ["a", "p", "g"]},
    "exclude_channel_subtype":  {"dtype": str,"default": ["w", "c"], "is_list": True,"allowed": ["r", "t", "n", "o", "w", "c", "h", "l", "a", "m", "b", "s", "x"]},
}

SCHEMA_TLC: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "x_lims":      {"dtype": float, "default": [],   "is_list": True},
    "x_tick":      {"dtype": float, "default": None, "is_list": False},   

    "use_non_rangecor":       {"dtype": bool,    "default": False,      "is_list": False},
    "use_last":               {"dtype": bool,    "default": False,      "is_list": False},
    "normalization_region":   {"dtype": float,   "default": [1.8, 2.2], "is_list": True},
    
    "smooth":                           {"dtype": bool,  "default": False,       "is_list": False},
    "smoothing_constant_window":        {"dtype": float, "default": 0.1,         "is_list": False, "min": 0.05, "max": 10.},
    "smoothing_exponential":            {"dtype": bool,  "default": False,       "is_list": False},

    "channels":                 {"dtype": str,"default": [], "is_list": True},
    "exclude_telescope_type":   {"dtype": str,"default": [], "is_list": True,"allowed": ["n", "f", "x", "m", "g", "y", "l", "h", "z"]},
    "exclude_channel_type":     {"dtype": str,"default": [], "is_list": True,"allowed": ["p", "c", "t", "v", "r", "a", "f"]},
    "exclude_acquisition_mode": {"dtype": str,"default": [], "is_list": True,"allowed": ["a", "p", "g"]},
    "exclude_channel_subtype":  {"dtype": str,"default": [], "is_list": True,"allowed": ["r", "t", "n", "o", "w", "c", "h", "l", "a", "m", "b", "s", "x"]},
}


SCHEMA_TLC_RIN: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "x_lims":      {"dtype": float, "default": [],   "is_list": True},
    "x_tick":      {"dtype": float, "default": None, "is_list": False},   

    "use_non_rangecor":       {"dtype": bool,    "default": False,      "is_list": False},
    "use_last":               {"dtype": bool,    "default": False,      "is_list": False},
    "normalization_region":   {"dtype": float,   "default": [1.8, 2.2], "is_list": True},
    
    "smooth":                           {"dtype": bool,  "default": False,       "is_list": False},
    "smoothing_constant_window":        {"dtype": float, "default": 0.1,         "is_list": False, "min": 0.05, "max": 10.},
    "smoothing_exponential":            {"dtype": bool,  "default": False,       "is_list": False},

    "channels":                 {"dtype": str,"default": [], "is_list": True},
    "exclude_telescope_type":   {"dtype": str,"default": [], "is_list": True,"allowed": ["n", "f", "x", "m", "g", "y", "l", "h", "z"]},
    "exclude_channel_type":     {"dtype": str,"default": [], "is_list": True,"allowed": ["p", "c", "t", "v", "r", "a", "f"]},
    "exclude_acquisition_mode": {"dtype": str,"default": [], "is_list": True,"allowed": ["a", "p", "g"]},
    "exclude_channel_subtype":  {"dtype": str,"default": [], "is_list": True,"allowed": ["r", "t", "n", "o", "w", "c", "h", "l", "a", "m", "b", "s", "x"]},
}

SCHEMA_PCB: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "x_lims_calibration":   {"dtype": float,   "default": [0., 15.], "is_list": True,  "min": 0.,  "max": 50.},
    "x_lims_rayleigh":      {"dtype": float,   "default": [0., 31.], "is_list": True,  "min": 0.,  "max": 50.},
    "x_tick_calibration":   {"dtype": float,   "default": 0.5,       "is_list": False, "min": 0.1, "max": 10.},   
    "x_tick_rayleigh":      {"dtype": float,   "default": 2.,        "is_list": False, "min": 0.1, "max": 10.},   
    "y_lims_calibration":   {"dtype": float,   "default": [],      "is_list": True},
    "y_lims_rayleigh":      {"dtype": float,   "default": [],      "is_list": True},

    "calibration_region": {"dtype": float,   "default": [2., 4.], "is_list": True,  "min": 0.,  "max": 50., "size": 2},
    "rayleigh_region":    {"dtype": float,   "default": [6., 8.], "is_list": True, "size": 2},
    
    "smooth":                           {"dtype": bool,  "default": True,        "is_list": False},
    "smoothing_range":                  {"dtype": float, "default": [0.05, 31.], "is_list": True,  "min": 0.,   "max": 100., "size": 2},
    "smoothing_constant_window":        {"dtype": float, "default": 0.5,         "is_list": False, "min": 0.05, "max": 10.},
    "smoothing_progressive_window":     {"dtype": float, "default": [],          "is_list": True,  "min": 0.05, "max": 10.,  "size": 2},
    "smoothing_exponential":            {"dtype": bool,  "default": False,       "is_list": False},

}

SCHEMA = {'qck': SCHEMA_QCK,
          'ray': SCHEMA_RAY,
          'tlc': SCHEMA_TLC,
          'tlc_rin': SCHEMA_TLC_RIN,
          'pcb': SCHEMA_PCB}

recognized_sections = {'qck': "quicklooks",
                       'ray': "rayleigh_fit",
                       'tlc': "quadrant_telecover",
                       'tlc_rin': "ring_telecover",
                       'pcb': "polarization_calibration"}

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
        if len(value) != meta["size"]:
            raise ConfigError(f"{name} must have size: {meta['size']}\nThe following values were provided: {value}")      
        
def _check_limits(name: str, value: Any, meta: Dict[str, Any]) -> None:
    def check_one(v: Any):

        if v is None or v is [] or isinstance(v, str):
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

# -------------------------------------------------------------------
# Expansion & computed defaults
# -------------------------------------------------------------------

def _fill_with_defaults(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    """Fill empty entries with default values"""

    for key in SCHEMA.keys():
        if key not in parser_args.keys():
            parser_args[key] = {}
        for sec_key, meta in SCHEMA[key].items():
            arg = parser_args[key].get(sec_key)
            default = meta["default"]
            if arg == None or arg == []:
                parser_args[key][sec_key] = default
    
    return(parser_args)

def _special_checks(parser_args: Dict[str, Any]) -> None:
    
    for key in SCHEMA.keys():
        for sec_key, meta in SCHEMA[key].items():
            if "smoothing_constant_window" in parser_args[key].keys() and \
                "smoothing_progressive_window" in parser_args[key].keys():
                if parser_args[key]["smoothing_constant_window"] != None and \
                    parser_args[key]["smoothing_progressive_window"] != []:
                    raise ConfigError(f"Section {key}: smoothing_constant_window and smoothing_progressive_window cannot be provided at the same time. Please provide only one of them (the other can be left empty)")


# -------------------------------------------------------------------
# Presence checks
# -------------------------------------------------------------------

def read_ini_file(filepath: str)  -> Dict[str, Any]:
    
    if filepath != None:

        """Read, convert, expand from scalar defaults, compute simple defaults, validate, and return dict."""
        config = configparser.ConfigParser(allow_no_value=True, strict=True)
    
        read_files = config.read(filepath, encoding="utf-8")
    
        if not read_files:
            print(f"--Warning: INI file not found or unreadable: {filepath}\n Default values will be used")
        else:
            for section in config.sections():
                if section not in [recognized_sections[k] for k in recognized_sections.keys()]:
                    raise ConfigError(f"{section} section is not recognised. Please revise the settings file. Recognized sections: {[recognized_sections[k] for k in recognized_sections.keys()]}")
            
            for section in recognized_sections.values():
                if section not in config.sections():
                    print(f"--Warning: {section} section not found in the settings file")
    
        parser_args: Dict[str, Any] = {}
        
        for key in SCHEMA.keys():
            parser_args[key] = {}
            section = recognized_sections[key]
            for sec_key, meta in SCHEMA[key].items():
                if section in config.sections():
                    if sec_key in config[section]:
                        raw = config[section][sec_key]
                        if meta["is_list"]:
                            parser_args[key][sec_key] = _convert_list(raw, meta, sec_key)
                        else:
                            parser_args[key][sec_key] = _convert_scalar(raw, meta["dtype"], sec_key)
                    else:
                        if meta["is_list"]:
                            parser_args[key][sec_key] = []
                        else:
                            parser_args[key][sec_key] = meta["default"]
                else:
                    if meta["is_list"]:
                        parser_args[key][sec_key] = []
                    else:
                        parser_args[key][sec_key] = meta["default"] 
    else:
        parser_args: Dict[str, Any] = {}
        print(f"--Warning: No settings file was provided. Default values will be assumed")

    return parser_args

# -------------------------------------------------------------------
# Public API
# -------------------------------------------------------------------

def parse_atlas_settings_file(filepath: str, debug: bool = False) -> Dict[str, Any]:

    print_header(f"Parsing settings file:\n{filepath}")
    
    # 1) Parse & convert
    parser_args = read_ini_file(filepath)    
    
    # 2) Fill with default values if empty
    parser_args = _fill_with_defaults(parser_args)
    
    # 3) Range & allowed checks (only on non-empty values)
    for key in SCHEMA.keys():
        for sec_key, meta in SCHEMA[key].items():
            _check_size(sec_key, parser_args[key].get(sec_key), meta)
            _check_limits(sec_key, parser_args[key].get(sec_key), meta)
            _check_allowed(sec_key, parser_args[key].get(sec_key), meta)
    

    # 4) Check the smoothing window options
    _special_checks(parser_args)
    
    # 5) Sort keys alphabetically
    for key in SCHEMA.keys():
        parser_args[key] = dict(sorted(parser_args[key].items()))
        
    if debug:
        pprint(parser_args)
        
    return parser_args
