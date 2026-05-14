#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:16:34 2025

@author: nikos
"""

from __future__ import annotations
import configparser
from typing import Any, Dict, List, Optional, Union, Mapping, Iterable, Sequence
from pathlib import Path
import numpy as np
import os, sys
from utils.caller_utils import autodetect_paths
from helper_functions.printouts import print_header, endpoint
from pprint import pprint
from utils.error_classes import ConfigError, CustomWarning
from datetime import datetime
import re

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
    'ray', 
    'pcb', 
    'tlc', 
    'tlc_rin', 
    'drk', 
    'trg', 
    'dtm', 
    'nsf', 
    'ray_pcb', 
    'pcb_aux', 
    'cam'
    ]

qck_tests = [
    'ray', 
    'pcb', 
    'tlc', 
    'tlc_rin', 
    'drk', 
    'ray_pcb', 
    'pcb_aux', 
    ]
recognised_slicers = [
    'ray', 
    'pcb_p45', 
    'pcb_m45', 
    'tlc_north', 
    'tlc_east', 
    'tlc_south', 
    'tlc_west', 
    'tlc_rin_inner', 
    'tlc_rin_outer', 
    'trg', 
    'dtm', 
    'nsf', 
    'ray_pcb', 
    'pcb_aux_p45', 
    'pcb_aux_m45', 
    'drk', 
    'drk_ray', 
    'drk_pcb', 
    'drk_tlc', 
    'drk_tlc_rin', 
    'drk_trg', 
    'drk_dtm', 
    'drk_nsf', 
    'drk_ray_pcb', 
    'drk_pcb_aux'
    ]
height_units = ['m_asl', 'm_agl', 'km_asl', 'km_agl']
pressure_units = ['Pa', 'hPa', 'atm']
temperature_units = ['K', 'C', 'Cx10']
humidity_units = ['percent','fraction']

SCHEMA: Dict[str, Dict[str, Any]] = {
    # -------------------- [System] --------------------
    "main_data_folder":          {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "dir"},
    "scc_station_id":            {"dtype": str, "default": None, "is_list": False, "category": "optional"},
    "scc_configuration_id":      {"dtype": str, "default": None, "is_list": False, "category": "optional"},
    "data_identifier":           {"dtype": str, "default": None, "is_list": False, "category": "optional"},
    "export_hoi_cfg":            {"dtype": str, "default": "0",  "is_list": False, "category": "optional", "allowed": ["0", "1", "2"]},
   
    "parent_folder":             {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "dir"},
    "atlas_configuration_file":  {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "file"},
    "atlas_settings_file":       {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "file"},
    "radiosonde_folder":         {"dtype": str, "default": None, "is_list": False, "category": "optional", "check_path": "dir"},

    "quick_run":           {"dtype": bool, "default": False,    "is_list": False, "category": "optional"},
    "process":             {"dtype": str,  "default": qa_tests, "is_list": True,  "category": "optional", "allowed": qa_tests + ['off']},
    "process_qck":         {"dtype": str,  "default": qck_tests,"is_list": True,  "category": "optional", "allowed": qa_tests + ['off']},
    "use_range":           {"dtype": bool, "default": True,     "is_list": False, "category": "optional"},
    "dpi":                 {"dtype": int,  "default": 300,      "is_list": False, "category": "optional"},
    "output_folder":       {"dtype": str,  "default": None,     "is_list": False, "category": "optional"},
    "expert_analyst":      {"dtype": str,  "default": None,     "is_list": False, "category": "optional"},
    "debug_signals":       {"dtype": bool, "default": False,    "is_list": False, "category": "optional"},
    "export_netcdf":       {"dtype": bool, "default": True,     "is_list": False, "category": "optional"},
    "export_all":          {"dtype": bool, "default": False,    "is_list": False, "category": "optional"},
    
    "channels":                 {"dtype": str,"default": [], "is_list": True, "category": "optional"},
    "exclude_telescope_type":   {"dtype": str,"default": [], "is_list": True, "category": "optional","allowed": ["n", "f", "x", "m", "g", "y", "l", "h", "z"]},
    "exclude_channel_type":     {"dtype": str,"default": [], "is_list": True, "category": "optional","allowed": ["p", "c", "t", "v", "r", "a", "f"]},
    "exclude_acquisition_mode": {"dtype": str,"default": [], "is_list": True, "category": "optional","allowed": ["a", "p", "g"]},
    "exclude_channel_subtype":  {"dtype": str,"default": [], "is_list": True, "category": "optional","allowed": ["r", "t", "n", "o", "w", "c", "h", "l", "a", "m", "b", "s", "x"]},

    "max_height_agl":               {"dtype": float, "default": 32.,  "is_list": False, "category": "optional"},
    "low_shot_threshold":           {"dtype": float, "default": 0.9,  "is_list": False, "category": "optional", "min": 0., "max": 0.999},
    "trim_overflows":               {"dtype": int,   "default": 0,    "is_list": False, "category": "optional", "allowed": [0,1,2,3]},
    "ray_averaging_rate":           {"dtype": str,   "default": None, "is_list": False, "category": "optional"},
    "ray_averaging_threshold":      {"dtype": float, "default": 0.5,  "is_list": False, "category": "optional", "min": 0., "max": 1.},
    "ray_qck_averaging_rate":       {"dtype": str,   "default": None, "is_list": False, "category": "optional"},
    "ray_qck_averaging_threshold":  {"dtype": float, "default": 0.5,  "is_list": False, "category": "optional", "min": 0., "max": 1.},
    "max_adjacent_overflows":       {"dtype": int,   "default": 5,    "is_list": False, "category": "optional", "min": 0},
    "slice_measurement":            {"dtype": str,   "default": [],   "is_list": True,  "category": "optional"},
    "exclude_measurement":          {"dtype": str,   "default": [],   "is_list": True,  "category": "optional"},
    
    "ray":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "nrm":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "pcb":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "tlc":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "tlc_rin":      {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "drk":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "trg":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "dtm":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "nsf":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "ray_pcb":      {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "nrm_pcb":      {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "pcb_aux":      {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},
    "cam":          {"dtype": str,"default": None, "is_list": False, "category": "optional", "check_path": "relative"},


    "files_per_sector":         {"dtype": int,   "default": None,     "is_list": False, "category": "optional", "min": 0},
    "files_per_ring":           {"dtype": int,   "default": None,     "is_list": False, "category": "optional", "min": 0},
    "files_per_set":            {"dtype": int,   "default": None,     "is_list": False, "category": "optional", "min": 0},
    "rsonde_skip_header":       {"dtype": int,   "default": 1,        "is_list": False, "category": "optional", "min": 0},
    "rsonde_skip_footer":       {"dtype": int,   "default": 0,        "is_list": False, "category": "optional", "min": 0},
    "rsonde_delimiter":         {"dtype": str,   "default": "S",      "is_list": False, "category": "optional", "allowed": ["S", "C"]},
    "rsonde_column_index":      {"dtype": int,   "default": [2,1,3,5], "is_list": True, "category": "optional", "min": 0, "size": [3, 4]},
    "rsonde_column_units":      {"dtype": str,   "default": ['m_asl', 'hPa', 'C', 'percent'], "is_list": True, "category": "optional", "size": [3, 4], "allowed": height_units + pressure_units + temperature_units + humidity_units},
    "rsonde_station_latitude":  {"dtype": float, "default": None,     "is_list": False, "category": "optional", "min": -90,  "max": 90},
    "rsonde_station_longitude": {"dtype": float, "default": None,     "is_list": False, "category": "optional", "min": -180, "max": 180},
    "rsonde_station_altitude":  {"dtype": float, "default": None,     "is_list": False, "category": "optional", "min": 0,  "max": 5000},
    "cloudnet_station_name":    {"dtype": str,   "default": None,     "is_list": False, "category": "optional"},
    "rsonde_station_name":      {"dtype": str,   "default": None,     "is_list": False, "category": "optional"},
    "rsonde_station_wmo_id":    {"dtype": int,   "default": None,     "is_list": False, "category": "optional", "min": 0},
}

explicit_path_keys = {"parent_folder", 
                      "atlas_configuration_file", 
                      "atlas_settings_file", 
                      "radiosonde_folder"}

tlc_subfolders = ["north", "east", "south", "west"]
tlc_rin_subfolders = ["inner", "outer"]
pcb_subfolders = ["p45", "m45", "+45","-45"]

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
    def check_one(v: Any):
        if v is None or isinstance(v, str):
            return
        if "min" in meta and v < meta["min"]:
            raise ConfigError(f"{name}: value {v} below minimum {meta['min']}")
        if "max" in meta and v > meta["max"]:
            raise ConfigError(f"{name}: value {v} above maximum {meta['max']}")
            
    if meta.get("is_list") and "size" in meta and len(value) > 0:
        if len(value) not in meta["size"]:
            if isinstance(meta["size"], list):
                raise ConfigError(f"{name} must have one of the following sizes: {meta['size']}\nThe following values were provided: {value}")
            else:
                raise ConfigError(f"{name} must have size: {meta['size']}\nThe following values were provided: {value}")      
        
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

def _fill_with_defaults(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    """Fill empty entries with default values"""

    for key, meta in SCHEMA.items():
        arg = parser_args.get(key)
        default = meta["default"]
        if arg == None or arg == []:
            parser_args[key] = default
    
    return(parser_args)
            

def _compute_emitted_wavelength_if_missing(parser_args: Dict[str, Any]) -> None:

    det = parser_args.get("detected_wavelength")
    em  = parser_args.get("emitted_wavelength")
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

    pattern = re.compile(r'^\d{1,2}(?:min|h)$')

    def is_valid_time_format(s: str) -> bool:
        return bool(pattern.match(s))

    if ray_averaging_time != None:    
        if not is_valid_time_format(ray_averaging_time):        
            raise ConfigError(f"The provided ray_averaging_time format was not understood: {ray_averaging_time} Allowed formats: 'xmin' or 'xH' where x is an up to 2 digit positive integer that corresponds either the number of integers or the number of hours. For exacmple: 10min or 3H")
    else:
        parser_args["ray_averaging_time"] = 'all'
        
    if ray_averaging_time != None:    
        if not is_valid_time_format(ray_qck_averaging_time):        
            raise ConfigError(f"The provided ray_qck_averaging_time format was not understood: {ray_qck_averaging_time} Allowed formats: 'xmin' or 'xH' where x is an up to 2 digit positive integer that corresponds either the number of integers or the number of hours. For exacmple: 10min or 3H")
    else:
        parser_args["ray_qck_averaging_time"] = 'raw'
    
    rsonde_column_units = parser_args.get("rsonde_column_units")
    
    if rsonde_column_units[0] not in height_units:
        raise ConfigError(f"The height units provided in rsonde_column_units are wrong {rsonde_column_units} Allowed values: {height_units}")
    if rsonde_column_units[1] not in pressure_units:
        raise ConfigError(f"The pressure units provided in rsonde_column_units are wrong {rsonde_column_units} Allowed values: {pressure_units}")
    if rsonde_column_units[2] not in temperature_units:
        raise ConfigError(f"The temperature units provided in rsonde_column_units are wrong {rsonde_column_units} Allowed values: {temperature_units}")
    if rsonde_column_units[3] not in humidity_units:
        raise ConfigError(f"The humidity units provided in rsonde_column_units are wrong {rsonde_column_units} Allowed values: {humidity_units}")
    
    slice_measurement = parser_args.get("slice_measurement")
    exclude_measurement = parser_args.get("exclude_measurement")

    def wrong_format(name: str, identifier: str, value: str):
        wrong_format_text = f"The format of the provided {name} {identifier} {value} is wrong. Please use the hhmm format where hh ranges from 00 to 23 and mm ranges from 00 to 59."
        return(wrong_format_text)
            
    example_text = "Only sets of 3 strings are acceptable where the first one is the test identifier (e.g. 'tlc'), the second is the start time in hhmm format (e.g. 2330) and the third is the stop time in hhmm format (e.g. 0100)"
            
    if len(slice_measurement) > 0:
        if len(slice_measurement) % 3 != 0:
            raise ConfigError(f"The provided slice_measurement parameter is wrong {slice_measurement}\n{example_text}")

        for i in range(0, len(slice_measurement), 3):
            if slice_measurement[i] not in recognised_slicers:
                raise ConfigError(f"slice_measurement: unrecognized qa_test ID provided: {slice_measurement[i]}\nRecognised IDs: {recognised_slicers}")                

        for i in range(1, len(slice_measurement), 3):
            if len(slice_measurement[i]) != 4:
                raise ConfigError(wrong_format("slice_measurement", "start time", slice_measurement[i]))
            elif slice_measurement[i][:2] not in [str(x).zfill(2) for x in np.arange(0,24,1)]:
                raise ConfigError(wrong_format("slice_measurement", "start time", slice_measurement[i]))
            elif slice_measurement[i][2:] not in [str(x).zfill(2) for x in np.arange(0,59,1)]:
                raise ConfigError(wrong_format("slice_measurement", "start time", slice_measurement[i]))
                
        for i in range(2, len(slice_measurement), 3):
            if len(slice_measurement[i]) != 4:
                raise ConfigError(wrong_format("slice_measurement", "stop time", slice_measurement[i]))
            elif slice_measurement[i][:2] not in [str(x).zfill(2) for x in np.arange(0,24,1)]:
                raise ConfigError(wrong_format("slice_measurement", "stop time", slice_measurement[i]))
            elif slice_measurement[i][2:] not in [str(x).zfill(2) for x in np.arange(0,59,1)]:
                raise ConfigError(wrong_format("slice_measurement", "stop time", slice_measurement[i]))
                
    if exclude_measurement != None:
        if len(exclude_measurement) % 3 != 0:
            raise ConfigError(f"The provided exclude_measurement parameter is wrong {exclude_measurement}\n{example_text}")
        
        for i in range(0, len(exclude_measurement), 3):
            if exclude_measurement[i] not in recognised_slicers:
                raise ConfigError(f"slice_measurement: unrecognized qa_test ID provided: {exclude_measurement[i]}\nRecognised IDs: {recognised_slicers}")                
        
        for i in range(1, len(exclude_measurement), 3):
            if len(exclude_measurement[i]) != 4:
                raise ConfigError(wrong_format("exclude_measurement", "start time", exclude_measurement[i]))
            elif exclude_measurement[i][:2] not in [str(x).zfill(2) for x in np.arange(0,24,1)]:
                raise ConfigError(wrong_format("exclude_measurement", "start time", exclude_measurement[i]))
            elif exclude_measurement[i][2:] not in [str(x).zfill(2) for x in np.arange(0,59,1)]:
                raise ConfigError(wrong_format("exclude_measurement", "start time", exclude_measurement[i]))
                
        for i in range(2, len(exclude_measurement), 3):
            if len(exclude_measurement[i]) != 4:
                raise ConfigError(wrong_format("exclude_measurement", "stop time", exclude_measurement[i]))
            elif exclude_measurement[i][:2] not in [str(x).zfill(2) for x in np.arange(0,24,1)]:
                raise ConfigError(wrong_format("exclude_measurement", "stop time", exclude_measurement[i]))
            elif exclude_measurement[i][2:] not in [str(x).zfill(2) for x in np.arange(0,59,1)]:
                raise ConfigError(wrong_format("exclude_measurement", "stop time", exclude_measurement[i]))
    
    CustomWarning("Please note that since update 0.6.0 the more generic slice_measurement parameter can be provided instead of slice_rayleigh. The slice_rayleigh parameter will be ignored even if provided")



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

def _main_data_folder_check(parser_args: Dict[str, Any], filepath: str) -> None:
    
    main_data_folder = parser_args.get("main_data_folder")
    
    if main_data_folder is None:
        print("-- The main_data_folder was not provided. By default it will be assigned to the parent folder of the call_atlas.ini file.\n")
        parser_args['main_data_folder'] = os.path.dirname(filepath)

    else:
        for key in explicit_path_keys:
            if parser_args.get(key) != None:
                CustomWarning(f"Both {key} and the main_data_folder were provided. Automatically detected paths are always overridden by explicit paths")
    return(parser_args)
    
def _absolute_paths_exist_check(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    
    for key, meta in SCHEMA.items():
        path = parser_args.get(key)
        if meta.get("check_path") in ["dir", "file"] and path != None:
            path = os.path.normpath(path)
            if not os.path.exists(path):
                raise ConfigError(f"{key} is provided but does not point to an existing path: {path}")
                
            if meta.get("check_path") == "dir":
                if not os.path.isdir(path):
                    raise ConfigError(f"{key} is provided but does not point to an existing directory: {path}")
            
            if meta.get("check_path") == "file":
                if not os.path.isfile(path):
                    raise ConfigError(f"{key} is provided but does not point to an existing file: {path}")
            
            parser_args[key] = path
    
    return parser_args

def _relative_paths_exist_check(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    
    parent_folder = parser_args.get("parent_folder")
    if parent_folder == None:
        raise ConfigError("parent_folder is empty. It should have been already assigned at this stage")
 
    for key, meta in SCHEMA.items():
        rel_path = parser_args.get(key)
        if meta.get("check_path") == "relative":
            parser_args[f"abs_{key}"] = None
            if key != "drk":
                parser_args[f"abs_drk_{key}"] = None
            if rel_path != None:
                path = os.path.normpath(os.path.join(parent_folder, rel_path))
                path_drk = os.path.normpath(os.path.join(parent_folder, f"drk_{rel_path}"))
                if not os.path.exists(path):
                    raise ConfigError(f"{key} is provided but does not point to an existing path: {path}")
                if not os.path.isdir(path):
                    raise ConfigError(f"{key} is provided but does not point to an existing directory: {path}")
                parser_args[f"abs_{key}"] = path
                if key == 'drk':
                    continue
                if os.path.exists(path_drk):
                    parser_args[f"abs_drk_{key}"] = path_drk
                    
            else:
                path = os.path.normpath(os.path.join(parent_folder, key))
                path_drk = os.path.normpath(os.path.join(parent_folder, f"drk_{key}"))
                if os.path.exists(path):
                    parser_args[f"abs_{key}"] = path
                if key == 'drk':
                    continue
                if os.path.exists(path_drk):
                    parser_args[f"abs_drk_{key}"] = path_drk
                        
    return parser_args

def _rename_folder(parser_args: Dict[str, Any], old_key: str, new_key: str, new_name: str, version_warning: str = "") -> Dict[str, Any]:
    
    if parser_args[old_key] != None and parser_args[new_key] != None:
        raise ConfigError(f"{old_key} and {new_key} folders cannot be present at the same time. {version_warning}")

    elif parser_args[old_key] != None and parser_args[new_key] == None:
        CustomWarning(f"folder with the {old_key} suffix was detected. {version_warning}")

        old = Path(parser_args[old_key])
        new = old.parent / new_name
        
        old.rename(new)    
    
        parser_args[new_key] = new
    
    del parser_args[old_key]

    return(parser_args)

def _special_path_handling(parser_args: Dict[str, Any]) -> Dict[str, Any]:
    
    version_warning_nrm = "\nSince update 0.6.0 the ray suffix was introduced as a replacement for nrm for Rayleigh measurements. Any folders with the nrm suffix will be renamed using the ray suffix\n"
    version_warning_pcb = "\nSince update 0.6.0 the p45 and m45 folders have been introduced as replacements to the +45 and -45 folders. Any subfolders named +45 (-45) will be automatically renamed to p45 (m45), respectively\n"

    # Rename any existing folders with the nrm suffix to ray
    parser_args = _rename_folder(parser_args, old_key = "abs_nrm", new_key = "abs_ray", new_name = "ray", version_warning = version_warning_nrm)
    parser_args = _rename_folder(parser_args, old_key = "abs_drk_nrm", new_key = "abs_drk_ray", new_name = "drk_ray", version_warning = version_warning_nrm)
    parser_args = _rename_folder(parser_args, old_key = "abs_nrm_pcb", new_key = "abs_ray_pcb", new_name = "ray_pcb", version_warning = version_warning_nrm)
    parser_args = _rename_folder(parser_args, old_key = "abs_drk_nrm_pcb", new_key = "abs_drk_ray_pcb", new_name = "drk_ray_pcb", version_warning = version_warning_nrm)

    # ray and ray_pcb measurements will be the same if only one of them is provided
    path_drk = parser_args.get("abs_drk")
    path_ray = parser_args.get("abs_ray")
    path_ray_pcb = parser_args.get("abs_ray_pcb")
    path_pcb = parser_args.get("abs_pcb")
    
    if path_ray == None and path_ray_pcb != None:
        parser_args["abs_ray"] = path_ray_pcb
        parser_args["abs_ray_drk"] = parser_args["abs_drk_ray_pcb"]
    
    if path_ray_pcb == None and path_ray != None and path_pcb != None:
        parser_args["abs_ray_pcb"] = path_ray
        parser_args["abs_drk_ray_pcb"] = parser_args["abs_drk_ray"]
        
    # use common drk measurement if not provided explicitely for a QA test
    for key in qa_tests:
        if key == 'drk':
            continue
        if parser_args.get(f"abs_drk_{key}") == None and path_drk != None and parser_args.get(f"abs_{key}") != None:
            parser_args[f"abs_drk_{key}"] = path_drk
    
    if parser_args["quick_run"]:
        timestamp = 'quick_run'
    else:
        now = datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S")
    if parser_args.get("output_folder") is None:
        fpath = parser_args["parent_folder"]
        parser_args["output_folder"] = os.path.join(fpath, f"analysis_{timestamp}")
        
    else:
        fname = f"{os.path.basename(parser_args['parent_folder'])}_{timestamp}"
        parser_args["output_folder"] = os.path.join(parser_args["output_folder"], fname)
        
    os.makedirs(parser_args["output_folder"], exist_ok = True)

    # Create the tlc, tlc_rin, pcb, and pcb_aux related subolders
    _expand_subfolders(parser_args, base_key = "abs_tlc", subfolders = tlc_subfolders)
    _expand_subfolders(parser_args, base_key = "abs_tlc_rin", subfolders = tlc_rin_subfolders)
    _expand_subfolders(parser_args, base_key = "abs_pcb", subfolders = pcb_subfolders)
    _expand_subfolders(parser_args, base_key = "abs_pcb_aux", subfolders = pcb_subfolders)

    # Rename any existing +45 and -45 folders to p45 and m45
    parser_args = _rename_folder(parser_args, old_key = "abs_pcb_+45", new_key = "abs_pcb_p45", new_name = "p45", version_warning = version_warning_pcb)
    parser_args = _rename_folder(parser_args, old_key = "abs_pcb_-45", new_key = "abs_pcb_m45", new_name = "m45", version_warning = version_warning_pcb)
    parser_args = _rename_folder(parser_args, old_key = "abs_pcb_aux_+45", new_key = "abs_pcb_aux_p45", new_name = "p45", version_warning = version_warning_pcb)
    parser_args = _rename_folder(parser_args, old_key = "abs_pcb_aux_-45", new_key = "abs_pcb_aux_m45", new_name = "m45", version_warning = version_warning_pcb)
           
    all_none = all([parser_args[key] == None for key in parser_args.keys() if key.startswith("abs")])
    if all_none:
        endpoint(1)
            
    return parser_args

def _expand_subfolders(
    d: dict,
    base_key: str,
    subfolders: Iterable[str] | Mapping[str, str],
    *,
    must_exist: bool = True,          # if True, put None if subfolder doesn't exist
    normalize_to_path: bool = False,   # convert values to pathlib.Path
    delete_base: bool = True,        # remove the original base_key after expansion
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

    # Normalize subfolders into (suffix, foldername) pairs
    items = subfolders.items() if isinstance(subfolders, Mapping) \
           else ((name, name) for name in subfolders)

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
    pattern: str = "*",             # which files to consider (glob pattern)
    overwrite: bool = False,        # if True, replace existing files at destination
) -> None:
    """
    If cfg[files_per_key] is not None, distribute files from cfg[base_key] into subfolders
    named <sector> (e.g., north/east/south/west) in sequential blocks of size files_per_sector.
    Enforces:
      - base folder must exist (when files_per_sector is not None)
      - base folder must not be empty (no files) when files_per_sector is not None
      - total files must be divisible by files_per_sector * len(sectors)
    Moves files in the order of sorted names to keep behavior deterministic.
    """

    files_per_sector = cfg.get(files_per_key, None)
    if files_per_sector is None:
        return  # nothing to do

    # Validate base folder
    base_val = cfg.get(base_key, None)
    if base_val is None:
        raise DistributionError(
            f"{base_key} is None but {files_per_key} is set to {files_per_sector}."
        )

    base = Path(base_val)
    if not base.exists() or not base.is_dir():
        raise DistributionError(f"{base_key} points to a non-existent directory: {base}")

    # Gather files in the base folder (exclude directories)
    files = sorted([p for p in base.glob(pattern) if p.is_file() and p.parent == base])

    if not files:
        raise DistributionError(
            f"No files detected in {base} (pattern='{pattern}') while {files_per_key}={files_per_sector}."
        )

    k = len(sectors)
    group = files_per_sector * k
    if len(files) % group != 0:
        raise DistributionError(
            f"The {files_per_key} was provided but the file count {len(files)} in {base} is not divisible by {files_per_key} * folders "
            f"({files_per_sector} * {k} = {group})."
        )

    # Ensure sector subfolders exist
    dest_dirs = []
    for sect in sectors:
        d = base / sect
        d.mkdir(parents=True, exist_ok=True)
        dest_dirs.append(d)

    # Move files in blocks: first N -> north, next N -> east, next N -> south, next N -> west, then repeat
    for i, f in enumerate(files):
        sector_idx = (i // files_per_sector) % k
        dest_dir = dest_dirs[sector_idx]
        target = dest_dir / f.name
        if overwrite:
            f.replace(target)     # overwrites if exists
        else:
            if target.exists():
                raise DistributionError(f"Destination already has file: {target}")
            f.rename(target)

    # Optionally: record the sector paths back into cfg (handy later)
    for sect, d in zip(sectors, dest_dirs):
        cfg[f"{base_key}_{sect}"] = d

def assert_pairs_unique(recorder_channel_id, laser_id):
    if len(recorder_channel_id) != len(laser_id):
        raise ConfigError("Duplicate recorder_channel_id values found {recorder_channel_id}. laser_id is mandatory in this case")
    seen = set()
    dups = []
    for i, pair in enumerate(zip(recorder_channel_id, laser_id)):
        if pair in seen:
            dups.append((i, pair))
        else:
            seen.add(pair)
    if dups:
        raise ConfigError(f"Duplicate recorder_channel_id values found {recorder_channel_id}. laser_id was provided {laser_id} but it is the same for at least one of the duplicate recorder_channel_id values")
        
def read_ini_file(filepath: str)  -> Dict[str, Any]:
    
    """Read, convert, expand from scalar defaults, compute simple defaults, validate, and return dict."""
    config = configparser.ConfigParser(allow_no_value=True, strict=True)

    read_files = config.read(filepath, encoding="utf-8")

    if not read_files:
        raise ConfigError(f"INI file not found or unreadable: {filepath}\nMake sure the encoding is utf-8")

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

                  
        if found == False:
            if meta["is_list"]:
                parser_args[key] = []
            else:
                parser_args[key] = meta["default"]
    
    return parser_args

def _get_mtype(d:Dict[str, Any]) -> Dict[str, Any]:
    
    to_add = {}
    for key in d.keys():
        if key.startswith("abs_drk"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "drk"
        elif key.startswith("abs_ray"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "nrm"
        elif key.startswith("abs_tlc"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "tlc"
        elif key.startswith("abs_pcb_p45"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "pcb_p45"
        elif key.startswith("abs_pcb_m45"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "pcb_m45"
        elif key.startswith("abs_trg"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "nrm"
        elif key.startswith("abs_dtm"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "nrm"
        elif key.startswith("abs_cam"):
            to_add[f"mtype_{key.removeprefix('abs_')}"] = "cam"
        
    d.update(to_add)

def _export_hoi_cfg_check(parser_args: Dict[str, Any]) -> None:
    export_hoi_cfg = parser_args["export_hoi_cfg"]
    atlas_configuration_file = parser_args["atlas_configuration_file"]
    if export_hoi_cfg in ["1", "2"] and atlas_configuration_file != None:
        raise ConfigError(f"export_hoi_cfg was set to {export_hoi_cfg} but atlas_configuration_file was provided in the explicit_paths section. Either change export_hoi_cfg to 0 or leave atlas_configuration_file empty")
    
# -------------------------------------------------------------------
# Public API
# -------------------------------------------------------------------

def parse_call_atlas_ini(filepath: str, debug: bool = False) -> Dict[str, Any]:

    print_header(f"Parsing the ATLAS initialization file\n{filepath}")
    
    # 1) Parse & convert
    parser_args = read_ini_file(filepath)         
      
    # 2) Enforce mandatory/recommended
    _enforce_mandatory_and_recommended(parser_args)

    # 3) Check if the main_data_folder exists, use any explicitely defined path instead of trying to automatically detect it later 
    parser_args = _main_data_folder_check(parser_args, filepath)

    # 4) Fill with default values if empty
    parser_args = _fill_with_defaults(parser_args)
    
    # 5) Check if export_hoi_cfg is conflicting atlas_configuration_file
    _export_hoi_cfg_check(parser_args)
    
    # 6) Check if all provided paths exist
    parser_args = _absolute_paths_exist_check(parser_args)
    
    # 7) Fill in the explicit paths if they are not available with the autotedected paths
    parser_args = autodetect_paths(parser_args)
    
    # 8) Check if the relative QA test folder paths exists and add the corresponding absolute paths to the dictionary - add also the corresponding dark paths per test
    parser_args = _relative_paths_exist_check(parser_args)
    
    # 9) Distribute tlc files in the correct sector if files_per_sector or files_per_sring are provided
    _distribute_files_into_sectors(parser_args, base_key = 'abs_tlc', files_per_key = 'files_per_sector', sectors = tlc_subfolders)
    _distribute_files_into_sectors(parser_args, base_key = 'abs_tlc_rin', files_per_key = 'files_per_sring', sectors = tlc_rin_subfolders)

    # 10) Fill absolute paths from provided relative paths and existing folders
    parser_args = _special_path_handling(parser_args)
    
    # 11) Set the measurement type according to the absolute paths 
    _get_mtype(parser_args)
    
    # 13) Range & allowed checks (only on non-empty values)
    for name, meta in SCHEMA.items():
        _check_size(name, parser_args.get(name), meta)
        _check_limits(name, parser_args.get(name), meta)
        _check_allowed(name, parser_args.get(name), meta)

    # 14) Special handling - check the format of slice_measurement and exclude_measurement
    _special_checks(parser_args)
    
    # 15) Sort keys alphabetically
    parser_args = dict(sorted(parser_args.items()))

    if debug == True:
        pprint(parser_args)
    
    return parser_args

