#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:22:31 2022

@author: nick
"""

import numpy as np
import xarray as xr
from utils.printouts import print_header, endpoint
from itertools import combinations
import pandas as pd
from typing import Any, Dict, Tuple
from utils.__parse_config_file__ import _compute_background_bins_if_missing, \
    _compute_emitted_wavelength_if_missing, _compute_dead_time_if_missing, \
        _compute_GH_if_missing, _check_background_limits, _normalize_entries
from utils.__parse_config_file__ import SCHEMA, SYSTEM_KEYS, CHANNEL_KEYS, POL_CAL_KEYS
from utils.error_classes import ConfigError, CustomWarning

def unique_cross_key_pairs(d1: Dict[str, Any], d2: Dict[str, Any]):
    keys = d1.keys()  # same as d2.keys()
    for k1, k2 in combinations(keys, 2):
        yield (k1, k2)

def assert_metadata(d: Dict[str, Any], key1: str, key2: str):
    v1 = d[key1]
    v2 = d[key2]
    if not v1.equals(v2):
        CustomWarning(f"{key1} and {key2} files have diferent metadata.")

def load_metadata(config_info: Dict[str, Any], 
                  metadata: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:

    print_header("Loading file metadata to configuration")

    unique_system_info = list(unique_cross_key_pairs(metadata["system_info"], metadata["system_info"]))
    unique_channel_info = list(unique_cross_key_pairs(metadata["channel_info"], metadata["channel_info"]))
    
    for (key1, key2) in unique_system_info:
        assert_metadata(metadata["system_info"], key1 = key1, key2 = key2)

    for (key1, key2) in unique_channel_info:
        assert_metadata(metadata["channel_info"], key1 = key1, key2 = key2)
        
    cfg = {}
    for meas_key in metadata["time_info"].keys():
        
        cfg[meas_key] = config_info.copy()
        
        src_sys = metadata["system_info"][meas_key]
        src_chs = metadata["channel_info"][meas_key]
        
        ext_ids = cfg[meas_key]["recorder_channel_id"]
        int_ids = src_chs.index.values
        
        print(f"Dataset: {meas_key}")
        if not set(ext_ids).issubset(int_ids):
            missing_ind = [i for i, v in enumerate(ext_ids) if v not in int_ids]
            missing_ids = [ext_ids[i] for i, x in enumerate(ext_ids) if i in missing_ind]
            CustomWarning(f"The provided recorder_channel_id is not a subset of the raw file header IDs. Unrecognised channels will be excluded: \n-- unrecognised recorder_channel_id \n{missing_ids}\n--header_channel_id \n{int_ids}\n")
            available_ids = [ext_ids[i] for i, x in enumerate(ext_ids) if i not in missing_ind]
            if len(available_ids) == 0:
                endpoint(2)
        else:
            missing_ind = []
            available_ids = ext_ids
                
        for key in config_info.keys():

            val = config_info[key]
            
            if val is None and not src_sys.empty:
                if key in src_sys.index:
                    if src_sys.loc[key] is not None:
                        print(f"Loading: {key}")
                        cfg[meas_key][key] = str(src_sys.loc[key])

            if isinstance(val, list):
                cfg[meas_key][key] = [v for i, v in enumerate(val) if i not in missing_ind]

                if all([x is None for x in val]) and not src_chs.empty and key in src_chs.columns:
                    print(f"Loading: {key}")
                    cfg[meas_key][key] = [str(src) if src is not None else None for src in src_chs.loc[available_ids,key].values]

    rec  = config_info["recorder_channel_id"].copy()
    for (key1, key2) in list(unique_cross_key_pairs(cfg,cfg)):
        if cfg[key1] != {} and cfg[key2] != {}:
            rec1 = cfg[key1]["recorder_channel_id"]
            rec2 = cfg[key2]["recorder_channel_id"]
            
            if rec1 != rec2:
                CustomWarning(f"Datasets {key1} and {key2} contain a different number of channels. The configuration parameters will be taken from the dataset with the most channels:\n--{key1}: {rec1}\n--{key2}: {rec2}")
                
                if (len(rec1) > len(rec2)) and (len(rec1) >= len(rec)):
                    config_info = cfg[key1]
                elif (len(rec2) > len(rec1)) and (len(rec2) >= len(rec)):
                    config_info = cfg[key2]
            
            if rec1 == rec2:
                config_info = cfg[key2]
                
    return(config_info)
   
def remove_unrecognised_channels(caller_info: Dict[str, Any], config_info: Dict[str, Any], 
                                 profiles: Dict[str, Any], metadata: Dict[str, Any]) -> Tuple[Dict[str, Any],Dict[str, Any]]:
    
    print_header("Removing unrecognised recorder channel IDs")

    recorder_channel_id = config_info["recorder_channel_id"]
            
    for key in profiles.keys():

        profiles[key] = profiles[key].loc[dict(channel = recorder_channel_id)]
        
        metadata["shots"][key] = metadata["shots"][key].loc[dict(channel = recorder_channel_id)]
            
    return(profiles, metadata)

def get_atlas_channel_id(config_info: Dict[str, Any], profiles: Dict[str, Any], 
                         metadata: Dict[str, Any]) -> Dict[str, Any]:

    print_header("Creating altas channel ID")

    dtw = config_info["detected_wavelength"]
    tel_type = config_info["telescope_type"]
    ch_type = config_info["channel_type"]
    ch_stype = config_info["channel_subtype"]
    acq_mode = config_info["acquisition_mode"]
    
    wvl = [str(round(float(w))).zfill(4) for w in dtw]
    
    atlas_channel_id = [f"{w}{t}{c}{a}{s}" for w, t, c, a, s in 
                        zip(wvl, tel_type, ch_type, acq_mode, ch_stype)]
    
    config_info["atlas_channel_id"] = atlas_channel_id
    
    for key in profiles.keys():
        signals = profiles[key].copy()
        shots = metadata["shots"][key].copy()
        
        signals = signals.assign_coords(channel = atlas_channel_id)
        shots = shots.assign_coords(channel = atlas_channel_id)
        
        profiles[key] = signals.copy()
        metadata["shots"][key] = shots.copy()
    
    return(config_info)


def load_pol_cal_defaults(config_info):
    
    print_header("Filling pol. cal. parameters with defaults")

    atlas_channel_id = config_info["atlas_channel_id"]
    recorder_channel_id = config_info["recorder_channel_id"]
    
    if config_info["ch_r"] == [] or config_info["ch_t"] == []:
    
    
        ch_r = [rec_id for at_id, rec_id in zip(atlas_channel_id, recorder_channel_id)
                if at_id[6] == "r"]
        ch_t = [rec_id for at_id, rec_id in zip(atlas_channel_id, recorder_channel_id) 
                if at_id[6] == "t"]
        
        ch_r_atlas = [at_id for at_id in atlas_channel_id if at_id[6] == "r"]
        ch_t_atlas = [at_id for at_id in atlas_channel_id if at_id[6] == "t"]
    
        clear_r = [f"{r[:5]}{r[7:]}" for r in ch_r_atlas]
        clear_t = [f"{t[:5]}{t[7:]}" for t in ch_t_atlas]
    
        com_index = [clear_t.index(r) for r in clear_r]
    
        ch_r = [ch_r[ind] for ind in com_index]
        ch_t = [ch_t[ind] for ind in com_index]
    
        ch_r_atlas = [ch_r_atlas[ind] for ind in com_index]
        ch_t_atlas = [ch_t_atlas[ind] for ind in com_index]
        
        config_info["ch_r"] = ch_r_atlas
        config_info["ch_t"] = ch_t_atlas
        
    else:
        ch_r = config_info["ch_r"]
        ch_t = config_info["ch_t"]
        
        index_r = [recorder_channel_id.index(r) for r in ch_r]
        index_t = [recorder_channel_id.index(t) for t in ch_t]
        
        ch_r_atlas = [atlas_channel_id[ind] for ind in index_r]
        ch_t_atlas = [atlas_channel_id[ind] for ind in index_t]

        for i, r in enumerate(ch_r_atlas):
            if r[6] != "r":
                raise ConfigError(f"Provided ch_r value: {ch_r[i]} ({r}) does not correspond to a reflected channel ")

        for i, t in enumerate(ch_t_atlas):
            if t[6] != "t":
                raise ConfigError(f"Provided ch_t: {ch_t[i]} ({t}) does not correspond to a transmitted channel ")
        
        config_info["ch_r"] = ch_r_atlas
        config_info["ch_t"] = ch_t_atlas
        

        
    pairs = [f"{r}_{t}" for r, t in zip(ch_r_atlas,ch_t_atlas)]
    
    config_info["pol_cal_pairs"] = pairs
    
    if config_info["pol_cal_pairs"] != []:
        for key in POL_CAL_KEYS:
            if key not in ["ch_r", "ch_t"]:
                meta = SCHEMA[key]
                if meta["is_list"]:
                    config_info[key] = len(config_info["pol_cal_pairs"]) * [str(meta["default"])]
                else:
                    config_info[key] = str(meta["default"])
    
    config_info = dict(sorted(config_info.items()))
    
    return(config_info)

def special_config_checks(config_info):
    
    necesary_channel_parameters = [
        "recorder_channel_id", 
        "acquisition_mode",
        "range_resolution",
        "laser_repetition_rate",
        "bins"
        ]
    
    necesary_system_parameters = [
        "station_id",
        "configuration_id",
        "station_altitude",
        "zenith_angle",
        ]
    
    # Check if the data aqcuisition range is empty for any analogue channel
    for i, ch in enumerate(config_info["recorder_channel_id"]):
        for key in necesary_channel_parameters:
            val = config_info[key][i]
            if val is None:
                raise Exception(f"Channel parameter {key} is empty for ch {ch}: {val}\n"+\
                                "Please provide the corresponding parameter in the config file if it is not provided in the raw file header")

    for key in necesary_system_parameters:
        val = config_info[key]
        if val is None:
            raise Exception(f"System parameter {key} is empty: {val}\n"+\
                            "Please provide the corresponding parameter in the config file if it is not provided in the raw file header")

    for i, ch in enumerate(config_info["recorder_channel_id"]):
        daq_range = config_info["data_acquisition_range"][i]
        aq_mode = config_info["acquisition_mode"][i]
        
        if aq_mode == "a" and (daq_range in [None, np.inf, np.nan]):
            raise Exception(f"data_acquisition_range is empty for analogue ch {ch}: {daq_range}")
        
    _normalize_entries(config_info)
    
    # Fill in the dead_time if missing
    config_info = _compute_dead_time_if_missing(config_info)       

    # Fill in the emitted_wavelength values based on the detected_wavelength if missing
    config_info = _compute_emitted_wavelength_if_missing(config_info)       

    # Fill in the background_low_bin and background_high_bin values if they are empty based on the zero bin
    config_info = _compute_background_bins_if_missing(config_info)      
    
    # Check if the background_low_bin and background_high_bin are within each channels bin range
    _check_background_limits(config_info)
            
    # Fill in the G and H factors if they are empty based on the channel type
    config_info = _compute_GH_if_missing(config_info)
    
    return(config_info)

def bring_to_correct_type(config_info):

    for key in config_info.keys():
        if key in SCHEMA:
            dtype = SCHEMA[key]["dtype"]
            entries = config_info[key]
            if dtype != str:
                if isinstance(entries, list):
                    config_info[key] = [dtype(s) if s is not None else None for s in entries]
                elif isinstance(entries, (int, float, str, bool)):
                    config_info[key] = dtype(entries)  
    
    return(config_info)
                
def store_updated_metadata(
    config_info: Dict[str, Any],
    metadata: Dict[str, Any],
) -> Dict[str, Any]:

    system_info = {}
    channel_info = pd.DataFrame(index=config_info["atlas_channel_id"])
    pol_cal_info = pd.DataFrame(index=config_info["pol_cal_pairs"])

    metadata["pol_cal_info"] = {}

    for key in config_info:

        if key in SYSTEM_KEYS:
            system_info[key] = config_info[key]

        elif key in CHANNEL_KEYS:
            channel_info.loc[:, key] = [v for v in config_info[key]]

        elif key in POL_CAL_KEYS:
            pol_cal_info.loc[:, key] = [v for v in config_info[key]]

    system_parameters = list(system_info.keys())
    system_values = np.array(
        list(system_info.values()),
        dtype=object,
    )

    for key in metadata["time_info"].keys():

        time_info = metadata["time_info"][key]

        metadata["system_info"][key] = xr.DataArray(
            system_values,
            dims=["parameters"],
            coords={"parameters": system_parameters},
        )

        metadata["channel_info"][key] = xr.DataArray(
            channel_info.T.values,
            dims=["parameters", "channel"],
            coords={
                "parameters": channel_info.columns.values,
                "channel": channel_info.index.values,
            },
        )

        metadata["time_info"][key] = xr.DataArray(
            time_info.T.values,
            dims=["parameters", "time"],
            coords={
                "parameters": time_info.columns.values,
                "time": time_info.index.values,
            },
        )

        metadata["pol_cal_info"][key] = xr.DataArray(
            pol_cal_info.T.values,
            dims=["parameters", "pairs"],
            coords={
                "parameters": pol_cal_info.columns.values,
                "pairs": pol_cal_info.index.values,
            },
        )

    return metadata
    