#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:22:31 2022

@author: nick
"""

import numpy as np
import sys
import datetime as dt
import xarray as xr
from helper_functions.printouts import print_header
from itertools import combinations
import pandas as pd
from typing import Any, Dict, List, Tuple
from utils.__parse_config_file__ import _compute_background_bins_if_missing, \
    _compute_emitted_wavelength_if_missing, _compute_dead_time_if_missing, \
        _compute_GH_if_missing
from utils.__parse_config_file__ import SCHEMA, SYSTEM_KEYS, CHANNEL_KEYS, POL_CAL_KEYS
from utils.error_classes import ConfigError
from utils.time_conversions import iso_to_datetimes

def unique_cross_key_pairs(d1: Dict[str, Any], d2: Dict[str, Any]):
    keys = d1.keys()  # same as d2.keys()
    for k1, k2 in combinations(keys, 2):
        yield (k1, k2)

def assert_metadata(d: Dict[str, Any], key1: str, key2: str, sec: str):
    v1 = d[key1][sec]
    v2 = d[key2][sec]
    if not v1.equals(v2):
        print(f"--Warning: {key1} and {key2} files have diferent metadata.")
        
def hhmm_to_datetime(hhmm: str, base: pd.Timestamp) -> pd.Timestamp:
    """Convert hhmm string to datetime on the date of base timestamp."""
    return pd.Timestamp(year=base.year,
                        month=base.month,
                        day=base.day,
                        hour=int(hhmm[:2]),
                        minute=int(hhmm[2:]))

def make_interval(start_str: str, stop_str: str, base: pd.Timestamp):
    start_dt = hhmm_to_datetime(start_str, base)
    stop_dt = hhmm_to_datetime(stop_str, base)

    # If stop is "before" start → assume it's the next day
    if stop_dt <= start_dt:
        stop_dt += pd.Timedelta(days=1)

    return start_dt, stop_dt

def load_metadata(config_info: Dict[str, Any], 
                  metadata: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:

    print_header("Loading file metadata to configuration")

    for (key1, key2) in list(unique_cross_key_pairs(metadata,metadata)):
        assert_metadata(metadata, key1 = key1, key2 = key2, sec = "system_info")
        assert_metadata(metadata, key1 = key1, key2 = key2, sec = "channel_info")
        
    cfg = {}
    for meas_key in metadata.keys():
        
        cfg[meas_key] = config_info.copy()
        
        src_sys = metadata[meas_key]["system_info"]
        src_chs = metadata[meas_key]["channel_info"]
        
        ext_ids = cfg[meas_key]["recorder_channel_id"]
        int_ids = src_chs.index.values
        
        if not set(ext_ids).issubset(int_ids):
            missing_ind = [i for i, v in enumerate(ext_ids) if v not in int_ids]
            missing_ids = [ext_ids[i] for i, x in enumerate(ext_ids) if i in missing_ind]
            print(f"{meas_key} dataset")
            print(f"--Warning: The provided recorder_channel_id is not a subset of the raw file header IDs. Unrecognised channels will be excluded: \n-- unrecognised recorder_channel_id \n{missing_ids}\n--header_channel_id \n{int_ids}\n")
            available_ids = [ext_ids[i] for i, x in enumerate(ext_ids) if i not in missing_ind]
            if len(available_ids) == 0:
                sys.exit("Endpoint 2: No channels to process! Check the provided channel_recorder_id values. ATLAS terminates here")
        else:
            missing_ind = []
            available_ids = ext_ids
                
        for key in config_info.keys():

            val = config_info[key]
            
            if val == None and not src_sys.empty:
                if key in src_sys.index:
                    if src_sys.loc[key] != None:
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
                print(f"--Warning: Datasets {key1} and {key2} contain a different number of channels. The configuration parameters will be taken from the dataset with the most channels:\n--{key1}: {rec1}\n--{key2}: {rec2}")
                
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
        
        metadata[key]["shots"] = metadata[key]["shots"].loc[dict(channel = recorder_channel_id)]
            
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
                        zip(wvl, tel_type, ch_type, ch_stype, acq_mode)]
    
    config_info["atlas_channel_id"] = atlas_channel_id
    
    for key in profiles.keys():
        signals = profiles[key].copy()
        shots = metadata[key]["shots"].copy()
        
        signals = signals.assign_coords(channel = atlas_channel_id)
        shots = shots.assign_coords(channel = atlas_channel_id)
        
        profiles[key] = signals.copy()
        metadata[key]["shots"] = shots.copy()
    
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
    
    necesary_parameters = ["recorder_channel_id", 
                           "acquisition_mode",
                           "range_resolution",
                           "laser_repetition_rate"]
    
    # Check if the data aqcuisition range is empty for any analogue channel
    for i, ch in enumerate(config_info["recorder_channel_id"]):
        for key in necesary_parameters:
            val = config_info[key][i]
            if val == None:
                raise Exception(f"{key} is empty for ch {ch}: {val}")

    for i, ch in enumerate(config_info["recorder_channel_id"]):
        daq_range = config_info["data_acquisition_range"][i]
        aq_mode = config_info["acquisition_mode"][i]
        
        if aq_mode == "a" and (daq_range in [None, np.inf, np.nan]):
            raise Exception(f"data_acquisition_range is empty for analogue ch {ch}: {daq_range}")
        
    # Fill in the dead_time if missing
    config_info = _compute_dead_time_if_missing(config_info)       

    # Fill in the emitted_wavelength values based on the detected_wavelength if missing
    config_info = _compute_emitted_wavelength_if_missing(config_info)       

    # Fill in the background_low_bin and background_high_bin values if they are empty based on the zero bin
    config_info = _compute_background_bins_if_missing(config_info)       
            
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
                
def store_updated_metadata(config_info: Dict[str, Any], 
                           metadata: Dict[str, Any]) -> Dict[str, Any]:
    
    system_info  = pd.Series()
    channel_info = pd.DataFrame(index = config_info["atlas_channel_id"])
    pol_cal_info = pd.DataFrame(index = config_info["pol_cal_pairs"])

    for key in config_info:
        if key in SYSTEM_KEYS:
            system_info.loc[key] = config_info[key]
        elif key in CHANNEL_KEYS:
            channel_info.loc[:,key] = [v for v in config_info[key]]
        elif key in POL_CAL_KEYS:
            pol_cal_info.loc[:,key] = [v for v in config_info[key]]
        
    for key in metadata.keys():

        time_info = metadata[key]["time_info"]
        
        metadata[key]["system_info"] = xr.DataArray(system_info.values,
                                                    dims = ["parameters"],
                                                    coords = [system_info.index.values])
        metadata[key]["channel_info"] = xr.DataArray(channel_info.T.values,
                                                     dims = ["parameters", "channel"],
                                                     coords = [channel_info.columns.values, channel_info.index.values])
        
        metadata[key]["time_info"] = xr.DataArray(time_info.T.values,
                                                  dims = ["parameters", "time"],
                                                  coords = [time_info.columns.values, time_info.index.values])
        
        metadata[key]["pol_cal_info"] = xr.DataArray(pol_cal_info.T.values,
                                                     dims = ["parameters", "pairs"],
                                                     coords = [pol_cal_info.columns.values, pol_cal_info.index.values])
        
        
    return(metadata)


def slice_and_exclude(caller_info: Dict[str, Any], profiles: Dict[str, Any], 
                      metadata: Dict[str, Any], profile_masks: Dict[str, Any]) -> \
    Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    
    print_header("Slicing and excluding measurement parts")

    mask_time = {}
        
    slicer = caller_info["slice_measurement"]
    debug_signals = caller_info["debug_signals"]

    if slicer != None:
        slice_meas = [slicer[i] for i in range(0,len(slicer),3)]
        slice_start_time = [slicer[i] for i in range(1,len(slicer),3)]
        slice_stop_time = [slicer[i] for i in range(2,len(slicer),3)]
    
    for i, key in enumerate(slice_meas):
        if key in profiles.keys():
            time_info = metadata[key]["time_info"]
            start_times = iso_to_datetimes(time_info.sel({"parameters":"start_time"}).values)
            
            start_dt, stop_dt =  make_interval(
                start_str = slice_start_time[i], 
                stop_str = slice_stop_time[i], 
                base = start_times[0]
                )
            
            time = profiles[key]['time']
            mask_interval = (time >= start_dt) & (time <= stop_dt)
            if key in mask_time.keys():
                mask_time[key] = (mask_time[key] | mask_interval)
            else:
                mask_time[key] = mask_interval
    
    slicer = caller_info["exclude_measurement"]
    if slicer != None:
        slice_meas = [slicer[i] for i in range(0,len(slicer),3)]
        slice_start_time = [slicer[i] for i in range(1,len(slicer),3)]
        slice_stop_time = [slicer[i] for i in range(2,len(slicer),3)]
    
    for i, key in enumerate(slice_meas):
        if key in profiles.keys():
            time_info = metadata[key]["time_info"]
            start_times = iso_to_datetimes(time_info.sel({"parameters":"start_time"}).values)

            start_dt, stop_dt =  make_interval(
                start_str = slice_start_time[i], 
                stop_str = slice_stop_time[i], 
                base = start_times[0]
                )
            
            time = profiles[key]['time']
            mask_interval = (time < start_dt) | (time > stop_dt)
            if key in mask_time.keys():
                mask_time[key] = (mask_time[key] & mask_interval)
            else:
                mask_time[key] = mask_interval
    
    del_keys = []
    for key in mask_time.keys():
        if not mask_time[key].any().values:
            del_keys.append(key)
        else:
            profiles[key] = profiles[key].where(mask_time[key]).copy()
            if debug_signals:
                profile_masks[key]["sliced"] = mask_time[key]

    
    for key in del_keys:
        print(f"--Warning: {key} measurement will not be processed because the combination of the following parameters remove all profiles:\n\n  slice_measurement: {caller_info['slice_measurement']}\n\n  exclude_measurement: {caller_info['exclude_measurement']}\n")
        del profiles[key]
        del metadata[key]
        del profile_masks[key]

    if profiles == {}:
        sys.exit("Endpoint 3: Slicing and excluding measurement parts removed all measurements. No signals to process. ATLAS terminates here")
            
    return(profiles, metadata, profile_masks)

def screen_low_shots(profiles: Dict[str, Any],  metadata: Dict[str, Any], 
                     profile_masks: Dict[str, Any], shot_limit: float = 0.9) -> \
    Dict[str, Any]:

    """Replaces values in all bins with nans if the numer of shots is
    lower than 10% of the maximum """
    
    shot_limit_percent = round(100. *  shot_limit)
    print_header(f"Removing profiles with too few shots ({shot_limit_percent}%)")

    for key in profiles.keys():
        
        signal = profiles[key].copy()
        filename = metadata[key]["time_info"].sel({"parameters": "filename"})
        shots = metadata[key]["shots"]
               
        mask_low_shots = xr.where((shots / shots.max("time")  < shot_limit) & \
            (shots.max("time") > 20), True, False).compute()
            
        any_low_shots = mask_low_shots.any("time")
            
        ch_low_shots = any_low_shots["channel"].values[any_low_shots.values]
        
        if len(ch_low_shots) > 0:
            print(f'{key} dataset')
            print(f"-- Warning: The following files have less shots than {shot_limit_percent}% of the naximum number of shots encountered and will be screened out:")
                          
            for ch in ch_low_shots:
                
                mask_time = mask_low_shots.sel({"channel" : ch}).values
                time_low_shots = mask_low_shots["time"][mask_time]
                
                for t in time_low_shots:
                    print(f'    Filename: {filename.loc[t].values} | Channel: {ch}')
            print("")
        
        profile_masks[key]["low_shots"] = mask_low_shots

        signal = signal.where(~mask_low_shots, np.nan)

        profiles[key] = signal.copy()
            
    return(profiles, profile_masks)

def initialize_debug_dictionaries(profiles: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    
    print_header("Initializing debugging-related data arrays")

    profile_masks = {}
    profile_db = {}
    
    for key in profiles.keys():
        profile_masks[key] = {}
        profile_db[key] = {}

    return(profile_masks, profile_db)
    
    