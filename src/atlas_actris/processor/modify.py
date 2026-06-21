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
from utils.parse_config_file import _compute_background_bins_if_missing, \
    _compute_emitted_wavelength_if_missing, _compute_dead_time_if_missing, \
        _compute_GH_if_missing, _check_background_limits, _normalize_entries
from utils.parse_config_file import (
    SCHEMA, SYSTEM_KEYS, CHANNEL_KEYS, POL_CAL_KEYS, WV_KEYS, TEMP_KEYS
    )
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
        assert_metadata(metadata["system_info"], key1=key1, key2=key2)

    for (key1, key2) in unique_channel_info:
        assert_metadata(metadata["channel_info"], key1=key1, key2=key2)

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

            CustomWarning(
                "The provided recorder_channel_id is not a subset of the raw file header IDs. "
                "Unrecognised channels will be excluded: \n"
                f"-- unrecognised recorder_channel_id \n{missing_ids}\n"
                f"--header_channel_id \n{int_ids}\n"
            )

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
                cfg[meas_key][key] = [
                    v for i, v in enumerate(val)
                    if i not in missing_ind
                ]

                if (
                    all([x is None for x in val])
                    and not src_chs.empty
                    and key in src_chs.columns
                ):
                    print(f"Loading: {key}")
                    cfg[meas_key][key] = [
                        str(src) if src is not None else None
                        for src in src_chs.loc[available_ids, key].values
                    ]

    valid_cfg_keys = [key for key, value in cfg.items() if value != {}]

    if len(valid_cfg_keys) == 1:
        return cfg[valid_cfg_keys[0]]

    rec = config_info["recorder_channel_id"].copy()

    for (key1, key2) in list(unique_cross_key_pairs(cfg, cfg)):
        if cfg[key1] != {} and cfg[key2] != {}:
            rec1 = cfg[key1]["recorder_channel_id"]
            rec2 = cfg[key2]["recorder_channel_id"]

            if rec1 != rec2:
                CustomWarning(
                    f"Datasets {key1} and {key2} contain a different number of channels. "
                    "The configuration parameters will be taken from the dataset with the most channels:\n"
                    f"--{key1}: {rec1}\n"
                    f"--{key2}: {rec2}"
                )

                if (len(rec1) > len(rec2)) and (len(rec1) >= len(rec)):
                    config_info = cfg[key1]
                elif (len(rec2) > len(rec1)) and (len(rec2) >= len(rec)):
                    config_info = cfg[key2]

            if rec1 == rec2:
                config_info = cfg[key2]

    return config_info

def find_reflected_transmitted_pairs(channels):
    """
    Find reflected/transmitted channel pairs.

    Channel ID convention:
        length = 8
        6th char, index 5: one of "t", "p", "c"
        8th char, index 7: "r" or "t"

    Pairing rule:
        r_channels = channels with 8th char "r"
        t_channels = channels with 8th char "t"

    Channels are paired when all characters are identical except
    the 6th and 8th characters.

    Returns
    -------
    r_channels : list[str]
        Reflected channels, used as numerator channels.

    t_channels : list[str]
        Transmitted channels, used as denominator channels.

    Notes
    -----
    The two returned lists are aligned by index:
        r_channels[i] pairs with t_channels[i]
    """

    allowed_6th = {"t", "p", "c"}
    allowed_8th = {"r", "t"}

    grouped = {}

    for ch in channels:
        if len(ch) != 8:
            continue

        if ch[5] not in allowed_6th:
            continue

        if ch[7] not in allowed_8th:
            continue

        # Keep chars 1-5 and 7 fixed.
        # Ignore chars 6 and 8.
        key = ch[:5] + ch[6]

        grouped.setdefault(key, {"r": [], "t": []})
        grouped[key][ch[7]].append(ch)

    r_channels = []
    t_channels = []

    for group in grouped.values():
        for ch_r in group["r"]:
            for ch_t in group["t"]:
                r_channels.append(ch_r)
                t_channels.append(ch_t)

    return r_channels, t_channels

def expand_with_loading_map(data, caller_info, *, copy=False, strict=True):
    """
    Expand a dictionary with alias keys from caller_info["loading_map"].

    If caller_info has no loading_map entry, or if loading_map is empty,
    no expansion is performed.

    Parameters
    ----------
    data : dict
        Dictionary to expand, e.g. profiles, metadata["time_info"], etc.

    caller_info : dict
        Caller information dictionary. If present and non-empty,
        caller_info["loading_map"] should map alias keys to source keys.

        Example:
            {
                "loading_map": {
                    "drk_ray": "drk",
                    "ray_pcb": "ray",
                }
            }

    copy : bool, default False
        If False, alias entries point to the same object as the source entry.
        If True, use copy.deepcopy() for each alias entry.

    strict : bool, default True
        If True, raise an error if a source key is missing or an alias key
        already exists in data.

    Returns
    -------
    expanded : dict
        New dictionary with original keys plus alias keys, or a shallow copy
        of the original dictionary if no expansion is needed.
    """

    import copy as _copy

    if not isinstance(data, dict):
        raise TypeError(f"data must be a dict. Got {type(data)}.")

    if not isinstance(caller_info, dict):
        raise TypeError(f"caller_info must be a dict. Got {type(caller_info)}.")

    loading_map = caller_info.get("loading_map")

    # No loading_map entry, None, or empty dict -> no expansion
    if not loading_map:
        return dict(data)

    if not isinstance(loading_map, dict):
        raise TypeError(
            f'caller_info["loading_map"] must be a dict. Got {type(loading_map)}.'
        )

    expanded = dict(data)

    for alias_key, source_key in loading_map.items():

        if source_key not in expanded:
            if strict:
                raise KeyError(
                    f'caller_info["loading_map"] points {alias_key!r} '
                    f"to missing source key {source_key!r}."
                )
            continue

        if alias_key in expanded:
            if strict:
                raise KeyError(
                    f'caller_info["loading_map"] alias key {alias_key!r} '
                    f"already exists in data."
                )
            continue

        if copy:
            expanded[alias_key] = _copy.deepcopy(expanded[source_key])
        else:
            expanded[alias_key] = expanded[source_key]

    return expanded

def expand_nested_with_loading_map(metadata, caller_info, *, copy=False, strict=False):
    """
    Expand second-level measurement dictionaries inside metadata.

    Expected structure
    ------------------
    metadata = {
        "system_info": {
            "ray": ...,
            "drk": ...,
        },
        "time_info": {
            "ray": ...,
            "drk": ...,
        },
        ...
    }

    Result
    ------
    metadata["system_info"]["ray_pcb"] = metadata["system_info"]["ray"]
    metadata["system_info"]["drk_ray"] = metadata["system_info"]["drk"]
    etc.

    Empty metadata groups remain empty.
    """

    if not isinstance(metadata, dict):
        raise TypeError(f"metadata must be a dict. Got {type(metadata)}.")

    loading_map = caller_info.get("loading_map") if isinstance(caller_info, dict) else None

    # No loading_map entry, None, or empty dict -> no expansion
    if not loading_map:
        return dict(metadata)

    expanded_metadata = dict(metadata)

    for group_name, group_data in metadata.items():

        # Only expand metadata groups that are dictionaries
        if not isinstance(group_data, dict):
            expanded_metadata[group_name] = group_data
            continue

        # Empty groups stay empty
        if not group_data:
            expanded_metadata[group_name] = {}
            continue

        expanded_metadata[group_name] = expand_with_loading_map(
            group_data,
            caller_info,
            copy=copy,
            strict=strict,
        )

    return expanded_metadata

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

def resolve_to_atlas_channel_id(ch, config_info):
    """
    Resolve a manually provided channel ID to atlas_channel_id.

    The provided channel may be:
        - atlas_channel_id
        - recorder_channel_id
        - scc_channel_id, if available

    Returns
    -------
    atlas_id : str
    """

    atlas_ids = config_info.get("atlas_channel_id", [])
    recorder_ids = config_info.get("recorder_channel_id", [])
    scc_ids = config_info.get("scc_channel_id", [])

    if ch in atlas_ids:
        return ch

    if ch in recorder_ids:
        return atlas_ids[recorder_ids.index(ch)]

    if scc_ids and ch in scc_ids:
        return atlas_ids[scc_ids.index(ch)]

    raise ConfigError(
        f"Provided channel '{ch}' was not found in atlas_channel_id, "
        f"recorder_channel_id, or scc_channel_id"
    )


def default_pol_cal_value(key, n_pairs):
    """
    Return the default value for a POL_CAL key.
    """

    meta = SCHEMA[key]

    if meta["is_list"]:
        return n_pairs * [str(meta["default"])]

    return str(meta["default"])


def extend_or_default_pol_cal_list(existing, default, n_target):
    """
    Keep manually provided values where available and fill the rest with defaults.

    This is used for POL_CAL list-like metadata other than ch_r/ch_t.
    """

    if existing is None:
        existing = []

    if not isinstance(existing, list):
        existing = [existing]

    out = [str(v) for v in existing[:n_target]]

    if len(out) < n_target:
        out.extend((n_target - len(out)) * [str(default)])

    return out

def check_reflected_transmitted_pairs(r_channels, t_channels):
    """
    Check reflected/transmitted channel pairs.

    Invalid pairs are skipped.
    Non-critical inconsistencies are reported as warnings.

    Parameters
    ----------
    r_channels : list[str]
        Channels expected to have last character "r".
    t_channels : list[str]
        Channels expected to have last character "t".

    Returns
    -------
    pairs : list[dict]
        Valid pairs, in the same style as find_reflected_transmitted_pairs.

    warnings : list[str]
        Warnings for skipped or suspicious pairs.
    """

    allowed_6th = {"t", "p", "c"}

    pairs = []
    warnings = []

    if len(r_channels) != len(t_channels):
        warnings.append(
            f"Different number of reflected and transmitted channels: "
            f"{len(r_channels)} vs {len(t_channels)}. "
            f"Only the first {min(len(r_channels), len(t_channels))} pairs will be checked."
        )

    for i, (ch_r, ch_t) in enumerate(zip(r_channels, t_channels)):

        pair_label = f"pair {i}: [{ch_r}, {ch_t}]"

        valid = True

        if len(ch_r) != 8:
            warnings.append(
                f"{pair_label}: skipped because reflected channel has length "
                f"{len(ch_r)}, expected 8."
            )
            valid = False

        if len(ch_t) != 8:
            warnings.append(
                f"{pair_label}: skipped because transmitted channel has length "
                f"{len(ch_t)}, expected 8."
            )
            valid = False

        if not valid:
            continue

        if ch_r[7] != "r":
            warnings.append(
                f"{pair_label}: skipped because reflected channel should end "
                f"with 'r', got '{ch_r[7]}'."
            )
            valid = False

        if ch_t[7] != "t":
            warnings.append(
                f"{pair_label}: skipped because transmitted channel should end "
                f"with 't', got '{ch_t[7]}'."
            )
            valid = False

        if ch_r[5] not in allowed_6th:
            warnings.append(
                f"{pair_label}: skipped because reflected channel has invalid "
                f"6th character '{ch_r[5]}', expected one of {sorted(allowed_6th)}."
            )
            valid = False

        if ch_t[5] not in allowed_6th:
            warnings.append(
                f"{pair_label}: skipped because transmitted channel has invalid "
                f"6th character '{ch_t[5]}', expected one of {sorted(allowed_6th)}."
            )
            valid = False

        if ch_r[:4] != ch_t[:4]:
            warnings.append(
                f"{pair_label}: skipped because wavelength differs: "
                f"'{ch_r[:4]}' != '{ch_t[:4]}'."
            )
            valid = False

        if not valid:
            continue

        if ch_r[4] != ch_t[4]:
            warnings.append(
                f"{pair_label}: warning, 5th character differs: "
                f"'{ch_r[4]}' != '{ch_t[4]}'."
            )

        if ch_r[6] != ch_t[6]:
            warnings.append(
                f"{pair_label}: warning, 7th character differs: "
                f"'{ch_r[6]}' != '{ch_t[6]}'."
            )

        pairs.append(
            {
                "numerator_channel": ch_r,
                "denominator_channel": ch_t,
                "ratio_id": f"{ch_r}_over_{ch_t}",
            }
        )

    return {
        "pairs": pairs,
        "warnings": warnings,
    }

def load_pol_cal_defaults(config_info):
    """
    Fill polarization calibration channel pairs and defaults.

    Final config_info["ch_r"] and config_info["ch_t"] are always
    atlas_channel_id.

    No pol_cal_pairs are created here.
    """

    print_header("Filling pol. cal. parameters with defaults")

    atlas_channel_id = config_info["atlas_channel_id"]

    # ------------------------------------------------------------------
    # 1. Resolve manually provided pairs to atlas_channel_id
    # ------------------------------------------------------------------
    ch_r_manual_raw = config_info.get("ch_r", []) or []
    ch_t_manual_raw = config_info.get("ch_t", []) or []

    if len(ch_r_manual_raw) != len(ch_t_manual_raw):
        raise ConfigError(
            "Provided ch_r and ch_t lists have different lengths: "
            f"{len(ch_r_manual_raw)} vs {len(ch_t_manual_raw)}"
        )

    ch_r_manual = [
        resolve_to_atlas_channel_id(ch, config_info)
        for ch in ch_r_manual_raw
    ]

    ch_t_manual = [
        resolve_to_atlas_channel_id(ch, config_info)
        for ch in ch_t_manual_raw
    ]

    # ------------------------------------------------------------------
    # 2. Sanity-check manual pairs
    # ------------------------------------------------------------------
    checked_manual = check_reflected_transmitted_pairs(
        r_channels=ch_r_manual,
        t_channels=ch_t_manual,
    )

    for warning in checked_manual["warnings"]:
        CustomWarning(warning)

    ch_r_manual = [
        pair["numerator_channel"]
        for pair in checked_manual["pairs"]
    ]

    ch_t_manual = [
        pair["denominator_channel"]
        for pair in checked_manual["pairs"]
    ]

    manual_pairs = set(zip(ch_r_manual, ch_t_manual))

    # ------------------------------------------------------------------
    # 3. Find all automatic pairs from atlas_channel_id
    # ------------------------------------------------------------------
    ch_r_auto, ch_t_auto = find_reflected_transmitted_pairs(atlas_channel_id)

    # ------------------------------------------------------------------
    # 4. Remove automatic pairs already provided manually
    # ------------------------------------------------------------------
    ch_r_auto_new = []
    ch_t_auto_new = []

    for r, t in zip(ch_r_auto, ch_t_auto):
        if (r, t) in manual_pairs:
            continue

        ch_r_auto_new.append(r)
        ch_t_auto_new.append(t)

    # ------------------------------------------------------------------
    # 5. Combine manual + automatic pairs
    # ------------------------------------------------------------------
    ch_r_all = ch_r_manual + ch_r_auto_new
    ch_t_all = ch_t_manual + ch_t_auto_new

    n_pairs = len(ch_r_all)

    config_info["ch_r"] = ch_r_all
    config_info["ch_t"] = ch_t_all

    # ------------------------------------------------------------------
    # 6. Fill the rest of the POL_CAL keys
    # ------------------------------------------------------------------
    for key in POL_CAL_KEYS:

        if key in ["ch_r", "ch_t", "pol_cal_pairs"]:
            continue

        meta = SCHEMA[key]

        if meta["is_list"]:
            config_info[key] = extend_or_default_pol_cal_list(
                existing=config_info.get(key, []),
                default=meta["default"],
                n_target=n_pairs,
            )
        else:
            value = config_info.get(key)

            if value in [None, ""]:
                config_info[key] = str(meta["default"])
            else:
                config_info[key] = str(value)

    config_info = dict(sorted(config_info.items()))

    return config_info

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
        print(i, ch, "bins =", config_info["bins"][i])
        for key in necesary_channel_parameters:
            print(i,ch,key)
            
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
    pol_cal_info = pd.DataFrame()
    water_vapour_info = pd.DataFrame()
    temperature_info = pd.DataFrame()

    # Always initialize metadata groups
    metadata["pol_cal_info"] = {}
    metadata["water_vapour_info"] = {}
    metadata["temperature_info"] = {}

    for key in config_info:

        if key in SYSTEM_KEYS:
            system_info[key] = config_info[key]

        elif key in CHANNEL_KEYS:
            channel_info.loc[:, key] = [v for v in config_info[key]]

        elif key in POL_CAL_KEYS:
            pol_cal_info.loc[:, key] = [v for v in config_info[key]]

        elif key in WV_KEYS:
            water_vapour_info.loc[:, key] = [v for v in config_info[key]]

        elif key in TEMP_KEYS:
            temperature_info.loc[:, key] = [v for v in config_info[key]]

    system_parameters = list(system_info.keys())
    system_values = np.array(
        list(system_info.values()),
        dtype=object,
    )

    has_pol_cal_info = pol_cal_info.index.size > 0
    has_water_vapour_info = water_vapour_info.index.size > 0
    has_temperature_info = temperature_info.index.size > 0

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

        if has_pol_cal_info:
            metadata["pol_cal_info"][key] = xr.DataArray(
                pol_cal_info.T.values,
                dims=["parameters", "pair"],
                coords={
                    "parameters": pol_cal_info.columns.values,
                    "pair": pol_cal_info.index.values,
                },
            )

        if has_water_vapour_info:
            metadata["water_vapour_info"][key] = xr.DataArray(
                water_vapour_info.T.values,
                dims=["parameters", "pair"],
                coords={
                    "parameters": water_vapour_info.columns.values,
                    "pair": water_vapour_info.index.values,
                },
            )

        if has_temperature_info:
            metadata["temperature_info"][key] = xr.DataArray(
                temperature_info.T.values,
                dims=["parameters", "pair"],
                coords={
                    "parameters": temperature_info.columns.values,
                    "pair": temperature_info.index.values,
                },
            )

    return metadata