#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 18:16:24 2025

@author: nikos
"""

import csv
from pathlib import Path
import configparser
from typing import Sequence, Hashable, Tuple, List, Dict, Any, Optional
from datetime import datetime
import numpy as np
import argparse
import requests
from utils.printouts import print_header
from pprint import pprint
from utils.error_classes import ConfigError, CustomWarning

sections_map = {
    "System": {
        "station_id": ("atlas_system", "station_id"),
        "version_id": ("atlas_system", "version_id"),
        "lidar_name": ("atlas_system", "lidar_name"),
        "configuration_name": ("atlas_system", "configuration_name"),
        "configuration_id": ("atlas_system", "configuration_id"),
        "station_altitude": ("atlas_system", "station_altitude"),
        "station_latitude": ("atlas_system", "station_latitude"),
        "station_longitude": ("atlas_system", "station_longitude"),
        "zenith_angle": ("atlas_system", "zenith_angle"),
    },
    "Channels": {
        "scc_channel_id": ("atlas_channels", "scc_channel_id"),
        "recorder_channel_id": ("atlas_channels", "recorder_channel_id"),
        "telescope_type": ("atlas_channels", "telescope_type"),
        "channel_type": ("atlas_channels", "channel_type"),
        "channel_subtype": ("atlas_channels", "channel_subtype"),
        "acquisition_mode": ("atlas_channels", "acquisition_mode"),
        "zero_bin": ("atlas_channels", "zero_bin"),
        "dead_time": ("atlas_channels", "dead_time"),
        "detected_wavelength": ("atlas_channels", "detected_wavelength"),
        "emitted_wavelength": ("atlas_channels", "emitted_wavelength"),
        "channel_bandwidth": ("atlas_channels", "channel_bandwidth"),
        "G": ("atlas_channels", "G"),
        "H": ("atlas_channels", "H"),
        "range_resolution": ("atlas_channels", "range_resolution"),
        "laser_repetition_rate": ("atlas_channels", "laser_repetition_rate"),
        # "laser_id": ("atlas_channels", "laser_id")
    },
}

internal_map = {
    "atlas_system": {
        "station_id": ("HoiSystems", "owner_id_id"),
        "version_id": ("HoiSystems", "lidar_version_id"),
        "lidar_name": ("HoiSystems", "name"),
        "configuration_name": ("HoiSystems", "configuration"),
        "configuration_id": ("HoiSystems", "id"),
        "station_altitude": ("HoiSystems", "height_asl"),
        "station_latitude": ("HoiSystems", "Latitude"),
        "station_longitude": ("HoiSystems", "Longitude"),
        "zenith_angle": ("HoiSystems", "lidar_pointing_angle"),
    },
    "atlas_channels": {
        "scc_channel_id": ("HoiChannels", "id"),
        "dead_time": ("HoiChannels", "dead_time"),
        "detected_wavelength": ("HoiChannels", "if_center"),
        "emitted_wavelength": ("HoiChannels", "emission_wavelength"),
        "channel_bandwidth": ("HoiChannels", "if_fwhm"),
        "range_resolution": ("HoiChannels", "raw_range_resolution"),
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process optional lists of IDs"
    )

    parser.add_argument(
        "-c",
        "--scc_configuration_id",
        type=str,
        nargs="?",          # one or more strings if provided
        default=None,       # None if not given
        help="The SCC configuration ID"
    )
    parser.add_argument(
        "-o",
        "--atlas_configuration_file",
        type=str,
        nargs="?",          # one or more strings if provided
        default=None,       # None if not given
        help="The path to the atlas configuration file"
    )
    parser.add_argument(
        "-e",
        "--export_hoi_cfg",
        type=str,
        nargs="?",          # one or more strings if provided
        default="1",       # None if not given
        help="The path to the exported ini file"
    )
    parser.add_argument(
        "-f",
        "--output_folder",
        type=str,
        nargs="?",          # one or more strings if provided
        default="./scc_hoi",       # None if not given
        help="The path to the exported ini file"
    )
    parser.add_argument(
        "--scc_compatible_format",
        "--scc_format",
        action="store_true",
        default=False,
        help=(
            "If provided, use the SCC channel IDs directly as recorder_channel_id "
            "values and skip the interactive recorder_channel_id prompt."
        ),
    )
        
    args = parser.parse_args()
    scc_configuration_id = args.scc_configuration_id
    atlas_configuration_file = args.atlas_configuration_file
    export_hoi_cfg = args.export_hoi_cfg
    output_folder = args.output_folder
    
    validate_string(scc_configuration_id, var_name = "scc_configuration_id")
    
    if atlas_configuration_file == None:
        atlas_configuration_file = f"./configurations/config_file_{scc_configuration_id}.ini"
        args.atlas_configuration_file = atlas_configuration_file
        Path(atlas_configuration_file).parent.mkdir(exist_ok=True)
        CustomWarning(f"The scc_config_file_path was not provided. The following path was selected by default: {atlas_configuration_file}")
    
    validate_string(atlas_configuration_file, var_name = "atlas_configuration_file")
    validate_string(export_hoi_cfg, var_name = "export_hoi_cfg", allowed_values = ["0", "1", "2"])
    validate_string(output_folder, var_name = "output_folder")

    return args


def validate_string(
    value,
    *,
    allowed_values: Optional[Sequence[str]] = None,
    var_name: str = "value"
) -> str:
    """
    Validate that `value` is a string (or numpy.str_).
    - Raises TypeError if not a string.
    - Raises ValueError if None or not in allowed_values.
    - Returns the validated string.
    """
    import numpy as np

    # Reject None
    if value is None:
        raise ValueError(f"{var_name} must not be None")

    # Type check
    if not isinstance(value, (str, np.str_)):
        raise TypeError(f"{var_name} must be a string, got {type(value).__name__}")

    # Allowed values check
    if allowed_values is not None and value not in allowed_values:
        raise ValueError(
            f"{var_name} must be one of {allowed_values}, got '{value}'"
        )

    return str(value)


def validate_bool(value, *, var_name: str = "value") -> bool:
    """Validate and normalize bool-like values."""

    if isinstance(value, bool):
        return value

    if isinstance(value, (str, np.str_)):
        value_str = str(value).strip().lower()
        if value_str in ["true", "1", "yes", "y", "on"]:
            return True
        if value_str in ["false", "0", "no", "n", "off", ""]:
            return False

    if isinstance(value, (int, np.integer)) and value in [0, 1]:
        return bool(value)

    raise ValueError(
        f"{var_name} must be a boolean or one of true/false, 1/0, yes/no; "
        f"got {value!r}"
    )


def download_with_auth(url: str, folder: str, username: str, password: str, filename: str = None) -> str:
    folder_path = Path(folder)
    folder_path.mkdir(parents=True, exist_ok=True)

    if filename is None:
        filename = url.split("/")[-1] or "downloaded_file"
    file_path = folder_path / filename

    response = requests.get(url, auth=(username, password), stream=True)
    if response.status_code == 401:
        raise PermissionError("Unauthorized: check your username/password")
    response.raise_for_status()

    with open(file_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    return str(file_path)

    
def _get_scc_file(output_folder: str, scc_configuration_id: str) -> str:

    output_folder = Path(output_folder) / "scc_hoi"
    output_folder.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    filename = f"scc_config_{scc_configuration_id}_{timestamp}.csv"
    
    url = f"https://scc.imaa.cnr.it/admin/database/hoisystems/export/{scc_configuration_id}"
    
    scc_file_path = download_with_auth(url = url, folder = output_folder, filename = filename, username = "scc_user", password = "sccforever!")

    return(scc_file_path)

def map_ids(values: List[Hashable]) -> Tuple[List[int], Dict[Hashable, int]]:
    """
    Assign an incremental number to each unique ID in the order
    of first appearance. All duplicates (even scattered) share the same number.

    Parameters
    ----------
    values : List[Hashable]
        The list of IDs (can contain duplicates).

    Returns
    -------
    result : List[int]
        A list of integers (same length as `values`), where each element
        is the group number assigned to that ID.
    mapping : Dict[Hashable, int]
        A dictionary mapping each unique ID to its assigned group number.
    """
    mapping: Dict[Hashable, int] = {}
    next_num: int = 1
    result: List[int] = []

    for val in values:
        if val not in mapping:
            mapping[val] = next_num
            next_num += 1
        result.append(str(mapping[val]))

    return result, mapping

def label_groups_edge_fixed(
    a: Sequence[Hashable],
    b: Sequence[float],
) -> Tuple[List[str], Dict[Any, str]]:
    """
    Label each unique value in `a` based on its unique corresponding b value.
    Smallest b -> 'n', largest b -> 'f', and middle labels depend on the count:
      1: n
      2: n, f
      3: n, o, f
      4: n, o, e, f
      5: n, o, p, e, f
      6: n, o, p, d, e, f
    Assumes each unique value in `a` maps to exactly one unique b.
    """
    if len(a) != len(b):
        raise ValueError("a and b must have the same length")

    # Build mapping a_val -> unique b, verifying consistency
    group_b: Dict[Any, float] = {}
    for av, bv in zip(a, b):
        if av in group_b and group_b[av] != bv:
            raise ValueError(f"Inconsistent b for group {av}: {group_b[av]} vs {bv}")
        group_b[av] = bv

    k = len(group_b)
    if k == 0:
        return [], {}
    if k > 6:
        raise ValueError(f"At most 6 unique entries supported; got {k}")

    # Label sets by count
    label_by_count = {
        1: ['n'],
        2: ['n','f'],
        3: ['n','o','f'],
        4: ['n','o','e','f'],
        5: ['n','o','p','e','f'],
        6: ['n','o','p','d','e','f'],
    }
    labels = label_by_count[k]

    # Order groups by b
    ordered_groups = sorted(group_b.keys(), key=lambda v: group_b[v])

    # Assign labels in that order
    group_to_label = {grp: lab for grp, lab in zip(ordered_groups, labels)}

    # Apply to full list
    labeled = [group_to_label[av] for av in a]
    return labeled, group_to_label

def _copy_atlas_keys(data: Dict[str, Dict[str, List]], atlas_section: str) -> Dict[str, Dict[str, List]]:
    
    data[atlas_section] = {}
    for key in internal_map[atlas_section].keys():    
        orig_section = internal_map[atlas_section][key][0]
        orig_key = internal_map[atlas_section][key][1]
        data[atlas_section][key] = data[orig_section][orig_key]

        # Keep the same existence checks as for all other mapped parameters.
        # Only after the value has been copied successfully, replace empty
        # zenith_angle entries by the ATLAS default vertical-pointing value.
        if atlas_section == "atlas_system" and key == "zenith_angle":
            values = data[atlas_section][key]
            if isinstance(values, list):
                data[atlas_section][key] = [
                    "0." if v is None or str(v).strip() == "" else v
                    for v in values
                ]
            elif values is None or str(values).strip() == "":
                data[atlas_section][key] = "0."
    
    return(data)

def _get_zero_bin(data:Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    
    n_channels = len(data["HoiChannels"]["id"])
    first_signal_rangebin = data["HoiChannels"].get("first_signal_rangebin")
    trigger_delay = data["HoiChannels"].get("trigger_delay")
    raw_range_resolution = data["HoiChannels"]["raw_range_resolution"]
    if (first_signal_rangebin, trigger_delay) != (None, None):
            zero_bin = []
            for i in range(n_channels):
                if int(first_signal_rangebin[i]) >= 1:
                    zero_bin.append(f"-{first_signal_rangebin[i]}")
                else:
                    raw_time_resolution = 2. * float(raw_range_resolution[i]) * 1e9 / (3 * 1e8)
                    trigger_delay_bins = int(float(trigger_delay[i]) / float(raw_time_resolution))
                    zero_bin.append(f"{trigger_delay_bins}")
            data["atlas_channels"]["zero_bin"] = zero_bin
    else:
        data["atlas_channels"]["zero_bin"] = n_channels * ['']

    return(data)

def _get_channel_types(data:Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    
    n_channels = len(data["HoiChannels"]["id"])
    scat_mechanism = data["HoiChannels"].get("_scat_mechanism_id_id")
    emitted_wavelength = data["HoiChannels"]["emission_wavelength"]
    detected_wavelength = data["HoiChannels"]["if_center"]
    channel_id = data['HoiChannels']["id"]
    scat_map = {
        "elT":        {"channel_type":"t",       "channel_subtype":["x", "t", "r"], "telescope_type": ["x","y","z"]},
        "elTnr":      {"channel_type":"t",       "channel_subtype":["x", "t", "r"], "telescope_type": ["n","o","p"]},
        "elTfr":      {"channel_type":"t",       "channel_subtype":["x", "t", "r"], "telescope_type": ["d","e","f"]},
        "elTunr":     {"channel_type":"t",       "channel_subtype":["x", "t", "r"], "telescope_type": ["n","o","p"]},
        
        # "elPR":       {"channel_type":["p","c"], "channel_subtype":"r", "telescope_type": ["x","y","z"]},
        # "elPRnr":     {"channel_type":["p","c"], "channel_subtype":"r", "telescope_type": ["n","o","p"]},
        # "elPRfr":     {"channel_type":["p","c"], "channel_subtype":"r", "telescope_type": ["d","e","f"]},
        
        # "elPT":       {"channel_type":["p","c"], "channel_subtype":["t", "r"], "telescope_type": ["x","y","z"]},
        # "elPTnr":     {"channel_type":["p","c"], "channel_subtype":["t", "r"], "telescope_type": ["n","o","p"]},
        # "elPTfr":     {"channel_type":["p","c"], "channel_subtype":["t", "r"], "telescope_type": ["d","e","f"]},
        
        "elPP":       {"channel_type":"p",       "channel_subtype":["t", "r"], "telescope_type": ["x","y","z"]},
        "elPPnr":     {"channel_type":"p",       "channel_subtype":["t", "r"], "telescope_type": ["n","o","p"]},
        "elPPfr":     {"channel_type":"p",       "channel_subtype":["t", "r"], "telescope_type": ["d","e","f"]},
        
        "elCP":       {"channel_type":"c",       "channel_subtype":["t", "r"], "telescope_type": ["x","y","z"]},
        "elCPnr":     {"channel_type":"c",       "channel_subtype":["t", "r"], "telescope_type": ["n","o","p"]},
        "elCPfr":     {"channel_type":"c",       "channel_subtype":["t", "r"], "telescope_type": ["d","e","f"]},
        
        "vrRN2":      {"channel_type":["v","r"], "channel_subtype":"n", "telescope_type": ["x","y","z"]},
        "vrRN2nr":    {"channel_type":["v","r"], "channel_subtype":"n", "telescope_type": ["n","o","p"]},
        "vrRN2fr":    {"channel_type":["v","r"], "channel_subtype":"n", "telescope_type": ["d","e","f"]},
        
        "vrRH2O":     {"channel_type":"v", "channel_subtype":"w", "telescope_type": ["x","y","z"]},
        "vrRH2Onr":   {"channel_type":"v", "channel_subtype":"w", "telescope_type": ["n","o","p"]},
        "vrRH2Ofr":   {"channel_type":"v", "channel_subtype":"w", "telescope_type": ["d","e","f"]},
        
        "pRRlow":     {"channel_type":"r",       "channel_subtype":"l", "telescope_type": ["x","y","z"]},
        "pRRlownr":   {"channel_type":"r",       "channel_subtype":"l", "telescope_type": ["n","o","p"]},
        "pRRlowfr":   {"channel_type":"r",       "channel_subtype":"l", "telescope_type": ["d","e","f"]},
        
        "pRRhigh":    {"channel_type":"r",       "channel_subtype":"h", "telescope_type": ["x","y","z"]},
        "pRRhighnr":  {"channel_type":"r",       "channel_subtype":"h", "telescope_type": ["n","o","p"]},
        "pRRhighfr":  {"channel_type":"r",       "channel_subtype":"h", "telescope_type": ["d","e","f"]},
        }
    
    if scat_mechanism == None:
        channel_type = n_channels * ['']
        channel_subtype = n_channels * ['']
    elif "" in scat_mechanism:
        channel_type = n_channels * ['']
        channel_subtype = n_channels * ['']
    else:
        channel_type = []
        channel_subtype = []
        for i in range(n_channels):
            ch_type = scat_map[scat_mechanism[i]]["channel_type"]
            ch_stype = scat_map[scat_mechanism[i]]["channel_subtype"]
            if not isinstance(ch_type, list):
                channel_type.append(ch_type)
            else:
                if ch_type == ["v","r"]:
                    if np.abs(float(detected_wavelength[i]) - float(emitted_wavelength[i])) > 10:
                        channel_type.append("v")
                    else:
                        channel_type.append("r")
                        
            if not isinstance(ch_stype, list):
                channel_subtype.append(ch_stype)
            else:
                if ch_stype == ["t", "r"] or  ch_stype == ["x", "t", "r"]:
                    if data.get("PRODUCSTCHANNELS") not in [None,{}]:
                        channel_id_list = data["PRODUCSTCHANNELS"]["channelID"]
                        channel_stype_list = data["PRODUCSTCHANNELS"]["signalTypeID"]
                        ch_stype_ids = set([ch_tp for ch_id, ch_tp in zip(channel_id_list, channel_stype_list) if ch_id == channel_id[i]])
                        if np.any([ch_id in ['6', '10', '11'] for ch_id in ch_stype_ids]):
                            channel_subtype.append("r")
                        elif np.any([ch_id in ['7', '12', '13'] for ch_id in ch_stype_ids]):
                            channel_subtype.append("t")
                        else:
                            channel_subtype.append("x") 
                    else:
                        channel_subtype.append("_") 
    
    data["atlas_channels"]["channel_type"] = channel_type
    data["atlas_channels"]["channel_subtype"] = channel_subtype
    
    return(data)

def _get_telescope_type(data:Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    
    n_channels = len(data["HoiChannels"]["id"])
    telescope_ids = data["HoiChannels"].get("_telescope_id_id")
    telescope_unique_ids = data["HoiTelescopes"].get("id")
    full_overlap_height_m = data["HoiTelescopes"].get("full_overlap_height_m")
    if telescope_ids != None:
        if len(telescope_unique_ids) == 1:
            telescope_type = n_channels * ['x']
        else:
            if full_overlap_height_m == None:
                telescope_type = label_list(telescope_ids)
            elif "" in full_overlap_height_m:
                telescope_type, _ = label_list(telescope_ids)
            else:
                fovl = []
                for t_id in telescope_ids:
                    t_index = telescope_unique_ids.index(t_id)
                    fovl.append(float(full_overlap_height_m[t_index]))
                telescope_type, _ = label_groups_edge_fixed(telescope_ids, fovl) 
    
    data["atlas_channels"]["telescope_type"] = telescope_type
    
    return(data)

def _get_acquisition_mode(data:Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:

    detection_mode = data["HoiChannels"].get("_detection_mode_id_id")
    if detection_mode != None:
        acquisition_mode = np.array(detection_mode)
        acquisition_mode[acquisition_mode == 'an'] = 'a'
        acquisition_mode[acquisition_mode == 'pc'] = 'p'
        acquisition_mode = [str(aq) for aq in acquisition_mode]
    
    data["atlas_channels"]["acquisition_mode"] = acquisition_mode
    
    return(data)

def _get_laser_parameters(data:Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    laser_id, _ = map_ids(data['HoiChannels']["_laser_id_id"])
    # data['atlas_channels']["laser_id"] = laser_id
    
    laser_ids = data['HoiChannels']["_laser_id_id"]
    laser_repetition_rate = []
    for l_id in laser_ids:
        l_index = data["HoiLaser"]['id'].index(l_id)
        laser_repetition_rate.append(data["HoiLaser"]['repetition_rate'][l_index])
    
    data['atlas_channels']["laser_repetition_rate"] = laser_repetition_rate
    
    return(data)    

def _get_GH(data:Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    channel_id = data["HoiChannels"]["id"]
    if data["PolarizationCrosstalkParameter"] != {}:
        pol_ch_ids = data["PolarizationCrosstalkParameter"]["channel_id"]
        G_list = data["PolarizationCrosstalkParameter"]["g"]
        H_list = data["PolarizationCrosstalkParameter"]["h"]
        if pol_ch_ids != []:
            G = []
            H = []
            for ch in channel_id:
                if ch in pol_ch_ids:
                    p_index = len(pol_ch_ids) - 1 - pol_ch_ids[::-1].index(ch)
                    if G_list[p_index] == "0E-8":
                        G.append("0.")
                    else:
                        G.append(G_list[p_index])
                        
                    if H_list[p_index] == "0E-8":
                        H.append("0.")
                    else:
                        H.append(H_list[p_index])
                        
                else:
                    G.append("_")
                    H.append("_")
            data["atlas_channels"]["G"] = G
            data["atlas_channels"]["H"] = H
        else:
            data["atlas_channels"]["G"] = "_"
            data["atlas_channels"]["H"] = "_"
        
            
    return(data)

def _get_atlas_channel_id(data: Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    
    detected_wavelength = [str(int(round(float(dwl),0))).zfill(4) for dwl in data["HoiChannels"]["if_center"]]
    telescope_type = data["atlas_channels"]["telescope_type"]
    channel_type = data["atlas_channels"]["channel_type"]
    acquisition_mode = data["atlas_channels"]["acquisition_mode"]
    channel_subtype = data["atlas_channels"]["channel_subtype"]
    
    atlas_channel_id = [ai + bi + ci + di + ei for ai, bi, ci, di, ei in zip(detected_wavelength, telescope_type, channel_type, acquisition_mode, channel_subtype)]
    
    data["atlas_channels"]["atlas_channel_id"] = atlas_channel_id

    return(data)

def _atlas_id_sort_key(atlas_id: str) -> Tuple[str, str]:
    """Sort ATLAS channel IDs by telescope type first, then alphabetically."""
    atlas_id = str(atlas_id)
    telescope_type = atlas_id[4] if len(atlas_id) >= 5 else ""
    return telescope_type, atlas_id

def _sort_atlas_channels_by_atlas_id(data: Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    """
    Reorder all per-channel atlas_channels entries using atlas_channel_id.

    The primary key is the fifth atlas_id character, i.e. telescope type.
    The secondary key is the complete atlas_id in alphabetical order.
    Non-channel metadata entries are left untouched.
    """
    atlas_channels = data.get("atlas_channels")
    if not atlas_channels or "atlas_channel_id" not in atlas_channels:
        return data

    atlas_channel_id = atlas_channels["atlas_channel_id"]
    if not isinstance(atlas_channel_id, list):
        return data

    order = sorted(range(len(atlas_channel_id)), key=lambda i: _atlas_id_sort_key(atlas_channel_id[i]))

    for key, values in list(atlas_channels.items()):
        if isinstance(values, list) and len(values) == len(atlas_channel_id):
            atlas_channels[key] = [values[i] for i in order]

    return data

# def _get_recorder_channel_id(data: Dict[str, Dict[str, List[str]]], 
#                               scc_channel_id: Optional[List[str]], 
#                               recorder_channel_id: Optional[List[str]]) -> Dict[str, Dict[str, List[str]]]:
    
#     if scc_channel_id == None and recorder_channel_id != None:
#         if len(recorder_channel_id) != len(data["HoiChannels"]["id"]):
#             raise ConfigError("recorder_channel_id must have the same number of elements as the HOI channels. Subsets are onlyu supported if the scc_channel_id is also provided. Please use a subset a recorder_channel_id for each existing IDs:\n{data['HoiChannels']['id']} ")
    
#     if scc_channel_id != None:
#         if len(scc_channel_id) > len(data["HoiChannels"]["id"]):
#             raise ConfigError(f"scc_channel_id size cannot be larger than the existing scc channel IDs in the HOI. Please use a subset of the existing IDs:\n{data['HoiChannels']['id']}")
#         data["atlas_channels"] = reorder_data_by_external(data["HoiChannels"], scc_channel_id = scc_channel_id)

#     if recorder_channel_id != None:
#         data["atlas_channels"]["recorder_channel_id"] = recorder_channel_id
        
#     return(data)

def _get_dead_time(data: Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:

    dead_time = ["_" if v == "" else v for v in data["HoiChannels"]["dead_time"]]
    data["atlas_channels"]["dead_time"] = dead_time
    
    return(data)


def _exclude_undifined_scc_channels(data: Dict[str, Dict[str, List[str]]]) -> Dict[str, Dict[str, List[str]]]:
    atlas_channel_id = data["atlas_channels"]["atlas_channel_id"]
    scc_channel_id = data["atlas_channels"]["scc_channel_id"]
    mask = ["_" in ch for ch in atlas_channel_id]
    
    if any(mask):
        unidentified_scc_channels = [val for val, reject in zip(scc_channel_id, mask) if reject]
        unidentified_atlas_channels = [val for val, reject in zip(atlas_channel_id, mask) if reject]
        CustomWarning("The following SCC channels were not properly defined and will not be included in the configuration file:")
        for i in range(len(unidentified_scc_channels)):
            print(f"    -- {unidentified_scc_channels[i]}: {unidentified_atlas_channels[i]}")
    for key in data["atlas_channels"].keys():
        values = data["atlas_channels"][key]
        data["atlas_channels"][key] = [val for val, reject in zip(values, mask) if not reject]

    data["atlas_system"] = {k: v for k, v in data["atlas_system"].items() if not all(x == "_" for x in v)}
    data["atlas_channels"] = {k: v for k, v in data["atlas_channels"].items() if not all(x == "_" for x in v)}
    
    return(data)

def _get_recorder_channel_id(
    data: Dict[str, Dict[str, List[str]]],
    scc_compatible_format: bool = False,
) -> Dict[str, Dict[str, List[str]]]:
    data = _sort_atlas_channels_by_atlas_id(data)

    atlas_channel_id = data["atlas_channels"]["atlas_channel_id"]
    scc_channel_id = data["atlas_channels"]["scc_channel_id"]

    if scc_compatible_format:
        recorder_channel_id = [str(ch) for ch in scc_channel_id]
        data["atlas_channels"]["recorder_channel_id"] = recorder_channel_id

        print(
            "-- SCC-compatible format requested. "
            "Using scc_channel_id values as recorder_channel_id values."
        )
        print("Final channel list:")
        for scc_id, atlas_id, rec_id in zip(scc_channel_id, atlas_channel_id, recorder_channel_id):
            print(f"    -- {scc_id} [{atlas_id}] --> {rec_id}")

        if len(recorder_channel_id) != len(set(recorder_channel_id)):
            raise ConfigError(
                f"Duplicates detected in the SCC-derived recorder_channel_id values: "
                f"{recorder_channel_id}"
            )

        return data

    recorder_channel_id = []
    print("--Input: Provide recorder_channel_id values for each scc_channel_d: (leave empty and press Enter to exclude a channel): ")
    for scc_id, atlas_id in zip(scc_channel_id, atlas_channel_id):
        print(f"    -- {scc_id} [{atlas_id}]:", end = "", flush = True)
        rec_id = input()
        recorder_channel_id.append(rec_id)  # could be "" if the user just presses Enter

    if all(v == "" for v in recorder_channel_id):
        CustomWarning("No recorder_channel_id was provided. Please fill in the recorder_channel_id manually in the config exported file")
    else:
        data["atlas_channels"]["recorder_channel_id"] = recorder_channel_id

        mask = [ch != "" for ch in recorder_channel_id]

        for key in data["atlas_channels"].keys():
            values = data["atlas_channels"][key]
            data["atlas_channels"][key] = [val for val, keep in zip(values, mask) if keep]

        print("Final channel list:")
        for scc_id, atlas_id, rec_id in zip(scc_channel_id, atlas_channel_id, recorder_channel_id):
            if rec_id != "":
                print(f"    -- {scc_id} [{atlas_id}] --> {rec_id}")
            else:
                print(f"    -- {scc_id} [{atlas_id}] --> Excluded!")

        rec_ch_id = recorder_channel_id.copy()
        if "" in rec_ch_id:
            rec_ch_id.remove("")
            if len(rec_ch_id) != len(set(rec_ch_id)):
                raise ConfigError(f"Duplicates detected in the provided recorder_channel_id values: {recorder_channel_id}")

    return(data)
    
def ensure_file(path_str: str):
    """
    Ensure that the given file path exists.
    - Creates parent directories if they don't exist.
    - Creates an empty file if it doesn't exist.
    """
    path = Path(path_str)

    # make sure parent folders exist
    path.parent.mkdir(parents=True, exist_ok=True)

    # make sure file exists
    if not path.exists():
        path.touch()
        print(f"Configuration file does not exist, creating empty file: {path}")
        
def export_selected_to_ini(
    data: Dict[str, Dict[str, List]],
    filename: str,
    sections_map: Dict[str, Dict[str, Tuple[str, str]]],
    joiner: str = ",",
) -> None:
    """
    Write one or more *new* sections to an INI file, with per-key renaming.

    sections_map format (mapping style only):
        {
          "NewSectionName": {
              "target_key_1": ("SourceSection", "source_key_1"),
              "target_key_2": ("SourceSection", "source_key_2"),
              ...
          },
          ...
        }
    """
    data = _sort_atlas_channels_by_atlas_id(data)

    cfg = configparser.ConfigParser()

    # Preserve option/key case when writing the ATLAS config file.
    # ConfigParser lowercases option names by default, which would turn
    # schema keys like G and H into g and h.
    cfg.optionxform = str

    for new_section, mapping in sections_map.items():
        if not cfg.has_section(new_section):
            cfg.add_section(new_section)

        for target_key, (src_sec, src_key) in mapping.items():
            # Skip quietly if missing
            if src_sec not in data or src_key not in data[src_sec]:
                continue
            values = data[src_sec][src_key]
            val_str = joiner.join(str(v) for v in values) if isinstance(values, list) else str(values)
            cfg.set(new_section, target_key, val_str)

    ensure_file(filename)
    
    with open(Path(filename), "w", encoding="utf-8") as f:
        cfg.write(f)


def read_sectioned_csv(path: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Read an SCC HOI sectioned CSV export into:
        {section_name: {column_name: [values...]}}

    The SCC export is a CSV file split into sections:
        SectionName
        header1,header2,...
        row1...

        NextSection
        ...

    This reader deliberately uses csv.reader on the file object instead of
    parsing line-by-line. That makes it robust to quoted fields containing
    embedded newlines and commas, for example a multi-line HoiSystems
    description. It is also tolerant to rows with missing trailing columns,
    which some SCC exports omit instead of writing explicit empty fields.
    """

    data: Dict[str, Dict[str, List[str]]] = {}

    def _is_blank(row: List[str]) -> bool:
        return len(row) == 0 or all(str(cell).strip() == "" for cell in row)

    def _normalize_row_length(
        row: List[str],
        headers: List[str],
        section_name: str,
    ) -> List[str]:
        """Pad short rows and merge surplus fields into the last column."""
        if len(row) == len(headers):
            return row

        if len(row) < len(headers):
            return row + [""] * (len(headers) - len(row))

        # This should be rare when csv.reader is used correctly, but keeping
        # the fallback makes the parser tolerant to malformed rows with extra
        # separators. Preserve the extra content instead of silently dropping it.
        if len(headers) == 0:
            raise ValueError(
                f"Row detected before a header in section {section_name}: {row}"
            )

        return row[:len(headers) - 1] + [",".join(row[len(headers) - 1:])]

    with open(Path(path), "r", encoding="utf-8", errors="replace", newline="") as f:
        rows = list(csv.reader(f))

    i, n = 0, len(rows)

    while i < n:
        # Skip blank separator rows between sections.
        while i < n and _is_blank(rows[i]):
            i += 1

        if i >= n:
            break

        # Section name. SCC section names are stored as one-cell rows.
        section_row = rows[i]
        section_name = str(section_row[0]).strip() if len(section_row) > 0 else ""
        i += 1

        if section_name == "":
            continue

        if "CHECKSUM" in section_name.upper():
            break

        # Empty section: section name followed by EOF or a blank separator.
        if i >= n or _is_blank(rows[i]):
            data[section_name] = {}
            if i < n and _is_blank(rows[i]):
                i += 1
            continue

        headers = [str(h).strip() for h in rows[i]]
        i += 1

        # Remove a possible UTF-8 BOM from the first header if present.
        if headers:
            headers[0] = headers[0].lstrip("\ufeff")

        section_dict: Dict[str, List[str]] = {h: [] for h in headers}

        # Collect rows until the blank separator before the next section.
        while i < n and not _is_blank(rows[i]):
            row = _normalize_row_length(rows[i], headers, section_name)

            for h, val in zip(headers, row):
                section_dict[h].append(val)

            i += 1

        data[section_name] = section_dict

        # Skip the blank separator if any.
        if i < n and _is_blank(rows[i]):
            i += 1

    return data

def label_list(values: List[str]) -> Tuple[List[str], Dict[str, List[str]]]:
    labels = ['x', 'y', 'z']
    mapping = {}
    result = []
    for v in values:
        if v not in mapping:
            if len(mapping) >= len(labels):
                raise ValueError("More than 3 unique entries in list")
            mapping[v] = labels[len(mapping)]
        result.append(mapping[v])
    return result, mapping

def reorder_data_by_external(
    data: Dict[str, Any],
    scc_channel_id: List[str],
    key_id: str = "id",
) -> Dict[str, Any]:
    """
    Reorder all one-to-one list fields in `data` to follow the order in `c`,
    where `c` is a subset/permutation of data[key_a].

    - Raises ValueError if `c` contains elements not in data[key_a].
    - Returns a new dict with the same keys; list fields aligned with `a`
      are filtered and reordered; other fields are copied unchanged.

    Example:
        data = {"a": ["A","B","C","D"], "b":[10,20,30,40], "meta":"x"}
        c = ["C","A","D"]
        -> {"a":["C","A","D"], "b":[30,10,40], "meta":"x"}
    """
    a = data[key_id]
    if not isinstance(a, list):
        raise TypeError(f"`data[{key_id!r}]` must be a list")

    # Build index lookup from a-value -> index
    pos = {val: i for i, val in enumerate(a)}

    # Validate c and build the new index order
    scc_channel_id_list = list(scc_channel_id)
    missing = [x for x in scc_channel_id_list if x not in pos]
    if missing:
        raise ValueError(f"Values in `scc_channel_id` not found in data[{key_id!r}]: {missing}")

    new_idx = [pos[x] for x in scc_channel_id_list]

    # Rebuild output dict: reorder all list fields of the same length as `a`
    out: Dict[str, Any] = {}
    for k, v in data.items():
        if isinstance(v, list) and len(v) == len(a):
            out[k] = [v[i] for i in new_idx]   # filtered & reordered
        else:
            out[k] = v                         # copy unchanged (e.g., metadata)

    return out

def export_scc_config(scc_configuration_id: str, 
                      atlas_configuration_file: str, 
                      export_hoi_cfg: str, 
                      output_folder: str,
                      scc_compatible_format: bool = False,
                      debug: bool = False) -> Dict[str, Dict[str, List[str]]]:
    
    if scc_configuration_id != None:
        validate_string(atlas_configuration_file, var_name = "atlas_configuration_file")
        validate_string(scc_configuration_id,     var_name = "scc_configuration_id")
        validate_string(output_folder,            var_name = "output_folder")
        validate_string(export_hoi_cfg,           var_name = "export_hoi_cfg", allowed_values = ["0", "1", "2"])
        scc_compatible_format = validate_bool(
            scc_compatible_format,
            var_name = "scc_compatible_format",
        )
    
        if export_hoi_cfg == "1" or (export_hoi_cfg == "2" and not Path(atlas_configuration_file).exists()):
            
            print_header(f"Exporting SCC HOI - Config ID:{scc_configuration_id}\nto: {atlas_configuration_file}")
            
            scc_file_path = _get_scc_file(output_folder = output_folder, 
                                          scc_configuration_id = scc_configuration_id)
        
            data = read_sectioned_csv(scc_file_path)
            
            data = _copy_atlas_keys(data, atlas_section = "atlas_system")
            data = _copy_atlas_keys(data, atlas_section = "atlas_channels")
                
            data = _get_zero_bin(data)
            
            data = _get_channel_types(data)
        
            data = _get_telescope_type(data)
        
            data = _get_acquisition_mode(data)
        
            data = _get_laser_parameters(data)
            
            data = _get_dead_time(data)
            
            data = _get_atlas_channel_id(data)
        
            data = _exclude_undifined_scc_channels(data)
        
            data = _get_GH(data)    
        
            data = _get_recorder_channel_id(
                data,
                scc_compatible_format = scc_compatible_format,
            )    
                
            if debug:
                pprint(data)
                
            export_selected_to_ini(data, atlas_configuration_file, sections_map)
            
            return(data)

