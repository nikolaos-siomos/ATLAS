#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

# from __master__ import main as atlas_master
import sys
from utils.__get_scc_config__ import export_scc_config
from utils.__parse_init_file__ import parse_call_atlas_ini
from utils.__parse_config_file__ import parse_atlas_config_file
from utils.__parse_settings_file__ import parse_atlas_settings_file
from readers.read_raw_lidar_files import infer_format, flexible_reader
from helper_functions.parse_caller_args import call_parser
from readers.reader_utils import special_path_rules
from trimming.modify import load_metadata, remove_unrecognised_channels, \
    slice_and_exclude, get_atlas_channel_id, load_pol_cal_defaults, \
        bring_to_correct_type, special_config_checks, store_updated_metadata, \
            screen_low_shots, initialize_debug_dictionaries
from trimming.handle_overflows import check_for_overflows
from processor.signal import detect_saturation    

# from helper_functions.caller_utils import autodetect_paths, prepare_master_args, export_report

# Get the input .ini file path of the ATLAS caller
cmd_args = call_parser()

# Parse the initialization file
caller_info = parse_call_atlas_ini(filepath = cmd_args['ini_file'])

# Export the config file from SCC HOI
scc_info = export_scc_config(scc_configuration_id = caller_info['scc_configuration_id'], 
                             atlas_configuration_file = caller_info['atlas_configuration_file'], 
                             export_hoi_cfg = caller_info['export_hoi_cfg'],
                             output_folder = caller_info['output_folder'])

# Parse the configuration file
config_info = parse_atlas_config_file(caller_info['atlas_configuration_file'])

# Parse the settings file
settings_info = parse_atlas_settings_file(caller_info['atlas_settings_file'])

# Inferring the raw file format
caller_info["raw_file_format"] = infer_format(caller_info, station_id = config_info["station_id"], debug=True)

# Modifying absolute paths for scc and polly_xt readers that have 2 measurements embeeded in one file
caller_info = special_path_rules(caller_info)
      
# Reading the files
profiles, metadata, loading_map = flexible_reader(caller_info)

# Transfer metadata from the raw file header to config_info
config_info = load_metadata(config_info, metadata)

# Special config checks related to parameters not provided with the config file nor the raw file header
config_info = special_config_checks(config_info)

# Remove all profiles corresponding to channels not provided in the config file
remove_unrecognised_channels(caller_info, config_info, profiles, metadata)

# Create the altas_channel_id now that all channel metadata are loaded and channels are screened
config_info = get_atlas_channel_id(config_info, profiles, metadata)

# Filling pol cal parameters of config_info with defaults if empty
config_info = load_pol_cal_defaults(config_info)

# Change the type of the config_info entries from str to the correct one according to the config file parser SCHEMA
config_info = bring_to_correct_type(config_info)

# Replace all system_info, channel_info, and pol_cal_info metadata with the common updated ones
metadata = store_updated_metadata(config_info, metadata)
    
# Initialize the profile_db
profile_masks, profile_db = initialize_debug_dictionaries(profiles)

# Slice and exclude parts of the measurements
profiles, metadata, profile_masks = slice_and_exclude(caller_info, profiles, metadata, profile_masks)

# Remove "incomplete" profiles with shots less than a certain percentage of the max shots
profiles, profile_masks = screen_low_shots(profiles, metadata, profile_masks)

# Handle overflow data (detect, interpolate, ignore)
profiles, metadata, profile_masks = check_for_overflows(caller_info, profiles, metadata, profile_masks)

# Detect saturation and clipping
profile_masks = detect_saturation(caller_info, profiles, metadata, profile_masks)



    # sig = signal.dead_time_correction(sig = sig.copy(), 
    #                                   dead_time = dead_time, 
    #                                   dead_time_cor_type = dead_time_cor_type)
    
    # sig, time_info = \
    #     signal.average_by_time(sig = sig.copy(),
    #                            time_info = time_info,
    #                            timescale = -1,
    #                            start_time = 'Raw_Data_Start_Time',
    #                            stop_time = 'Raw_Data_Stop_Time')
        
    # bgr = signal.background_calculation(sig = sig.copy(), 
    #                                     lower_bin = bg_low,
    #                                     upper_bin = bg_high)
    
    # sig = signal.trigger_correction(sig = sig.copy(), 
    #                            daq_trigger_offset = trd_bins)
    
    # sig = signal.trim_vertically(sig = sig.copy(), 
    #                              ground_alt = ground_alt,
    #                              zenith_angle = zenith_angle, 
    #                              alt_lim = 1E3 * alt_lim,
    #                              resol = resol)
    
    # ranges = signal.range_calculation(bins = sig.copy().bins.values, 
    #                                   resol = resol)
    
    # heights = signal.height_calculation(bins = sig.copy().bins.values, 
    #                                     resol = resol,
    #                                     zenith_angle = zenith_angle)
    
    # sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    # sig = signal.dark_correction(sig = sig.copy(), 
    #                              drk = sig_drk.copy())
    # sig = signal.range_correction(sig = sig, ranges = ranges)
