#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

# from __master__ import main as atlas_master
from utils.__get_scc_config__ import export_scc_config
from utils.__parse_init_file__ import parse_call_atlas_ini
from utils.__parse_config_file__ import parse_atlas_config_file
from utils.__parse_settings_file__ import parse_atlas_settings_file
from readers.read_raw_lidar_files import infer_format, flexible_reader
from helper_functions.parse_caller_args import call_parser
from readers.reader_utils import special_path_rules
from trimming.modify import load_metadata, remove_unrecognised_channels, \
    get_atlas_channel_id, load_pol_cal_defaults, \
        bring_to_correct_type, special_config_checks, store_updated_metadata
from helper_functions.printouts import endpoint
from processor.pipeline import Context, Processor

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
    
# Load all data needed for processing in a data class
ctx = Context(
    processing_info = {'caller_info':caller_info, 'loading_map':loading_map},
    settings_info = settings_info,
    starting_dataset = {'profile':profiles} | metadata,
    )

# Initialize the Processor class
processor = Processor(ctx)

# Calculate the profiles of ranges and heights per bin
processor.run(output_id = 'ranges_and_heights', input_id = 'init', stage_name = "height_and_range_calculation")

# Slice and exclude parts of the measurements
processor.run(output_id = 'sliced', input_id = 'ranges_and_heights', stage_name = "slice_and_exclude")

# Remove "incomplete" profiles with shots less than a certain percentage of the max shots
processor.run(output_id = 'shots_screened', input_id = 'sliced', stage_name = "screen_low_shots")

# Handle overflow data (detect, interpolate, ignore)
processor.run(output_id = 'overflows_checked', input_id = 'shots_screened', stage_name = "handling_overflows")

# Checkpoint: End of filtering data 
processor.copy(output_id = 'screening_complete', input_id = 'overflows_checked')

# Detect saturation and clipping
processor.run(output_id = 'saturation_detected', input_id = 'screening_complete', stage_name = "check_saturation")

# Unit conversion - raw counts to MHz for the photon channels
processor.run(output_id = 'photon_units_converted', input_id = 'screening_complete', stage_name = "photon_units_conversion")

# Dead time correction for photon channels
processor.run(output_id = 'dead_time_corrected', input_id = 'photon_units_converted', stage_name = "dead_time_correction")

# Solar background calculation - full temporal resolution
processor.run(output_id = 'background_raw_res', input_id = 'dead_time_corrected', stage_name = "background_calculation")

# Averaging profiles - single averages for all datasets and special handling for ray
processor.run(output_id = 'averaged', input_id = 'dead_time_corrected', stage_name = "averaging_by_time")
processor.run(output_id = 'averaged_low_res', input_id = 'dead_time_corrected', stage_name = "averaging_by_time_low_res")
processor.run(output_id = 'averaged_high_res', input_id = 'dead_time_corrected', stage_name = "averaging_by_time_high_res")

# Solar background calculation
processor.run(output_id = 'background', input_id = 'averaged', stage_name = "background_calculation")
processor.run(output_id = 'background_low_res', input_id = 'averaged_low_res', stage_name = "background_calculation")
processor.run(output_id = 'background_high_res', input_id = 'averaged_high_res', stage_name = "background_calculation")

# Background correction
processor.run(output_id = 'background_corrected', input_id = 'background', stage_name = "background_correction")
processor.run(output_id = 'background_corrected_low_res', input_id = 'background_low_res', stage_name = "background_correction")
processor.run(output_id = 'background_corrected_high_res', input_id = 'background_high_res', stage_name = "background_correction")

# Range correction
processor.run(output_id = 'range_corrected', input_id = 'background_corrected', stage_name = "range_correction")
processor.run(output_id = 'range_corrected_low_res', input_id = 'background_corrected_low_res', stage_name = "range_correction")
processor.run(output_id = 'range_corrected_high_res', input_id = 'background_corrected_high_res', stage_name = "range_correction")

# Trim signals and ranges/heights vertically
processor.run(output_id = 'vert_trimmed', input_id = 'range_corrected', stage_name = "trim_vertically")
processor.run(output_id = 'vert_trimmed_low_res', input_id = 'range_corrected_low_res', stage_name = "trim_vertically")
processor.run(output_id = 'vert_trimmed_high_res', input_id = 'range_corrected_high_res', stage_name = "trim_vertically")

# Dark correction
processor.run(output_id = 'dark_corrected', input_id = 'vert_trimmed', stage_name = "dark_correction")
processor.run(output_id = 'dark_corrected_low_res', input_id = 'vert_trimmed_low_res', stage_name = "dark_correction")
processor.run(output_id = 'dark_corrected_high_res', input_id = 'vert_trimmed_high_res', stage_name = "dark_correction")

# Checkpoint: End of preprocessing
processor.copy(output_id = 'preprocessing_complete', input_id = 'dark_corrected')
processor.copy(output_id = 'preprocessing_complete_low_res', input_id = 'dark_corrected_low_res')
processor.copy(output_id = 'preprocessing_complete_high_res', input_id = 'dark_corrected_high_res')

# Identify gluing region
processor.run(output_id = 'gluing_region', input_id = 'preprocessing_complete', stage_name = 'gluing_region')
processor.run(output_id = 'gluing_region_low_res', input_id = 'preprocessing_complete_low_res', stage_name = 'gluing_region')
processor.run(output_id = 'gluing_region_high_res', input_id = 'preprocessing_complete_high_res', stage_name = 'gluing_region')

# Glue
processor.run(output_id = 'glued', input_id = 'gluing_region', stage_name = 'gluing')
processor.run(output_id = 'glued_low_res', input_id = 'gluing_region_low_res', stage_name = 'gluing')
processor.run(output_id = 'glued_high_res', input_id = 'gluing_region_high_res', stage_name = 'gluing')

raise Exception

# 3.a) Reading meteorological profiles from ECMWF (cloudnet)
meteo, metadata = select_radiosonde_ecmwf(
    input_folder = ecmwf_rsonde_folder, 
    signal_times = sig.time.values, 
    metadata = metadata
    )

# 3.b) Reading meteorological profiles from Wyoming if cloudnet data are not available
if metadata['file_validity_flag'] == 'empty':
    meteo, metadata = select_radiosonde_wyoming(
        input_folder = wyoming_rsonde_folder, 
        signal_times = sig.time.values, 
        metadata = metadata
        )

    if metadata['file_validity_flag'] in ['empty', 'unusable']:
        endpoint(1)


import numpy as np
window = 10
input_stage = "raw"
 
for key in profile_db[input_stage].keys():

    zero_bin = metadata["channel_info"][key].loc["zero_bin"]
    
    max_zero_bin = zero_bin.max().values
    
    if max_zero_bin < -240:

        sig = profile_db[input_stage][key][:,:,:abs(max_zero_bin)]#.sel({"bins": slice(abs(max_zero_bin)-1)})
        
        roll = sig.rolling(bins = window)
        sig_wv = roll.mean()
        sig_wv_err = roll.std() / np.sqrt(window)
        err = sig_wv_err.min(dim="bins")
        
        wct_err = 4. * np.sqrt(2) * err
        
        wct = sig_wv.shift(bins=window) - sig_wv[:,:,-window]#.sel({"bins": slice(-window)})
        
        
        
        


    
    # sig = signal.dead_time_correction(sig = sig.copy(), 
    #                                   dead_time = dead_time, 
    #                                   dead_time_cor_type = dead_time_cor_type)
    
        
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
