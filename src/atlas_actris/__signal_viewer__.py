#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

# from __master__ import main as atlas_master
from utils.filtering import filter_channels
from utils.find_radiosonde import find_radiosonde
from utils.get_scc_config import export_scc_config
from utils.parse_init_file import parse_call_atlas_ini
from visualizer.signal_viewer import generate_line_plots
from utils.parse_config_file import parse_atlas_config_file
from utils.cleaners import ask_clean_cache, ask_clean_viewer
from utils.parse_settings_file import parse_atlas_settings_file

from utils.parse_caller_args import call_parser
from processor.pipeline import Context, Processor
from readers.reader_utils import special_path_rules
from readers.read_radiosondes import load_radiosonde
from readers.read_raw_lidar_files_cache import infer_format, flexible_reader

from utils.cookbook import (
    recipes,
    checkin_stages, 
    checkout_stages,
    run_linear_recipe
    )

from processor.modify import (
    load_metadata, 
    remove_unrecognised_channels,
    get_atlas_channel_id, 
    load_pol_cal_defaults,
    bring_to_correct_type, 
    special_config_checks, 
    store_updated_metadata,
    )


# Get the input .ini file path of the ATLAS caller
cmd_args = call_parser()

# Parse the initialization file
caller_info = parse_call_atlas_ini(filepath = cmd_args['ini_file'])

# Export the config file from SCC HOI
scc_info = export_scc_config(
    scc_configuration_id = caller_info['scc_configuration_id'], 
    atlas_configuration_file = caller_info['atlas_configuration_file'], 
    export_hoi_cfg = caller_info['export_hoi_cfg'],
    output_folder = caller_info['output_folder'],
    scc_compatible_format = caller_info['scc_compatible_format']
    )

# Parse the configuration file
config_info = parse_atlas_config_file(caller_info['atlas_configuration_file'])

# Parse the settings file
settings_info = parse_atlas_settings_file(caller_info['atlas_settings_file'])

# Inferring the raw file format
caller_info["raw_file_format"] = infer_format(caller_info, station_id = config_info["station_id"], debug=True)

# Modifying absolute paths for scc and polly_xt readers that have 2 measurements embeeded in one file
caller_info = special_path_rules(caller_info)
     
# Reading the files
profiles, metadata = flexible_reader(caller_info)

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

# Find radiosonde - Download if file does not exist
metadata = find_radiosonde(caller_info, metadata)

# Load radiosonde files
meteo, metadata = load_radiosonde(caller_info, metadata)

# Filter out channels
profiles, meatadata = filter_channels(caller_info, profiles, metadata)     

# Load all data needed for processing in a data class
ctx = Context(
    processing_info = {'caller_info':caller_info, 'settings_info':settings_info},
    starting_dataset = {'profile' : profiles, 'meteo' : meteo} | metadata,
    )

# Initialize the Processor class
processor = Processor(ctx)

for key, recipe in recipes.items():
    run_linear_recipe(
        processor, 
        recipe = recipe, 
        initial_input = checkin_stages[key],
        checkout_id = checkout_stages[key],
        )

for view_stage in caller_info['view_mean_signal_stages']:
    generate_line_plots(processor, stage = view_stage, db = 'profile_mean')

for view_stage in caller_info['view_signal_stages']:
    generate_line_plots(processor, stage = view_stage, db = 'profile')

# Commandline promt to clean cache or not
ask_clean_cache(caller_info)

# Commandline promt to clean signal_viewer output or not
ask_clean_viewer(caller_info)
