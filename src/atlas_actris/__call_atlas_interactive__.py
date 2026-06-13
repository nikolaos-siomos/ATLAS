#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

# from __master__ import main as atlas_master
from utils.cleaners import ask_clean_cache
from utils.caller_utils import export_report
from utils.find_radiosonde import find_radiosonde
from utils.__get_scc_config__ import export_scc_config
from utils.__parse_init_file__ import parse_call_atlas_ini
from utils.__parse_config_file__ import parse_atlas_config_file
from utils.__parse_settings_file__ import parse_atlas_settings_file

from visualizer.__quicklook_lazy__ import generate_quicklooks
from visualizer.__rayleigh_fit__ import generate_rayleigh_fit
from visualizer.__ring_telecover__ import generate_ring_telecover
from visualizer.__quicklook_lazy_vldr__ import generate_vldr_quicklooks
from visualizer.__quadrant_telecover__ import generate_quadrant_telecover
from visualizer.__polarization_calibration__ import generate_polarization_calibration

from utils.parse_caller_args import call_parser
from processor.pipeline import Context, Processor
from readers.reader_utils import special_path_rules
from readers.read_radiosondes import load_radiosonde
from readers.read_raw_lidar_files_cache import infer_format, flexible_reader

from utils.cookbook import (
    run_linear_recipe,
    screening_recipe, 
    preprocessing_recipe, 
    pol_cal_recipe
    )

from trimming.modify import (
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
    output_folder = caller_info['output_folder']
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

# # Expand profiles dictionary using the loading_map aliases
# profiles = expand_with_loading_map(profiles, caller_info)

# # Expand each entry of metadata dictionary using the loading_map aliases
# metadata = expand_nested_with_loading_map(metadata, caller_info)

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

# Load all data needed for processing in a data class
ctx = Context(
    processing_info = {'caller_info':caller_info, 'settings_info':settings_info},
    starting_dataset = {'profile' : profiles, 'meteo' : meteo} | metadata,
    )

# Initialize the Processor class
processor = Processor(ctx)

# Apply screening recipe
run_linear_recipe(
    processor, 
    recipe = screening_recipe, 
    initial_input = "init",
    checkout_id = "screening_complete",
    )

# Apply preprocessing recipe
run_linear_recipe(
    processor, 
    recipe = preprocessing_recipe, 
    initial_input = "screening_complete",
    checkout_id = "preprocessing_complete",
    )

run_linear_recipe(
    processor, 
    recipe = pol_cal_recipe, 
    initial_input = "preprocessing_complete",
    checkout_id = "pol_cal_complete",
    )

# Package measurements for quicklooks
processor.package(output_id = 'preprocessing_complete_qck', input_id = 'preprocessing_complete')

pol_cal__metadata = generate_polarization_calibration(
    data_pack = processor.export_test_from_stage("pol_cal_complete"),
    caller_info = processor.processing_info["caller_info"],
    settings_info = settings_info,
)

# Rayleigh fit test
rayleigh_fit__metadata = generate_rayleigh_fit(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info
    )

# Quadrant telecover test
quadrant_telecover__metadata = generate_quadrant_telecover(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info['tlc_qua']
    )

# Ring telecover test
ring_telecover__metadata = generate_ring_telecover(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info['tlc_rin']
    )

# Quicklooks
generate_quicklooks(
    data_pack = processor.export_test_from_stage('preprocessing_complete_qck'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info['qck'],
    )

generate_vldr_quicklooks(
    data_pack = processor.export_test_from_stage('pol_cal_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info['qck'],
    )

# Create report
export_report(caller_info)

# Commandline promt to clean cache or not
ask_clean_cache(caller_info)
