#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

# from __master__ import main as atlas_master
from utils.find_radiosonde import find_radiosonde
from utils.__get_scc_config__ import export_scc_config
from utils.__parse_init_file__ import parse_call_atlas_ini
from utils.__parse_config_file__ import parse_atlas_config_file
from utils.__parse_settings_file__ import parse_atlas_settings_file

from visualizer.__quicklook__ import generate_quicklooks
from visualizer.__rayleigh_fit__ import generate_rayleigh_fit
from visualizer.__quadrant_telecover__ import generate_quadrant_telecover
from visualizer.__ring_telecover__ import generate_ring_telecover

from utils.cookbook import run_linear_recipe
from utils.parse_caller_args import call_parser
from processor.pipeline import Context, Processor
from readers.reader_utils import special_path_rules
from readers.read_radiosondes import load_radiosonde
from readers.read_raw_lidar_files import infer_format, flexible_reader

from trimming.modify import (
    load_metadata, 
    remove_unrecognised_channels,
    get_atlas_channel_id, 
    load_pol_cal_defaults,
    bring_to_correct_type, 
    special_config_checks, 
    store_updated_metadata
    )

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

# Find radiosonde - Download if file does not exist
metadata = find_radiosonde(caller_info, metadata)

# Load radiosonde files
meteo, metadata = load_radiosonde(caller_info, metadata)

# Load all data needed for processing in a data class
ctx = Context(
    processing_info = {'caller_info':caller_info, 'loading_map':loading_map},
    settings_info = settings_info,
    starting_dataset = {'profile' : profiles, 'meteo' : meteo} | metadata,
    )

# Initialize the Processor class
processor = Processor(ctx)

screening_recipe = [
    ("ranges_and_heights", "height_and_range_calculation"),
    ("sliced", "slice_and_exclude"),
    ("shots_screened", "screen_low_shots"),
    ("overflows_checked", "handling_overflows"),
    ("saturation_detected", "check_saturation"),
]

preprocessing_recipe = [
    ("photon_units_converted", "photon_units_conversion"),
    ("dead_time_corrected", "dead_time_correction"),
    ("background_calculated", "background_calculation"),
    ("background_corrected", "background_correction"),
    ("range_corrected", "range_correction"),
    ("dark_corrected", "dark_correction"),
    ("vert_trimmed", "trim_vertically"),
    ("gluing_region_found", "gluing_region"),
    ("glued", "gluing"),
    ("noise_calculated", "signal_noise_calculation"),
    ("molecular_calculated", "molecular_calculations"),
]

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

# Averaging profiles - single averages for all datasets and special handling for ray
processor.run(output_id = 'averaged', input_id = 'preprocessing_complete', stage_name = "averaging_by_time")
processor.run(output_id = 'averaged_lr', input_id = 'preprocessing_complete', stage_name = "averaging_by_time_low_res")
processor.run(output_id = 'averaged_hr', input_id = 'preprocessing_complete', stage_name = "averaging_by_time_high_res")

# Package the measurements for quicklooks
processor.package(output_id = 'averaged_hr_qck', input_id = 'averaged_hr')

# Rayleigh fit test
rayleigh_fit__metadata = generate_rayleigh_fit(
    data_pack = processor.export_test_from_stage('averaged'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = processor.settings_info
    )


# Quadrant telecover test
quadrant_telecover__metadata = generate_quadrant_telecover(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = processor.settings_info['tlc_qua']
    )

# Ring telecover test
ring_telecover__metadata = generate_ring_telecover(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = processor.settings_info['tlc_rin']
    )


# Quicklooks
generate_quicklooks(
    data_pack = processor.export_test_from_stage('averaged_hr_qck'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = processor.settings_info['qck']
    )
