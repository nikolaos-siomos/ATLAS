#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

# from __master__ import main as atlas_master
import os
from utils.cleaners import ask_clean_cache
from utils.filtering import filter_channels
from utils.caller_utils import export_report
from utils.find_radiosonde import find_radiosonde
from utils.get_scc_config import export_scc_config
from utils.parse_init_file import parse_call_atlas_ini
from utils.parse_config_file import parse_atlas_config_file
from utils.parse_settings_file import parse_atlas_settings_file

from visualizer.qa_test_dark import generate_dark
from visualizer.qa_test_quicklook import generate_quicklooks
from visualizer.qa_test_background import generate_background
from visualizer.qa_test_rayleigh_fit import generate_rayleigh_fit
from visualizer.qa_test_ring_telecover import generate_ring_telecover
from visualizer.qa_test_quicklook_vldr import generate_vldr_quicklooks
from visualizer.qa_test_quadrant_telecover import generate_quadrant_telecover
from visualizer.qa_test_polarization_calibration import generate_polarization_calibration

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

from utils.export_processing_stage import (
    export_processor_stages,
    delete_all_exported_stages,
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

# Package measurements for quicklooks
processor.package(
    output_id = 'preprocessing_complete_qck', 
    input_id = 'preprocessing_complete'
    )

pol_cal_metadata = generate_polarization_calibration(
    data_pack = processor.export_test_from_stage("pol_cal_complete"),
    caller_info = processor.processing_info["caller_info"],
    settings_info = settings_info,
)

# Rayleigh fit test
rayleigh_fit_metadata = generate_rayleigh_fit(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info
    )

# Quadrant telecover test
quadrant_telecover_metadata = generate_quadrant_telecover(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info['tlc']
    )

# Ring telecover test
ring_telecover_metadata = generate_ring_telecover(
    data_pack = processor.export_test_from_stage('preprocessing_complete'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info['tlc_rin']
    )


#Dark test
dark_metadata = generate_dark(
    data_pack = processor.export_test_from_stage("dark_preprocessing_complete"),
    caller_info = processor.processing_info["caller_info"],
    settings_info = settings_info['drk'],
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
    settings_info = settings_info['qck_vldr'],
    )

# Background
generate_background(
    data_pack = processor.export_test_from_stage('preprocessing_complete_qck'),
    caller_info = processor.processing_info['caller_info'],
    settings_info = settings_info['bgd'],
    )

# Create report
export_report(caller_info)

# Commandline promt to clean cache or not
ask_clean_cache(caller_info)

# Delete all exported stages
delete_exported_answer = os.environ.get("ATLAS_DELETE_EXPORTED_ANSWER")

if delete_exported_answer is None:
    delete_all_exported_stages(
        output_folder=caller_info['output_folder'],
    )
elif delete_exported_answer.strip().lower() in ["y", "yes"]:
    delete_all_exported_stages(
        output_folder=caller_info['output_folder'],
        ask=False,
    )
else:
    print("-- Skipping deletion of all exported processing stages")

# Export latest stage
export_stage_answer = os.environ.get("ATLAS_EXPORT_STAGE_ANSWER")

if export_stage_answer is None:
    export_processor_stages(
        processor=processor,
        stage_names = caller_info['export_stages'],
        output_folder = caller_info['output_folder'],
        overwrite=True,
        ask=True,
        default_answer=False,
        print_estimated_size=True,
    )
elif export_stage_answer.strip().lower() in ["y", "yes"]:
    export_processor_stages(
        processor=processor,
        stage_names = caller_info['export_stages'],
        output_folder = caller_info['output_folder'],
        overwrite=True,
        ask=False,
        default_answer=False,
        print_estimated_size=True,
    )
else:
    print("-- Skipping export of processing stages")
