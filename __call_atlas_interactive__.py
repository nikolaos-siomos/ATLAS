#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

from __master__ import main as atlas_master
from helper_functions.parse_caller_args import call_parser, parse_ini_file, check_parser_args
from helper_functions.caller_utils import autodetect_paths, prepare_master_args, export_report

# Get the input .ini file path of the ATLAS caller
cmd_args = call_parser()

# Parse the fields from the initialization file
parser_args = parse_ini_file(path = cmd_args['ini_file'])

# The the input file fields
check_parser_args(parser_args)

# Infer the SCC configuration ID, SC lidar ID and data_idetnifier from the folder structure
parser_args = autodetect_paths(parser_args)

# Prepare the ATLAS arguments
mst_args = prepare_master_args(parser_args)

# Call ATLAS
atlas_master(mst_args)

# Export to an html draft QA report file
export_report(parser_args)

  
