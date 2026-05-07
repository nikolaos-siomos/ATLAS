#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

import sys
from .helper_functions.get_scc_config import parse_args, export_scc_config
from pprint import pprint

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]
        
    # Get the input .ini file path of the ATLAS caller
    cmd_args = parse_args()
    
    # Parse the fields from the initialization file
    data = export_scc_config(
        scc_configuration_id = cmd_args['scc_configuration_id'],
        atlas_configuration_file = cmd_args['atlas_configuration_file'],
        hoi_output_folder =cmd_args['hoi_output_folder'],
        verbose =cmd_args['verbose'],
        )

if __name__ == "__main__":
    
    main()