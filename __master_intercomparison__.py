#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

import os, warnings, glob, sys
from helper_functions import processing_chain
from helper_functions import parse_intercomparison_args
from helper_functions.parse_intercomparison_args import call_parser as parse_mst
from helper_functions.parse_intercomparison_args import check_parser as check_mst
from visualizer.readers.parse_cmp_args import call_parser as parse_cmp
from helper_functions import read_intercomparison_config

warnings.filterwarnings('ignore')

mst_args = parse_mst()
mst_args = check_mst(mst_args)

cmp_cfg = read_intercomparison_config.config(mst_args['settings_file'])

# Reset argument list
sys.argv = [sys.argv[0]]   

# Default intercomparison arguments
cmp_args = parse_cmp()

cmp_args = parse_intercomparison_args.substitute(org = cmp_args, rpl = cmp_cfg.cmp)

# Call Intercomparison sequence
processing_chain.intercomparison(cmp_args = cmp_args, 
                                 input_files = mst_args['input_files'],
                                 cmp_out = mst_args['output_folder'])