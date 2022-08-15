#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:48:11 2022

@author: nick
"""

import warnings, os, sys
from readers.parse_config import parse_config
from readers import read_files
from lidar_processing.short_prepro import rayleigh

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the command line argument information
args = parse_config()

lidar_info, channel_info, time_info, time_info_d, sig, sig_d, shots, shots_d = \
    read_files.short_reader(args['input_file'])
    
sig_pack = rayleigh(sig_raw = sig, shots = shots, lidar_info = lidar_info, \
                    channel_info = channel_info, time_info = time_info,\
                        external_info = args)