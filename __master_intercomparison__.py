#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob
from helper_functions import processing_chain

from scc_converter.readers.parse_args import call_parser as parse_cnv
from processor.readers.parse_args import call_parser as parse_prs
from visualizer.readers.parse_cmp_args import call_parser as parse_cmp
from visualizer.readers.parse_ray_args import call_parser as parse_ray
# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

channels_1 = ['xppr0355', 'xppt0532', 'xpar0355', 'xpat0532']
# channels_2 = ['nppr0355', 'nppr0532', 'npar0355', 'npar0532']
channels_2 = ['xppr0355', 'xppr0532', 'xpar0355', 'xpar0532']

input_folder_1 = '/mnt/DATA/Big_data/Databases/CARS/Systems/POLIS/220924/01/out/atlas_preprocessor'
# input_folder_2 = '/mnt/DATA/Big_data/Databases/CARS/Systems/ALPHA/220924/01/out/atlas_preprocessor'
input_folder_2 = '/mnt/DATA/Big_data/Databases/CARS/Systems/EMORAL/220924/01/out/atlas_preprocessor'

output_folder_com = '/mnt/DATA/Big_data/Databases/CARS/Intercomparisons/POLIS_ALPHA_EMORAL/220924/'

# Intercomparison arguments
cmp_args = parse_cmp()

cmp_args['channels_1'] = channels_1
cmp_args['channels_2'] = channels_2

cmp_args['x_lims'] = [0., 25]
cmp_args['x_tick'] = 2.

cmp_args['normalization_height'] = 7.
cmp_args['half_normalization_window'] = 500.

cmp_args['smooth'] = True 
cmp_args['smoothing_range'] = [1., 20.]
cmp_args['half_window'] = [5., 1000.]

cmp_args['output_folder'] = output_folder_com

# Call Intercomparison sequence
processing_chain.intercomparison(input_folder_1 = input_folder_1,
                                 input_folder_2 = input_folder_2,
                                 vis_args = cmp_args, 
                                 reprocess = True,
                                 skip = False)