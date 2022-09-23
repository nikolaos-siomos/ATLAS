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

isday = True

lidar_1 = 'POLIS'
lidar_2 = 'ALPHA'

channels_1 = ['xppr0355', 'xppt0532', 'xpar0355', 'xpat0532']
channels_2 = ['nppr0355', 'nppr0532', 'npar0355', 'npar0532']

# Force reprocessing
reprocess = {'scc_converter' : False,
             'atlas_preprocessor' : False,
             'atlas_visualizer' : True}

parent_folder_1 = '/mnt/DATA/Big_data/Databases/POLIS/220922/01'
parent_folder_2 = '/mnt/DATA/Big_data/Databases/ALPHA/220922'

output_folder_1 = os.path.join(parent_folder_1, 'out')
output_folder_2 = os.path.join(parent_folder_2, 'out')

output_folder_com = '/mnt/DATA/Big_data/Databases/Intercomparisons/POLIS_ALPHA/220922'

# Create the main output folder
if os.path.exists(parent_folder_1):
    os.makedirs(output_folder_1, exist_ok = True)
    
if os.path.exists(parent_folder_2):
    os.makedirs(output_folder_2, exist_ok = True)
    
# Converter arguments
cnv_args = parse_cnv()

cnv_args['trim_overflows'] = 2
cnv_args['files_per_sector'] = 3
cnv_args['radiosonde_folder'] = '/mnt/DATA/Big_data/Databases/POLIS/radiosondes'
cnv_args['rsonde_skip_header'] = 4
cnv_args['rsonde_delimiter'] = 'S'
cnv_args['rsonde_column_index'] = [2, 1, 3, 5]
cnv_args['rsonde_geodata'] = [44.50, 26.13, 91.0]
cnv_args['rsonde_column_units'] = ['m', 'hPa', 'C', 'percent']

cnv_args_1 = cnv_args.copy()
cnv_args_2 = cnv_args.copy()

cnv_args_1['parent_folder'] = parent_folder_1
cnv_args_1['output_folder'] = os.path.join(output_folder_1, 'scc_converter')
cnv_args_1['config_file'] = '/mnt/DATA/Big_data/Databases/POLIS/configs/polis_bucharest_all_channels.ini'
cnv_args_1['measurement_identifier'] = 'm'

cnv_args_2['parent_folder'] = parent_folder_2
cnv_args_2['output_folder'] = os.path.join(output_folder_2, 'scc_converter')
cnv_args_2['config_file'] = '/mnt/DATA/Big_data/Databases/ALPHA/configs/alpha_bucharest_all_channels.ini'
cnv_args_2['measurement_identifier'] = 'RM'

# ATLAS preprocessor arguments
prs_args = parse_prs()

prs_args['vertical_limit'] = 20.
if isday == True:
    prs_args['exclude_scattering_type'] = ['v', 'f']

prs_args_1 = prs_args.copy()
prs_args_2 = prs_args.copy()


prs_args_1['channels'] = channels_1
prs_args_1['output_folder'] = os.path.join(output_folder_1, 'atlas_preprocessor')

prs_args_2['channels'] = channels_2
prs_args_2['output_folder'] = os.path.join(output_folder_2, 'atlas_preprocessor')

# Common plotting arguments
vis_inp = prs_args['output_folder']
use_distance = True
dpi = 150

# Rayleigh fit arguments
ray_args = parse_ray()

ray_args['x_lims'] = [0., prs_args['vertical_limit']]
ray_args['x_tick'] = 2.

ray_args['reference_height'] = 7.
ray_args['half_reference_window'] = 500.

ray_args['smooth'] = True 
ray_args['smoothing_range'] = [1., prs_args['vertical_limit']]
ray_args['half_window'] = [5., 1000.]

ray_args['output_folder'] = output_folder_com

# Intercomparison arguments
cmp_args = parse_cmp()

cmp_args['channels_1'] = channels_1
cmp_args['channels_2'] = channels_2

cmp_args['x_lims'] = [0.2, 8]
cmp_args['x_tick'] = 1.

cmp_args['normalization_height'] = 3.
cmp_args['half_normalization_window'] = 200.

cmp_args['smooth'] = True 
cmp_args['smoothing_range'] = [1., prs_args['vertical_limit']]
cmp_args['half_window'] = [5., 1000.]

cmp_args['output_folder'] = output_folder_com

# Call lidar_1 scc_converter sequence
processing_chain.scc_converter(cnv_args = cnv_args_1, 
                               reprocess = reprocess['scc_converter'])

# Call lidar_1 Rayleigh Fit sequence - only preprocessing
processing_chain.QA_test(input_folder = cnv_args_1['output_folder'], 
                         prs_args = prs_args_1, 
                         vis_args = ray_args, 
                         test_type = 'ray', 
                         reprocess_prs = reprocess['atlas_preprocessor'], 
                         reprocess_vis = reprocess['atlas_visualizer'],
                         quicklook = False,
                         skip = True)

# Call lidar_2 scc_converter sequence
processing_chain.scc_converter(cnv_args = cnv_args_2, 
                               reprocess = reprocess['scc_converter'])

# Call lidar_2 Rayleigh Fit sequence - only preprocessing
processing_chain.QA_test(input_folder = cnv_args_2['output_folder'], 
                         prs_args = prs_args_2, 
                         vis_args = ray_args, 
                         test_type = 'ray', 
                         reprocess_prs = reprocess['atlas_preprocessor'], 
                         reprocess_vis = reprocess['atlas_visualizer'],
                         quicklook = False,
                         skip = True)

# Call Intercomparison sequence
processing_chain.intercomparison(input_folder_1 = prs_args_1['output_folder'],
                                 input_folder_2 = prs_args_2['output_folder'],
                                 vis_args = cmp_args, 
                                 reprocess = reprocess['atlas_visualizer'],
                                 skip = False)