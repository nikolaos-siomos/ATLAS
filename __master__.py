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
from visualizer.readers.parse_qck_args import call_parser as parse_qck
from visualizer.readers.parse_ray_args import call_parser as parse_ray
from visualizer.readers.parse_tlc_args import call_parser as parse_tlc
from visualizer.readers.parse_pcl_args import call_parser as parse_pcl
# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

isday = False

# Force reprocessing
process = {'rayleigh_fit' : True,
           'telecover' : True,
           'polarization_calibration' : True}

reprocess = {'scc_converter' : False,
             'atlas_preprocessor' : False,
             'atlas_visualizer' : True}

# Skip module
skip_visualizer = {'rayleigh_fit' : False,
                   'telecover' : False,
                   'polarization_calibration' : True,
                   'quicklook' : False}

# Control quicklooks
quicklook = {'rayleigh_fit' : True,
             'telecover' : False,
             'polarization_calibration' : False}

parent_folder = '/mnt/DATA/Big_data/Databases/CARS/Systems/RALI/220924/03'
output_folder = os.path.join(parent_folder, 'out')

# Create the main output folder
if os.path.exists(parent_folder):
    os.makedirs(output_folder, exist_ok = True)
    
# Converter arguments
cnv_args = parse_cnv()

cnv_args['parent_folder'] = parent_folder
cnv_args['output_folder'] = os.path.join(output_folder, 'scc_converter')
cnv_args['telecover_sectors_folder'] = os.path.join(parent_folder, 'tlc')
cnv_args['config_file'] = '/mnt/DATA/Big_data/Databases/CARS/Configurations/rali_bucharest_all_channels_no_pretrigger.ini'
cnv_args['measurement_identifier'] = 'RM'
cnv_args['trim_overflows'] = 2
cnv_args['files_per_sector'] = 1

cnv_args['radiosonde_folder'] = '/mnt/DATA/Big_data/Databases/CARS/Radiosondes/ino'
cnv_args['rsonde_skip_header'] = 4
cnv_args['rsonde_delimiter'] = 'S'
cnv_args['rsonde_column_index'] = [2, 1, 3, 5]
cnv_args['rsonde_geodata'] = [44.50, 26.13, 91.0]
cnv_args['rsonde_column_units'] = ['m', 'hPa', 'C', 'percent']

# ATLAS preprocessor arguments
prs_args = parse_prs()

prs_args['output_folder'] = os.path.join(output_folder, 'atlas_preprocessor')
prs_args['vertical_limit'] = 20.
if isday == True:
    prs_args['exclude_scattering_type'] = ['v', 'f']

# Common plotting arguments
vis_inp = prs_args['output_folder']
use_distance = True
dpi = 150

# Quicklook arguments
qck_args = parse_qck()

qck_args['z_max_zone'] = [0., 1.]
qck_args['y_lims'] = [0., 14.]
qck_args['use_distance'] = True

qck_args['smooth'] = True 
qck_args['smoothing_exponential'] = True
qck_args['smoothing_range'] = [0.150, 14.]
qck_args['half_window'] = [5., 250.]

qck_args['exclude_detection_mode'] = []
qck_args['exclude_scattering_type'] = ['v', 'r']
qck_args['output_folder'] = os.path.join(output_folder, 'atlas_visualizer', 'qck')

# Rayleigh fit arguments
ray_args = parse_ray()

ray_args['x_lims'] = [0., prs_args['vertical_limit']]
ray_args['x_tick'] = 2.

ray_args['reference_height'] = 7.
ray_args['half_reference_window'] = 500.

ray_args['smooth'] = True 
ray_args['smoothing_range'] = [1., prs_args['vertical_limit']]
ray_args['half_window'] = [100., 100.]

ray_args['exclude_detection_mode'] = []
ray_args['exclude_channel_subtype'] = ['w', 'c']
ray_args['output_folder'] = os.path.join(output_folder, 'atlas_visualizer', 'ray')

# Telecover
tlc_args = parse_tlc()

tlc_args['x_lims'] = [0., 3.]
tlc_args['x_tick'] = 0.5


tlc_args['normalization_height'] = 1.2
tlc_args['half_normalization_window'] = 100.

tlc_args['smooth'] = True 
tlc_args['smoothing_range'] = [0.2, 3.]
tlc_args['half_window'] = [50., 50.]

tlc_args['exclude_detection_mode'] = []
tlc_args['output_folder'] = os.path.join(output_folder, 'atlas_visualizer', 'tlc')

# Polarization calibration
pcl_args = parse_pcl()

pcl_args['channels_r_pcl'] = None
pcl_args['channels_t_pcl'] = None

pcl_args['transmission_ratio_pcl'] = None

pcl_args['calibration_height'] = 3.5 
pcl_args['half_calibration_window'] = 500.
pcl_args['rayleigh_height'] = 3.5 
pcl_args['half_rayleigh_window'] = 500.

pcl_args['smooth'] = True 
pcl_args['smoothing_range'] = [1., prs_args['vertical_limit']]
pcl_args['half_window'] = [500., 500.]

pcl_args['output_folder'] = os.path.join(output_folder, 'atlas_visualizer', 'pcl')

# Call scc_converter sequence
processing_chain.scc_converter(cnv_args = cnv_args, 
                               reprocess = reprocess['scc_converter'],
                               process_cnv = process)

# Call Rayleigh Fit sequence
processing_chain.QA_test(input_folder = cnv_args['output_folder'], 
                         prs_args = prs_args, 
                         vis_args = ray_args, 
                         test_type = 'ray', 
                         process = process['rayleigh_fit'],
                         reprocess_prs = reprocess['atlas_preprocessor'], 
                         reprocess_vis = reprocess['atlas_visualizer'],
                         quicklook = quicklook['rayleigh_fit'],
                         skip = skip_visualizer['rayleigh_fit'])

# Call Telecover Test sequence
processing_chain.QA_test(input_folder = cnv_args['output_folder'], 
                         prs_args = prs_args, 
                         vis_args = tlc_args, 
                         test_type = 'tlc', 
                         process = process['telecover'],
                         reprocess_prs = reprocess['atlas_preprocessor'], 
                         reprocess_vis = reprocess['atlas_visualizer'],  
                         quicklook = quicklook['telecover'],
                         skip = skip_visualizer['telecover'])

# Call Polarization Calibration sequence
processing_chain.QA_test(input_folder = cnv_args['output_folder'], 
                         prs_args = prs_args, 
                         vis_args = pcl_args, 
                         test_type = 'pcl', 
                         process = process['polarization_calibration'],
                         reprocess_prs = reprocess['atlas_preprocessor'], 
                         reprocess_vis = reprocess['atlas_visualizer'],
                         quicklook = quicklook['polarization_calibration'],
                         skip = skip_visualizer['polarization_calibration'])

# Call Quicklook sequence
processing_chain.quicklook(input_folder = prs_args['output_folder'], 
                           vis_args = qck_args, 
                           reprocess = reprocess['atlas_visualizer'],
                           skip = skip_visualizer['quicklook'])