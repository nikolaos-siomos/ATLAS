#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob
from scc_converter.__scc_converter__ import main as converter
from processor.__preprocessor__ import main as processor
from visualizer.__quicklook__ import main as qck
from visualizer.__rayleigh_fit__ import main as ray
from visualizer.__telecover__ import main as tlc
from visualizer.__polarization_calibration__ import main as pcl 
from scc_converter.readers.parse_args import call_parser as parse_cnv
from processor.readers.parse_args import call_parser as parse_prs
from visualizer.readers.parse_qck_args import call_parser as parse_qck
from visualizer.readers.parse_ray_args import call_parser as parse_ray
from visualizer.readers.parse_tlc_args import call_parser as parse_tlc
from visualizer.readers.parse_pcl_args import call_parser as parse_pcl
# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

# Force reprocessing
reprocess_scc_converter = False
reprocess_atlas_preprocessor = False

parent_folder = '/mnt/DATA/Big_data/Databases/POLIS/data/220916/01'
output_folder = os.path.join(parent_folder, 'out')

# Converter arguments
cnv_args = parse_cnv()

cnv_args['parent_folder'] = parent_folder
cnv_args['output_folder'] = os.path.join(output_folder, 'scc_converter')
cnv_args['config_file'] = '/mnt/DATA/Big_data/Databases/POLIS/data/configs/polis_bucharest_all_channels.ini'
cnv_args['measurement_identifier'] = 'm'
cnv_args['trim_overflows'] = 2
cnv_args['files_per_sector'] = 3

cnv_args['radiosonde_folder'] = '/mnt/DATA/Big_data/Databases/POLIS/data/radiosondes'
cnv_args['rsonde_skip_header'] = 4
cnv_args['rsonde_delimiter'] = 'S'
cnv_args['rsonde_column_index'] = [2, 1, 3, 5]
cnv_args['rsonde_geodata'] = [44.50, 26.13, 91.0]
cnv_args['rsonde_column_units'] = ['m', 'hPa', 'C', 'percent']

# ATLAS preprocessor arguments
prs_args = parse_prs()

prs_args['output_folder'] = os.path.join(output_folder, 'atlas_preprocessor')
prs_args['vertical_limit'] = 20.

quicklook_ray = True
quicklook_tlc = False
quicklook_pcl = False

# Common plotting arguments
vis_inp = prs_args['output_folder']
use_distance = True
dpi = 150

# Quicklook arguments
qck_args = parse_qck()

qck_args['z_max_zone'] = [0., 3.]
qck_args['y_lims'] = [0., prs_args['vertical_limit']]

qck_args['smooth'] = True 
qck_args['smoothing_range'] = [0.5, prs_args['vertical_limit']]
qck_args['half_window'] = [50., 500.]

qck_args['output_folder'] = os.path.join(output_folder, 'atlas_visualizer', 'qck')

# Rayleigh fit arguments
ray_args = parse_ray()

ray_args['x_lims'] = [0., prs_args['vertical_limit']]

ray_args['reference_height'] = 5.5 
ray_args['half_reference_window'] = 500.

ray_args['smooth'] = True 
ray_args['smoothing_range'] = [1., prs_args['vertical_limit']]
ray_args['half_window'] = [500., 500.]

ray_args['output_folder'] = os.path.join(output_folder, 'atlas_visualizer', 'ray')

# Telecover
tlc_args = parse_tlc()

tlc_args['normalization_height'] = 1.
tlc_args['half_normalization_window'] = 200.

tlc_args['smooth'] = True 
tlc_args['smoothing_range'] = [0.2, 4.]
tlc_args['half_window'] = [10., 500.]

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

# QA files
ray_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'ray_*.nc'))
tlc_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'tlc_*.nc'))
pcl_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'pcl_*.nc'))

# ATLAS
qck_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'qck_*ATLAS*.nc'))
ray_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'ray_*ATLAS*.nc'))
tlc_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'tlc_*ATLAS*.nc'))
pcl_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'pcl_*ATLAS*.nc'))

# Ececute scc_converter
if (len(ray_QA_files) == 0 and len(tlc_QA_files) == 0 and \
    len(pcl_QA_files) == 0) or reprocess_scc_converter == True:
    converter(cnv_args)

ray_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'ray_*.nc'))
tlc_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'tlc_*.nc'))
pcl_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'pcl_*.nc'))

# Excecute ATLAS preprocessor
if len(ray_ATLAS_files) == 0 or reprocess_atlas_preprocessor == True:
    for file in ray_QA_files:
        prs_args['input_file'] = file
        prs_args['quicklook'] = quicklook_ray
        processor(prs_args, __version__)

if len(tlc_QA_files) == 0 or reprocess_scc_converter == True:
    for file in tlc_QA_files:
        prs_args['input_file'] = file
        prs_args['quicklook'] = quicklook_tlc
        processor(prs_args, __version__)

if len(pcl_QA_files) == 0 or reprocess_scc_converter == True:
    for file in pcl_QA_files:
        prs_args['input_file'] = file
        prs_args['quicklook'] = quicklook_pcl
        processor(prs_args, __version__)

qck_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'qck_*ATLAS*.nc'))
ray_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'ray_*ATLAS*.nc'))
tlc_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'tlc_*ATLAS*.nc'))
pcl_ATLAS_files = glob.glob(os.path.join(prs_args['output_folder'], 'pcl_*ATLAS*.nc'))

# Excecute ATLAS visualizer
for file in qck_ATLAS_files:
    qck_args['input_file'] = file
    qck(qck_args)
    
for file in ray_ATLAS_files:
    ray_args['input_file'] = file
    ray(ray_args)
    
for file in tlc_ATLAS_files:
    tlc_args['input_file'] = file
    tlc(tlc_args)
        
for file in pcl_ATLAS_files:
    pcl_args['input_file'] = file
    pcl(pcl_args)
    