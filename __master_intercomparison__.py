#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob
from helper_functions import processing_chain
from helper_functions import parse_intercomparison_args
from helper_functions import read_intercomparison_config
from helper_functions.parse_master_args import call_parser as parse_mst
from helper_functions.parse_master_args import check_parser as check_mst
from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

mst_args = parse_mst()
mst_args = check_mst(mst_args)

mst_cfg = read_intercomparison_config.config(mst_args['settings_file'])

cmp_args = cmp_args.copy()   

# Call Intercomparison sequence
processing_chain.intercomparison(input_folder_1 = input_folder_1,
                                 input_folder_2 = input_folder_2,
                                 vis_args = cmp_args, 
                                 reprocess = True,
                                 skip = False)

def intercomparison(cmp_args, prs_files, reprocess = True):
    
    cmp_args = cmp_args.copy()   
    
    cmp_args['input_file'] = prs_files['ray'][0]   
    
    cmp_args = check_ray(cmp_args)
    
    cmp_out = cmp_args['output_folder']
    os.makedirs(cmp_out, exist_ok = True)
    
    if reprocess == True:
        for file in glob.glob(os.path.join(cmp_out, '*_cmp_*_ATLAS_*.nc')):
            os.remove(file)

    # Preprocessed files - shoulde be zero if reprocess is True
    cmp_file = glob.glob(os.path.join(cmp_out, '*_cmp_*_ATLAS_*.nc'))

    if len(cmp_file) > 1:
        raise Exception(f'More than one rayleigh fit files detected in folder {cmp_out}. Please make sure that only one rayleigh file exists in that folder ')
   
    # Excecute ATLAS visualizer
    if len(cmp_file) == 0:
        view_ray(cmp_args)
        __intercomparison__(cmp_args, __version__)

    return()

(input_folder_1, input_folder_2, vis_args, reprocess = True, skip = False):

    files = glob.glob(os.path.join(vis_args['output_folder'], 'cmp_*.png'))

    ATLAS_files_1 = \
        glob.glob(os.path.join(input_folder_1, 'ray_*_ATLAS_*.nc'))
        
    if len(ATLAS_files_1) > 1:
        raise Exception(f'Too many input rayleigh files ({ATLAS_files_1}) for the intercomparison for lidar_1. Please make sure there is only one rayleigh file in the folder ')

    ATLAS_files_2 = \
        glob.glob(os.path.join(input_folder_2, 'ray_*_ATLAS_*.nc'))

    if len(ATLAS_files_2) > 1:
        raise Exception(f'Too many input rayleigh files ({ATLAS_files_2}) for the intercomparison for lidar_2. Please make sure there is only one rayleigh file in the folder ')
        
    # Excecute ATLAS visualizer
    if (len(files) == 0 or reprocess == True) and skip == False:
                    
        # cleaner.files(vis_args['output_folder'], pattern = 'cmp_*', 
        #               extension = 'png')
        if len(ATLAS_files_1) != 0 and len(ATLAS_files_2) != 0:
            os.makedirs(vis_args['output_folder'], exist_ok = True)
            vis_args['input_files'] = [ATLAS_files_1[0], ATLAS_files_2[0]]
            cmp(vis_args)
    
    return