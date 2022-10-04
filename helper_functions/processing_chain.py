#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob
from helper_functions import cleaner
from scc_converter.__scc_converter__ import main as converter
from processor.__preprocessor__ import main as processor
from visualizer.__quicklook__ import main as qck
from visualizer.__rayleigh_fit__ import main as ray
from visualizer.__telecover__ import main as tlc
from visualizer.__polarization_calibration__ import main as pcl
from visualizer.__intercomparison__ import main as cmp

# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

def scc_converter(cnv_args, process_cnv, reprocess = True):
    
    # QA files
    ray_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'ray_*.nc'))
    tlc_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'tlc_*.nc'))
    pcl_QA_files = glob.glob(os.path.join(cnv_args['output_folder'], 'pcl_*.nc'))

    # Ececute scc_converter
    if (len(ray_QA_files) == 0 or reprocess == True) and \
        process_cnv['rayleigh_fit'] == True:
        
        if os.path.exists(cnv_args['parent_folder']):
            os.makedirs(cnv_args['output_folder'], exist_ok = True)
        
        cnv_ray_args = cnv_args.copy()
        cnv_ray_args['mode'] = 'R'
        
        cleaner.files(cnv_ray_args['output_folder'], pattern = 'ray_*', extension = 'nc')
        converter(cnv_ray_args)

    if (len(tlc_QA_files) == 0 or reprocess == True) and \
        process_cnv['telecover'] == True:
        
        if os.path.exists(cnv_args['parent_folder']):
            os.makedirs(cnv_args['output_folder'], exist_ok = True)       
            
            cnv_tlc_args = cnv_args.copy()
            cnv_tlc_args['mode'] = 'T'
                    
            cleaner.files(cnv_tlc_args['output_folder'], pattern = 'tlc_*', extension = 'nc')
            converter(cnv_tlc_args)
    
    if (len(pcl_QA_files) == 0 or reprocess == True) and \
        process_cnv['polarization_calibration'] == True:
                
        if os.path.exists(cnv_args['parent_folder']):
            os.makedirs(cnv_args['output_folder'], exist_ok = True)
            
        cnv_pcl_args = cnv_args.copy()
        cnv_pcl_args['mode'] = 'C'
        
        if cnv_pcl_args['rayleigh_filename'] == None:
            ray_QA_file = glob.glob(os.path.join(cnv_args['output_folder'], 'ray_*.nc'))
            if len(ray_QA_file) >= 1:
                cnv_pcl_args['rayleigh_filename'] = ray_QA_file[0]
                if len(ray_QA_file) > 1:
                    print(f'-- Warning: More than 1 rayleigh QA files available to link with the calibration. File {ray_QA_file[0]} was automatically selected. Please make sure this is the correct file')
            else:
                print(f"-- No rayleigh QA file detected inside {cnv_args['output_folder']}. Please either process a rayleigh measurement or manually provide the file in that folder in order to proceed")
            
        cleaner.files(cnv_pcl_args['output_folder'], pattern = 'pcl_*', extension = 'nc')
        converter(cnv_pcl_args)
    
    return

def QA_test(input_folder, prs_args, vis_args, test_type, process = True,
            reprocess_prs = True, reprocess_vis = True, quicklook = True,
            skip = False):
    
    # Find QA files and ascociated telecover files
    prs_files = \
        glob.glob(os.path.join(prs_args['output_folder'], 
                               f'{test_type}_*_ATLAS_*.nc'))
        
    if os.path.exists(prs_args['output_folder']):
        prs_qck_files = \
            [os.path.join(prs_args['output_folder'], 
                          f'qck_{os.path.basename(prs_file)[4:]}') 
                          for prs_file in prs_files]
    else:
        prs_qck_files = []
    
        
    cnv_files = glob.glob(os.path.join(input_folder, f'{test_type}_*.nc'))
    
    # Excecute ATLAS preprocessor
    if (len(prs_files) == 0 or reprocess_prs == True) and process == True:
        
        os.makedirs(prs_args['output_folder'], exist_ok = True)
            
        for file in prs_files: 
            if os.path.exists(file): 
                os.remove(file)
        for file in prs_qck_files:
            if os.path.exists(file): 
                os.remove(file)
        for file in cnv_files:
            prs_args['input_file'] = file
            prs_args['quicklook'] = quicklook
            processor(prs_args, __version__)

    prs_files_out = \
        glob.glob(os.path.join(prs_args['output_folder'], f'{test_type}_*ATLAS*.nc'))

    vis_files = glob.glob(os.path.join(vis_args['output_folder'], '{test_type}_*.png'))
        
    if os.path.exists(vis_args['output_folder']):
        vis_qck_files = \
            [os.path.join(prs_args['output_folder'], 
                          f'qck_{os.path.basename(vis_file)[4:]}') 
                          for vis_file in vis_files]
    else:
        vis_qck_files = []
        
    if (len(vis_files) == 0 or reprocess_vis == True) and skip == False and \
        process == True:
            
        for file in vis_files:  
            if os.path.exists(file): 
                os.remove(file)
        for file in vis_qck_files:  
            if os.path.exists(file): 
                os.remove(file)
        for file in prs_files_out:
            vis_args['input_file'] = file
            if test_type == 'ray':
                os.makedirs(vis_args['output_folder'], exist_ok = True)
                ray(vis_args)
            elif test_type == 'tlc':
                os.makedirs(vis_args['output_folder'], exist_ok = True)
                tlc(vis_args)
            elif test_type == 'pcl': 
                os.makedirs(vis_args['output_folder'], exist_ok = True)
                pcl(vis_args)
            else:
                raise Exception(f'Test type {test_type} not recognized please use one of ray, tlc or pcl ')
    
    return

def quicklook(input_folder, vis_args, reprocess = True, skip = False):

    files = glob.glob(os.path.join(vis_args['output_folder'], 'qck_*.png'))

    ATLAS_files = \
        glob.glob(os.path.join(input_folder, 'qck_*_ATLAS_*.nc'))

    # Excecute ATLAS visualizer
    if (len(files) == 0 or reprocess == True) and skip == False:
                    
        # cleaner.files(vis_args['output_folder'], pattern = 'qck_*', 
        #               extension = 'png')
        for file in ATLAS_files:
            os.makedirs(vis_args['output_folder'], exist_ok = True)
            vis_args['input_file'] = file
            qck(vis_args)
    
    return

def intercomparison(input_folder_1, input_folder_2, vis_args, reprocess = True, skip = False):

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