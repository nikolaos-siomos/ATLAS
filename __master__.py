#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob, sys
from helper_functions import processing_chain
from helper_functions import read_master_config
from helper_functions.parse_master_args import call_parser as parse_mst
from helper_functions.parse_master_args import check_parser as check_mst
import numpy as np

# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

mst_args = parse_mst()
mst_args = check_mst(mst_args)

mst_cfg = read_master_config.config(mst_args['settings_file'])

isday = mst_args['isday']
quick_run = mst_args['quick_run']
process = mst_args['process']
process_qck = mst_args['process_qck']

if not np.isin(np.array(process_qck),np.hstack((process,"off"))).all():
    raise Exception(f"-- Error: The provided process_qck option {process_qck} is not a subset of the provided process option {process}. Please revise the general section of the settings file or use the default values")
    
processing = {'ray' : True,
              'tlc' : True,
              'pcb' : True}
    
quicklook = {'ray' : True,
             'tlc' : True,
             'pcb' : True}
    

if 'ray' not in process: processing['ray'] = False
if 'tlc' not in process: processing['tlc'] = False
if 'pcb' not in process: processing['pcb'] = False
if 'off' in process: 
    processing['ray'] = False
    processing['tlc'] = False
    processing['pcb'] = False


if 'ray' not in process_qck: quicklook['ray'] = False
if 'tlc' not in process_qck: quicklook['tlc'] = False
if 'pcb' not in process_qck: quicklook['pcb'] = False
if 'off' in process_qck: 
    quicklook['ray'] = False
    quicklook['tlc'] = False
    quicklook['pcb'] = False

if quick_run == True:
    reprocess = {'converter' : False,
                 'preprocessor' : False,
                 'visualizer' : True}
else:
    reprocess = {'converter' : True,
                 'preprocessor' : True,
                 'visualizer' : True}

# Reset argument list
sys.argv = [sys.argv[0]]   

# Call converter sequence
cnv_files = processing_chain.converter(mst_args = mst_args,
                                       mst_cfg = mst_cfg,  
                                       cnv_out = mst_args['converter_out'],
                                       processing = processing,
                                       reprocess = reprocess['converter'])

# Call preprocessor sequence
prs_files = processing_chain.preprocessor(mst_args = mst_args,
                                          mst_cfg = mst_cfg,   
                                          cnv_files = cnv_files,
                                          prs_out = mst_args['preprocessor_out'],
                                          quicklook = quicklook,
                                          processing = processing,
                                          reprocess = reprocess['preprocessor'],
                                          isday = isday)

# Call Rayleigh Fit sequence
processing_chain.ray_test(mst_cfg = mst_cfg, 
                          prs_files = prs_files, 
                          ray_out = mst_args['visualizer_out'],
                          processing = processing,
                          reprocess = reprocess['visualizer'])


# Telecover arguments
processing_chain.tlc_test(mst_cfg = mst_cfg,  
                          prs_files = prs_files, 
                          tlc_out = mst_args['visualizer_out'],
                          processing = processing,
                          reprocess = reprocess['visualizer'])

# Polarization Calibration arguments
processing_chain.pcb_test(mst_cfg = mst_cfg,  
                          prs_files = prs_files,
                          pcb_out = mst_args['visualizer_out'],
                          processing = processing,
                          reprocess = reprocess['visualizer'])

# Quicklook arguments
processing_chain.quicklook(mst_cfg = mst_cfg,  
                           prs_files = prs_files, 
                           qck_out = mst_args['visualizer_out'],
                           quicklook = quicklook,
                           reprocess = reprocess['visualizer'])
