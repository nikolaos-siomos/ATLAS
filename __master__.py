#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob, sys
from helper_functions import processing_chain
from helper_functions import parse_master_args
from helper_functions import read_master_config
from helper_functions.parse_master_args import call_parser as parse_mst
from scc_converter.readers.parse_args import call_parser as parse_cnv
from processor.readers.parse_args import call_parser as parse_prs
from visualizer.readers.parse_qck_args import call_parser as parse_qck
from visualizer.readers.parse_ray_args import call_parser as parse_ray
from visualizer.readers.parse_tlc_args import call_parser as parse_tlc
from visualizer.readers.parse_pcb_args import call_parser as parse_pcb
from helper_functions.parse_master_args import check_parser as check_mst

# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

mst_args = parse_mst()
mst_args = check_mst(mst_args)

mst_cfg = read_master_config.config(mst_args['settings_file'])

if 'isday' in  mst_cfg.gen.keys(): isday = mst_cfg.gen['isday']
else: isday = False

if 'newdata' in  mst_cfg.gen.keys(): newdata = mst_cfg.gen['newdata']
else: newdata = True

if 'visualize' in  mst_cfg.gen.keys(): newdata = mst_cfg.gen['visualize']
else: visualize = True

if 'process' in  mst_cfg.gen.keys(): process = mst_cfg.gen['process']
else: process = ['ray', 'tlc', 'pcb']

if visualize == False:
    processing = {'ray' : False,
                  'tlc' : False,
                  'pcb' : False}

else:
    processing = {'ray' : True,
                  'tlc' : True,
                  'pcb' : True}
    
quicklook = {'ray' : True,
             'tlc' : False,
             'pcb' : False}
    

if 'ray' not in process: processing['ray'] = False
if 'tlc' not in process: processing['tlc'] = False
if 'pcb' not in process: processing['pcb'] = False

if 'process_qck' in mst_cfg.gen.keys():
    process_qck = mst_cfg.gen['process_qck']
    if 'ray' in process_qck: quicklook['ray'] = True
    if 'tlc' in process_qck: quicklook['tlc'] = True
    if 'pcb' in process_qck: quicklook['pcb'] = True

if newdata == False:
    reprocess = {'converter' : False,
                 'preprocessor' : False,
                 'visualizer' : True}
else:
    reprocess = {'converter' : True,
                 'preprocessor' : True,
                 'visualizer' : True}

# Reset argument list
sys.argv = [sys.argv[0]]   
 
# Converter arguments
cnv_args = parse_cnv()

cnv_args = parse_master_args.substitute(org = cnv_args, rpl = mst_args)
cnv_args = parse_master_args.substitute(org = cnv_args, rpl = mst_cfg.cnv)

# Call converter sequence
cnv_files = processing_chain.converter(cnv_args = cnv_args, 
                                       cnv_out = mst_args['converter_out'],
                                       processing = processing,
                                       reprocess = reprocess['converter'])

# Preprocessor arguments
prs_args = parse_prs()

prs_args = parse_master_args.substitute(org = prs_args, rpl = mst_cfg.prs)

if isday == True:
    if prs_args['exclude_channel_type'] == None:
        prs_args['exclude_channel_type'] = ['v', 'f']
    else:
        prs_args['exclude_scattering_type'] = \
            prs_args['exclude_channel_type'] + ['v', 'f']

# Call preprocessor sequence
prs_files = processing_chain.preprocessor(prs_args = prs_args, 
                                          cnv_files = cnv_files,
                                          prs_out = mst_args['preprocessor_out'],
                                          quicklook = quicklook,
                                          processing = processing,
                                          reprocess = reprocess['preprocessor'])

# Rayleigh fit arguments
ray_args = parse_ray()
ray_args = parse_master_args.substitute(org = ray_args, rpl = mst_cfg.ray)

# Call Rayleigh Fit sequence
processing_chain.ray_test(ray_args = ray_args, 
                          prs_files = prs_files, 
                          ray_out = mst_args['visualizer_out'],
                          processing = processing,
                          reprocess = reprocess['visualizer'])

# Telecover arguments
tlc_args = parse_tlc()
tlc_args = parse_master_args.substitute(org = tlc_args, rpl = mst_cfg.tlc)

# Telecover arguments
processing_chain.tlc_test(tlc_args = tlc_args, 
                          prs_files = prs_files, 
                          tlc_out = mst_args['visualizer_out'],
                          processing = processing,
                          reprocess = reprocess['visualizer'])

# Polarization Calibration arguments
pcb_args = parse_pcb()
pcb_args = parse_master_args.substitute(org = pcb_args, rpl = mst_cfg.pcb)

# Polarization Calibration arguments
processing_chain.pcb_test(pcb_args = pcb_args, 
                          prs_files = prs_files,
                          pcb_out = mst_args['visualizer_out'],
                          processing = processing,
                          reprocess = reprocess['visualizer'])

# Quicklook arguments
qck_args = parse_qck()
qck_args = parse_master_args.substitute(org = qck_args, rpl = mst_cfg.qck)

# Telecover arguments
processing_chain.quicklook(qck_args = qck_args, 
                           prs_files = prs_files, 
                           qck_out = mst_args['visualizer_out'],
                           quicklook = quicklook,
                           reprocess = reprocess['visualizer'])
