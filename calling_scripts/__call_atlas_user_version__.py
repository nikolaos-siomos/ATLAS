#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:43:26 2023

@author: nikos, vf
"""

import os, sys, glob
import numpy as np
from datetime import datetime
from PIL import Image
import base64

# -----------------------------------------------------------------------------
# Dictionary with Configuration IDs
# -----------------------------------------------------------------------------
# Station ID (e.g. 'aky', 'brc') 
st_id = ''

# Configuration ID (in the format <lidar_ID>_<version_ID>_<config_ID>)
config_id = ''

# Date 'yyyymmdd' format (e.g. 20240101)
date = ''

# Subfolder name (to be provided if the parent folder is a subfolder of the <config_id>_<date> folder)
subfolder = ''

# File format, choose among: 'licel', 'polly_xt', 'licel_matlab', 'polly_xt_first', 'licel_old2rack'
# licel:          Official Licel format
# polly_xt:       All Polly_XT formats except 1st generation systems
# licel_matlab:   brc, run
# polly_xt_first: evo
# licel_old2rack: cbw
file_format = ''

#Enable/Disable True run (generate only plots) - choose among: True or false
quick_run = False

# Process only a specific test - choose among 'ray', 'tlc', 'tcl_rin', 'pcb'
process = ''

# Process quicklooks for a specific test (must be a subset of process, off or empty string)
process_qck = ''

# Select part of the rayleigh measurement (hhmm format)
slice_rayleigh = ['', '']

# Shortcut of the anaylsing expert (vf, ns, lb, aa)
expert_analyst = ''

# Export all channels in the html
export_all = False

# -----------------------------------------------------------------------------
# User-specific paths (set only once)
# -----------------------------------------------------------------------------
# ATLAS folder 
# atlas_folder = '/home/nikos/Nextcloud/Programs/git/atlas_dev/'
atlas_folder = ''

# Main folder
main_data_folder = ''

# Configuration file
config_file = ''

# Settings file
settings_file = ''

# Radiosonde folder
rsonde_folder = ''

# -----------------------------------------------------------------------------
# Explicit folders per test (non-default for ATLAS)
# Provide only if it is necessary to deviate from the default folder structure
# Otherwise leave as ''
# -----------------------------------------------------------------------------
# Provide the name of the normal folder inside the parent folder (e.g. nrm_02)
nrm = ''

# Provide the name of the pol. cal folder inside the parent folder (e.g. pcb_off)
pcb = ''

# Provide the name of the quadrant telecover folder inside the parent folder (e.g. tlc_raman)
tlc = ''

# Provide the name of the ring telecover folder inside the parent folder (e.g. tlc_rin_01)
tlc_rin = ''

# # Provide the name of the zero bin folder inside the parent folder (e.g. trg_01)
# zeb = ''

# Provide the name of the ring telecover folder inside the parent folder (e.g. tlc_rin_01)
drk = ''

# Parent folder
if len(date) > 0:
    parent_folder = os.path.join(main_data_folder, date, subfolder)
else:
    parent_folders = np.sort(glob.glob(os.path.join(main_data_folder, '*')))
    if len(parent_folders) == 0:
        raise Exception(f'-- Error: No folders found inside the main_data_folder {os.path.join(main_data_folder)}')
    elif len(parent_folders) > 1: 
        raise Exception(f'-- Error: More than one folders found inside the main_data_folder {os.path.join(main_data_folder)}:\n{parent_folders}\n Please provide the date variable in order to select one of them')
    else:
        parent_folder = parent_folders[0]
        parent_folder = os.path.join(parent_folder, subfolder)
        
# Master file
master_file = os.path.join(atlas_folder, '__master__.py')

# -----------------------------------------------------------------------------
# Setup ATLAS arguments
# -----------------------------------------------------------------------------
sys.path.insert(0, atlas_folder)
from __master__ import main as atlas_master
from helper_functions.parse_master_args import call_parser as parse_mst
from helper_functions.parse_master_args import check_parser as check_mst

mst_args = parse_mst()

mst_args['parent_folder']       = parent_folder
mst_args['settings_file']       = settings_file
mst_args['config_file']         = config_file
mst_args['radiosonde']          = rsonde_folder
mst_args['file_format']         = file_format

if len(nrm) > 0:
    mst_args['rayleigh_folder'] = os.path.join(parent_folder, nrm)

if len(pcb) > 0:
    mst_args['pcb_cal_p45_folder'] = os.path.join(parent_folder, pcb, '+45')
    mst_args['pcb_cal_m45_folder'] = os.path.join(parent_folder, pcb, '-45')

if len(tlc) > 0:
    mst_args['telecover_sectors_folder']  = os.path.join(parent_folder, tlc)
    
if len(tlc_rin) > 0:
    mst_args['telecover_rings_folder'] = os.path.join(parent_folder, tlc_rin)

# if len(zeb) > 0:
#     mst_args['zero_bin_folder'] = os.path.join(parent_folder, trg)
    
if len(drk) > 0:
    mst_args['dark_folder'] = os.path.join(parent_folder, drk)

if len(process) > 0:
    mst_args['process'] = process

if len(process_qck) > 0:
    mst_args['process_qck'] = process_qck
        
if quick_run:
    mst_args['quick_run'] = True

if len(slice_rayleigh[0]) != 0 and len(slice_rayleigh[1]) != 0:
    mst_args['slice_rayleigh'] = [slice_rayleigh[0], slice_rayleigh[1]]

mst_args = check_mst(mst_args)
atlas_master(mst_args)

# -----------------------------------------------------------------------------
# Export to HTML
# -----------------------------------------------------------------------------
# Plot landscape or portrait?
plot_orientation = 'portrait'

if plot_orientation == 'landscape':
    plot_width = 1000
elif plot_orientation == 'portrait':
    plot_width = 710

# Get actual date_time for marking the version of the report file
jetzt = str(datetime.now())
now = jetzt[2:4] + jetzt[5:7] + jetzt[8:10] + '-' + jetzt[11:13] + jetzt[14:16]

# write results and plots to html_file
# create filename and open for write
folder_name = os.path.basename(os.path.dirname(parent_folder))
html_file = os.path.join(parent_folder,f'QA_report_{st_id}_{folder_name}_{expert_analyst}_{now}.html')
odt_file = os.path.join(parent_folder,f'QA_report_{st_id}_{folder_name}_{expert_analyst}_{now}.odt')

with open(html_file, 'w') as f:
    # write html header
    f.write('<!DOCTYPE html><html><head>')
    f.write('<style> #t1 { -moz-tab-size: 4; tab-size: 4; } </style>')

    f.write('</head><body>')
    f.write('\n')
    # write the body of the file
    # Add a Title to the document
    # f.write('<p style="font-family: Helvetica,Arial,sans-serif; font-size: xx-large; font-weight: normal;">')
    f.write(f'<h1>QA report_{st_id}_{config_id}_{date}_{expert_analyst}_{now}</h1>')
    f.write('\n')

    # Quicklook: Get names of qck plots
    qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_0532?pa*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_0532?ta*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_0532?pp*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_0532?tp*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_1064?ta*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_1064?pp*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_1064?tp*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_0355?ta*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_0355?pp*.png'))
    if len(qck_images) == 0:
        qck_images = glob.glob(os.path.join(parent_folder,'plots','*_qck_ray_0355?tp*.png'))
    
    f.write('<h1>Quicklooks</h1>')
    f.write('\n')
    ch_qck = np.empty(len(qck_images), dtype = object)
    ch_qck_scc = np.empty(len(qck_images), dtype = object)
    for i in range(len(qck_images)):
        im = qck_images[i]
        with open(im, "rb") as image_file:
            data_uri = base64.b64encode(open(im, 'rb').read()).decode('utf-8')
        ch_qck[i] = Image.open(im).text['channel']
        ch_qck_scc[i] = Image.open(im).text['scc_id']
        x = im.find("_ATLAS")
        f.write(f'<h2>{ch_qck[i]}</h2>')
        f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
        f.write('\n')

    # Rayleigh-Fit plots
    ray_images = np.sort(glob.glob(os.path.join(parent_folder,'plots','*_ray_*.png')))
    ray_images = [item for item in ray_images if 'qck_ray' not in item]
    ray_images = [item for item in ray_images if 'mask' not in item]

    f.write('<h1>Rayleigh Fit</h1>')
    f.write('\n')
    ch_ray = np.empty(len(ray_images), dtype = object)
    ch_ray_scc = np.empty(len(ray_images), dtype = object)
    for i in range(len(ray_images)):
        im = ray_images[i]
        with open(im, "rb") as image_file:
            data_uri = base64.b64encode(open(im, 'rb').read()).decode('utf-8')
        ch_ray[i] = Image.open(im).text['atlas_channel_id']
        ch_ray_scc[i] = Image.open(im).text['scc_channel_id']
        if ch_ray[i][6] == 'p' or (ch_ray[i][6] == 'a' and float(ch_ray[i][:4]) > 900) or export_all == True:
            x = im.find("_ATLAS")
            f.write(f'<h2>{ch_ray[i]}</h2>')
            f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
            f.write('\n')

    # Telecover Plots
    tlc_images = np.sort(glob.glob(os.path.join(parent_folder,'plots','*_tlc_*.png')))
    tlc_images = [item for item in tlc_images if 'qck_tlc' not in item]

    f.write('<h1>Telecover</h1>')
    f.write('\n')

    ch_tlc = np.empty(len(tlc_images), dtype = object)
    ch_tlc_scc = np.empty(len(tlc_images), dtype = object)
    for i in range(len(tlc_images)):
        im = tlc_images[i]
        with open(im, "rb") as image_file:
            data_uri = base64.b64encode(open(im, 'rb').read()).decode('utf-8')
        ch_tlc[i] = Image.open(im).text['channel']
        ch_tlc_scc[i] = Image.open(im).text['scc_id']
        if ch_tlc[i][6] == 'a' or (file_format in ['polly_xt', 'polly_xt_first'])  or export_all == True:
            x = im.find("_ATLAS")
            f.write(f'<h2>{ch_tlc[i]}</h2>')
            f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
            f.write('\n')

    # Polarization Calibration Plots
    pcb_images = np.sort(glob.glob(os.path.join(parent_folder,'plots','*_pcb_*.png')))
    pcb_images = [item for item in pcb_images if 'qck_pcb' not in item]

    f.write('<h1>Polarization Calibration</h1>')
    f.write('\n')

    ch_pcb_r = np.empty(len(pcb_images), dtype = object)
    ch_pcb_t = np.empty(len(pcb_images), dtype = object)
    ch_pcb_scc_r = np.empty(len(pcb_images), dtype = object)
    ch_pcb_scc_t = np.empty(len(pcb_images), dtype = object)
    for i in range(len(pcb_images)):
        im = pcb_images[i]
        with open(im, "rb") as image_file:
            data_uri = base64.b64encode(open(im, 'rb').read()).decode('utf-8')
        ch_pcb_r[i] = Image.open(im).text['channel_r']
        ch_pcb_t[i] = Image.open(im).text['channel_t']
        ch_pcb_scc_r[i] = Image.open(im).text['scc_id_r']
        ch_pcb_scc_t[i] = Image.open(im).text['scc_id_t']
        if ch_pcb_r[i][6] == 'p' or (ch_pcb_r[i][6] == 'a' and float(ch_pcb_r[i][:4]) > 900) or export_all == True:
            x = im.find("_ATLAS")
            f.write(f'<h2>{ch_pcb_r[i]} to {ch_pcb_t[i]}</h2>')
            f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
            f.write('\n')

    f.write('<h1>Channel List Rayleigh</h1>')
    f.write('\n')
    f.write('<h2>ATLAS ID</h2>')
    for ch in ch_ray:
            f.write(f'{ch}')        
            f.write('<br>')
    f.write('\n')
    f.write('<h2>SCC ID</h2>')
    for ch in ch_ray_scc:
        f.write(f'{ch}')        
        f.write('<br>')
    
    f.write('<h1>Channel List Pol. Cal.</h1>')
    f.write('\n')
    f.write('<h2>ATLAS ID</h2>')
    for i in range(len(ch_pcb_r)):
        f.write(f'{ch_pcb_r[i]} / {ch_pcb_r[i]}')        
        f.write('<br>')
    f.write('<h2>SCC ID</h2>')
    for i in range(len(ch_pcb_scc_r)):
        f.write(f'{ch_pcb_scc_r[i]} / {ch_pcb_scc_t[i]}')        
        f.write('<br>')
    # finish the html page
    f.write('</body></html>')
f.close()


# # Unblock to convert the html-file in an odt-file (pandoc must be installed on the computer)

#  # Path to pandoc executable
# pandoc_executable = r'c:\util\pandoc-3.1.8\pandoc.exe'

#  # Path to pandoc reference file for the odt output format
# pandoc_reference = r'c:\Users\volker\AppData\Roaming\Pandoc\templates\template.odt'

# pandoc_command = pandoc_executable + ' -f html -t odt --dpi=135 ' +  html_file + ' -o ' + odt_file + ' --reference-doc ' + pandoc_reference
# os.system(pandoc_command)
