#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 13:51:14 2025

@author: nikos
"""

import datetime, os, glob, base64
import numpy as np
from PIL import Image

def make_filename(scc_station_id = '', data_identifier = '', scc_configuration_id = '', expert_analyst = '', extension = 'html'):
    
    time_now = str(datetime.datetime.now())
    time_now_text = time_now[2:4] + time_now[5:7] + time_now[8:10] + '-' + time_now[11:13] + time_now[14:16]
    
    parts = ['QA_report', scc_station_id, scc_configuration_id, data_identifier, expert_analyst, time_now_text]
    basename = "_".join(str(p) for p in parts if p is not None)

    filename = f'{basename}.{extension}'
    
    return(filename)

def QA_report(plots_folder, html_filename, photon_only = False, export_all = False):
    
    # -----------------------------------------------------------------------------
    # Export to HTML
    # -----------------------------------------------------------------------------
    # Plot landscape or portrait?
    plot_orientation = 'portrait'

    if plot_orientation == 'landscape':
        plot_width = 1000
    elif plot_orientation == 'portrait':
        plot_width = 710

    with open(html_filename, 'w') as f:
        # write html header
        f.write('<!DOCTYPE html><html><head>')
        f.write('<style> #t1 { -moz-tab-size: 4; tab-size: 4; } </style>')

        f.write('</head><body>')
        f.write('\n')
        # write the body of the file
        # Add a Title to the document
        # f.write('<p style="font-family: Helvetica,Arial,sans-serif; font-size: xx-large; font-weight: normal;">')
        f.write(f'<h1>{os.path.basename(html_filename)}</h1>')
        f.write('\n')

        # Quicklook: Get names of qck plots
        qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0532?pa*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0532?ta*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0532?pa*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_1064?ta*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_1064?pa*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0355?ta*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0355?pa*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_1064?tp*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_1064?pp*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0532?pp*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0532?tp*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0355?tp*.png'))
        if len(qck_images) == 0:
            qck_images = glob.glob(os.path.join(plots_folder,'*_qck_ray_0355?pp*.png'))
        
        f.write('<h1>Quicklooks</h1>')
        f.write('\n')
        ch_qck = np.empty(len(qck_images), dtype = object)
        ch_qck_scc = np.empty(len(qck_images), dtype = object)
        for i in range(len(qck_images)):
            im = qck_images[i]
            with open(im, "rb") as image_file:
                data_uri = base64.b64encode(open(im, 'rb').read()).decode('utf-8')
            ch_qck[i] = Image.open(im).text['atlas_channel_id']
            ch_qck_scc[i] = Image.open(im).text['scc_channel_id']
            # x = im.find("_ATLAS")
            f.write(f'<h2>{ch_qck[i]}</h2>')
            f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
            f.write('\n')

        # Rayleigh-Fit plots
        ray_images = np.sort(glob.glob(os.path.join(plots_folder,'*_ray_*.png')))
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
                # x = im.find("_ATLAS")
                f.write(f'<h2>{ch_ray[i]}</h2>')
                f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
                f.write('\n')

        # Telecover Plots
        tlc_images = np.sort(glob.glob(os.path.join(plots_folder,'*_tlc_*.png')))
        tlc_images = [item for item in tlc_images if 'qck_tlc' not in item]

        f.write('<h1>Telecover</h1>')
        f.write('\n')

        ch_tlc = np.empty(len(tlc_images), dtype = object)
        ch_tlc_scc = np.empty(len(tlc_images), dtype = object)
        for i in range(len(tlc_images)):
            im = tlc_images[i]
            with open(im, "rb") as image_file:
                data_uri = base64.b64encode(open(im, 'rb').read()).decode('utf-8')
            ch_tlc[i] = Image.open(im).text['atlas_channel_id']
            ch_tlc_scc[i] = Image.open(im).text['scc_channel_id']
            if ch_tlc[i][6] == 'a' or photon_only or export_all == True:
                # x = im.find("_ATLAS")
                f.write(f'<h2>{ch_tlc[i]}</h2>')
                f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
                f.write('\n')

        # Polarization Calibration Plots
        pcb_images = np.sort(glob.glob(os.path.join(plots_folder,'*_pcb_*.png')))
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
            ch_pcb_r[i] = Image.open(im).text['atlas_channel_id_r']
            ch_pcb_t[i] = Image.open(im).text['atlas_channel_id_t']
            ch_pcb_scc_r[i] = Image.open(im).text['scc_channel_id_r']
            ch_pcb_scc_t[i] = Image.open(im).text['scc_channel_id_t']
            if ch_pcb_r[i][6] == 'p' or (ch_pcb_r[i][6] == 'a' and float(ch_pcb_r[i][:4]) > 900) or export_all == True:
                # x = im.find("_ATLAS")
                f.write(f'<h2>{ch_pcb_r[i]} to {ch_pcb_t[i]}</h2>')
                f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
                f.write('\n')

        # f.write('<h1>Channel List Rayleigh</h1>')
        # f.write('\n')
        # f.write('<h2>ATLAS ID</h2>')
        # for ch in ch_ray:
        #         f.write(f'{ch}')        
        #         f.write('<br>')
        # f.write('\n')
        # f.write('<h2>SCC ID</h2>')
        # for ch in ch_ray_scc:
        #     f.write(f'{ch}')        
        #     f.write('<br>')
        
        # f.write('<h1>Channel List Pol. Cal.</h1>')
        # f.write('\n')
        # f.write('<h2>ATLAS ID</h2>')
        # for i in range(len(ch_pcb_r)):
        #     f.write(f'{ch_pcb_r[i]} / {ch_pcb_r[i]}')        
        #     f.write('<br>')
        # f.write('<h2>SCC ID</h2>')
        # for i in range(len(ch_pcb_scc_r)):
        #     f.write(f'{ch_pcb_scc_r[i]} / {ch_pcb_scc_t[i]}')        
        #     f.write('<br>')
            
        # finish the html page
        f.write('</body></html>')
    f.close()
    
    print('-- QA report exported!\n')

    return()

# # Unblock to convert the html-file in an odt-file (pandoc must be installed on the computer)

#  # Path to pandoc executable
# pandoc_executable = r'c:\util\pandoc-3.1.8\pandoc.exe'

#  # Path to pandoc reference file for the odt output format
# pandoc_reference = r'c:\Users\volker\AppData\Roaming\Pandoc\templates\template.odt'

# pandoc_command = pandoc_executable + ' -f html -t odt --dpi=135 ' +  html_file + ' -o ' + odt_file + ' --reference-doc ' + pandoc_reference
# os.system(pandoc_command)
