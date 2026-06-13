#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 13:51:14 2025

@author: nikos
"""

import datetime, os, glob, base64
import numpy as np
from PIL import Image
import html
import ast

def make_filename(scc_station_id = '', data_identifier = '', scc_configuration_id = '', expert_analyst = '', extension = 'html'):
    
    time_now = str(datetime.datetime.now())
    time_now_text = time_now[2:4] + time_now[5:7] + time_now[8:10] + '-' + time_now[11:13] + time_now[14:16]
    
    parts = ['QA_report', scc_station_id, scc_configuration_id, data_identifier, expert_analyst, time_now_text]
    basename = "_".join(str(p) for p in parts if p is not None)

    filename = f'{basename}.{extension}'
    
    return(filename)

def bigger_numeric_string(a, b):
    vals = [v for v in (a, b) if v != ""]
    return max(vals, key=int) if vals else ""

def channel_limit_table(f, data, photon_only, export_all):
    
    # Print a list of channels
    title_format = 'font-family: Liberation Serif; font-size: 14pt; font-weight: bold;'
    header_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center; font-weight: bold; background-color: #C5DBF0;'
    units_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center; font-weight: bold; background-color: #E0E0E0;'
    table_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center;'
    footnote_format = 'font-family: Liberation Serif; font-size: 7pt; text-align: center; font-style: italic;'
    
    header = {
        "atlas_channel_id":"ATLAS Channel ID",
        "scc_channel_id":"SCC Channel ID",
        "minimum_channel_height":"Minimum channel height",
        "maximum_channel_height":"Maximum channel height",
        "dead_time":"Dead time",
        "first_signal_rangebin":"First signal rangebin",
        "trigger_delay":"Trigger delay",
        "emitted_wavelength":"Emission Wavelength",
        "detected_wavelength":"Interference filter center",
        "channel_bandwidth":"Interference filter FWHM",
        "G":"G",
        "H":"H"
        }
    
    sub_header = {
        "atlas_channel_id":"",
        "scc_channel_id":"",
        "minimum_channel_height":"based on Telecover Test",
        "maximum_channel_height":"based on Rayleigh Fit",
        "dead_time":"Based on PI input",
        "first_signal_rangebin":"Based on zero bin test",
        "trigger_delay":"Based on zero bin test",
        "emitted_wavelength":"Based on PI input",
        "detected_wavelength":"Based on PI input",
        "channel_bandwidth":"Based on PI input",
        "G":"CARS or manufacturer based on PI input & GHK script",
        "H":"CARS or manufacturer based on PI input & GHK script"
        }
    

    units = {
        "atlas_channel_id":"",
        "scc_channel_id":"",
        "minimum_channel_height":"[m]",
        "maximum_channel_height":"[m]",
        "dead_time":"[ns]",
        "first_signal_rangebin":"[bins]",
        "trigger_delay":"[ns]",
        "emitted_wavelength":"[nm]",
        "detected_wavelength":"[nm]",
        "channel_bandwidth":"[nm]",
        "G":"",
        "H":""
        }

    
    f.write(f'<h2 style="{title_format}">Channel Range Limits - Related Optical Parameters</h2>\n')

    f.write('<table border="1">\n')
    
    # Header row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{header_format}">{html.escape(str(header[key]))}</td>')
    f.write("</tr>\n")

    # Subheader row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{footnote_format}">{html.escape(str(sub_header[key]))}</td>')
    f.write("</tr>\n")  
    
    # Unit row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{units_format}">{html.escape(str(units[key]))}</td>')
    f.write("</tr>\n")
    
    if data['ray']:
        table_metas = data['ray']
    elif data['tlc_qua']:
        table_metas = data['tlc_qua']
    elif data['tlc_rin']:
        table_metas = data['tlc_rin']
    else:
        table_metas = {} 
        
    # for key in header.keys():  
        
    #     table_metas.setdefault(key, "")
            
    for ch, meta in table_metas.items():
        
        atlas_to_scc_triggering(meta)
              
        if data['tlc_qua'].get(ch):
            minimum_channel_height_qua = data['tlc_qua'][ch].get('minimum_channel_height',"")
        else:
            minimum_channel_height_qua = ''
        
        if data['tlc_rin'].get(ch):
            minimum_channel_height_rin = data['tlc_rin'][ch].get('minimum_channel_height',"")
        else:
            minimum_channel_height_rin = ''
            
        meta["minimum_channel_height"] = bigger_numeric_string(
            minimum_channel_height_qua, 
            minimum_channel_height_rin
            )
        
        is_1064 = (ch[6] == 'a' and float(ch[:4]) > 900)
        
        if (ch[6] == 'a' and not is_1064):
            meta["maximum_channel_height"] = ""

        normalization_flag = meta.get("normalization_flag", "")
        if normalization_flag in ["", "external", "default"]:
            meta["maximum_channel_height"] = ""
                
        if ch[6] == 'p' and not photon_only and not export_all:
            meta["minimum_channel_height"] = ""                

        if ch[6] == 'a':
            meta["dead_time"] = ""
               
        # Table row
        f.write("<tr>")
        for key in header.keys():
            f.write(f'<td style="{table_format}">{html.escape(str(meta.get(key,"")))}</td>')
        f.write("</tr>\n")  
        
    f.write('</table>\n')
    f.write('<br clear="all" style="page-break-before:always;">\n')

def channel_background_table(f, data, photon_only, export_all):
    
    # Print a list of channels
    title_format = 'font-family: Liberation Serif; font-size: 14pt; font-weight: bold;'
    header_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center; font-weight: bold; background-color: #C5DBF0;'
    units_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center; font-weight: bold; background-color: #E0E0E0;'
    table_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center;'
    footnote_format = 'font-family: Liberation Serif; font-size: 7pt; text-align: center; font-style: italic;'

    header = {
        "atlas_channel_id":"ATLAS Channel ID",
        "scc_channel_id":"SCC Channel ID",
        "background_low_bin":"Background Low Bin",
        "background_high_bin":"Background High Bin",
        "background_mode":"Background Mode",
        "background_low":"Background Low",
        "background_high":"Background High",
        }
    
    sub_header = {
        "atlas_channel_id":"",
        "scc_channel_id":"SCC Channel ID",
        "background_low_bin":"based on the dark measurement and Rayleigh fit test",
        "background_high_bin":"based on the dark measurement and Rayleigh fit test",
        "background_mode":"checked by CARS",
        "background_low":"in bins/meters if the background Mode is set to Pre-Trigger/Far Field",
        "background_high":"in bins/meters if the background Mode is set to Pre-Trigger/Far Field",
        }
    

    units = {
        "atlas_channel_id":"",
        "scc_channel_id":"",
        "background_low_bin":"[bins]",
        "background_high_bin":"[bins]",
        "background_mode":"[Pre-Trigger or Far Field]",
        "background_low":"[bins or meters]",
        "background_high":"[bins or meters]",
        }

    
    f.write(f'<h2 style="{title_format}">Background Related Signal Parameters</h2>\n')

    f.write('<table border="1">\n')
    
    # Header row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{header_format}">{html.escape(str(header[key]))}</td>')
    f.write("</tr>\n")

    # Subheader row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{footnote_format}">{html.escape(str(sub_header[key]))}</td>')
    f.write("</tr>\n")  
    
    # Unit row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{units_format}">{html.escape(str(units[key]))}</td>')
    f.write("</tr>\n")
    
    if data['ray']:
        table_metas = data['ray']
    elif data['tlc_qua']:
        table_metas = data['tlc_qua']
    elif data['tlc_rin']:
        table_metas = data['tlc_rin']
    else:
        table_metas = {}
            
    for ch, meta in table_metas.items():
        
        atlas_to_scc_triggering(meta)
                    
        meta["background_low_bin"] = str(int(float(meta["background_low_bin"])))
        meta["background_high_bin"] = str(int(float(meta["background_high_bin"])))
        
        meta["background_low"] = str(int(float(meta["background_low"])))
        meta["background_high"] = str(int(float(meta["background_high"])))  
               
        # Table row
        f.write("<tr>")
        for key in header.keys():
            f.write(f'<td style="{table_format}">{html.escape(str(meta.get(key,"")))}</td>')
        f.write("</tr>\n")          
        
    f.write('</table>\n')
    f.write('<br clear="all" style="page-break-before:always;">\n')
   
    
def channel_polarization(f, data, photon_only, export_all):
    
    # Print a list of channels
    title_format = 'font-family: Liberation Serif; font-size: 14pt; font-weight: bold;'
    header_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center; font-weight: bold; background-color: #C5DBF0;'
    units_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center; font-weight: bold; background-color: #E0E0E0;'
    table_format = 'font-family: Liberation Serif; font-size: 10pt; text-align: center;'
    footnote_format = 'font-family: Liberation Serif; font-size: 7pt; text-align: center; font-style: italic;'


    header = {
        "atlas_channel_id_r":"Channel R ATLAS ID",
        "scc_channel_id_r":"Channel R SCC ID",
        "atlas_channel_id_t":"Channel T ATLAS ID",
        "scc_channel_id_t":"Channel T SCC ID",
        "vldr_residual":"Systematic VLDR error",
        "min_bsc_ratio":"Min (backscatter) ratio for PLDR<sup>(a)</sup>",
        "K":"K",
        "R_to_T_transmission_ratio":"ND filter transmission ratio",
        "R_to_T_transmission_ratio_error":"ND filter transmission ratio abs. uncertainty",
        "K_R_to_T_transmission_ratio":"K * ND filter transmission ratio",
        "K_R_to_T_transmission_ratio_error":"K * ND filter transmission ratio abs. uncertainty"
        }
    
    sub_header = {
        "atlas_channel_id_r":"",
        "scc_channel_id_r":"",
        "atlas_channel_id_t":"",
        "scc_channel_id_t":"",
        "vldr_residual":"estimated by CARS from the VLDR residual based on the pol. cal. analysis (&delta;<sup>V</sup><sub>res</sub> on ATLAS plots)",
        "min_bsc_ratio":"minimum backscatter ratio for calculating particle depolarization with an error less than 0.050, based on the pol. cal. analysis",
        "K":"calculated by CARS or manufacturer based on the GHK script and PI input",
        "R_to_T_transmission_ratio":"Pol. Cal. ND filter transmission ratio (R/T) determined by the PI or CARS",
        "R_to_T_transmission_ratio_error":"",
        "K_R_to_T_transmission_ratio":"",
        "K_R_to_T_transmission_ratio_error":""
        }
    

    units = {
        "atlas_channel_id_r":"",
        "scc_channel_id_r":"",
        "atlas_channel_id_t":"",
        "scc_channel_id_t":"",
        "vldr_residual":"",
        "min_bsc_ratio":"",
        "K":"",
        "R_to_T_transmission_ratio":"",
        "R_to_T_transmission_ratio_error":"",
        "K_R_to_T_transmission_ratio":"",
        "K_R_to_T_transmission_ratio_error":""
        }
    
    
    f.write(f'<h2 style="{title_format}">Polarization Correction Parameters and Associated Product Uncertainties</h2>\n')

    f.write('<table border="1">\n')
    
    # Header row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{header_format}">{header[key]}</td>')
    f.write("</tr>\n")

    # Subheader row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{footnote_format}">{sub_header[key]}</td>')
    f.write("</tr>\n")  
    
    # Unit row
    f.write("<tr>")
    for key in header.keys():
        f.write(f'<td style="{units_format}">{html.escape(str(units[key]))}</td>')
    f.write("</tr>\n")
    
    if data['pcb']:
        table_metas = data['pcb']
    else:
        table_metas = {}
    
    for pair, meta in table_metas.items():
        
        meta["vldr_residual"] = str(
            round(float(meta["vldr_residual"]), 4)
            )

        meta["min_bsc_ratio"] = str(
            round(float(meta["sr_limit"]), 3)
            )
        
        meta["R_to_T_transmission_ratio"] = str(
            round(float(meta["R_to_T_transmission_ratio"]), 4)
            )
        
        meta["K_R_to_T_transmission_ratio"] = str(
            round(float(meta["K"]) * float(meta["R_to_T_transmission_ratio"]), 4)
            )
        
        if meta["vldr_residual"] == 'nan':
            meta["vldr_residual"] = ""
            
        if meta["min_bsc_ratio"] == 'nan':
            meta["min_bsc_ratio"] = ""
            
        both_analogue = (meta["atlas_channel_id_r"][6] == 'a' and \
                   meta["atlas_channel_id_t"][6] == 'a')
            
        is_1064 = (both_analogue and \
                       float(meta["atlas_channel_id_t"][:4]) > 900)
        
        if not both_analogue or is_1064:
            # Table row
            f.write("<tr>")
            for key in header.keys():
                f.write(f'<td style="{table_format}">{html.escape(str(meta.get(key,"")))}</td>')
            f.write("</tr>\n")          
    
    f.write('</table>\n')
    f.write('<br clear="all" style="page-break-before:always;">\n')

    
def _encode_png_as_data_uri(path_png: str) -> str:
    """Return base64 string (without the data:image/... prefix) for embedding into HTML."""
    with open(path_png, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode("utf-8")

def qck_patterns(wl, qa):
    
    pattern = [
        f'*_qck_{qa}_{wl}?ta*.png',
        f'*_qck_{qa}_{wl}?pa*.png',
        f'*_qck_{qa}_{wl}?tp*.png',
        f'*_qck_{qa}_{wl}?pp*.png'
        ]
    
    return pattern

def _select_quicklook_images(plots_folder: str, qa):
    """Replicates the original quicklook selection cascade exactly."""

    qck_images_0355 = []
    qck_images_0532 = []
    qck_images_1064 = []
    
    patterns = qck_patterns(wl = '0355', qa = qa)
    for pattern in patterns:
        if len(qck_images_0355) == 0:
                qck_images_0355.extend(
                    glob.glob(os.path.join(plots_folder, pattern))
                    )
    
    patterns = qck_patterns(wl = '0532', qa = qa)
    for pattern in patterns:
        if len(qck_images_0532) == 0:
            qck_images_0532.extend(
                glob.glob(os.path.join(plots_folder, pattern))
                )
    
    patterns = qck_patterns(wl = '1064', qa = qa)
    for pattern in patterns:
        if len(qck_images_1064) == 0:
            qck_images_1064.extend(
                glob.glob(os.path.join(plots_folder, pattern))
                )

    return qck_images_0355 + qck_images_0532 + qck_images_1064

def _select_vldr_images(plots_folder: str):
    """Replicates the original quicklook selection cascade exactly."""

    vldr_images_355 = glob.glob(os.path.join(plots_folder, '*_qck_vldr_0355?va*.png'))
    if len(vldr_images_355) == 0:
        vldr_images_355 = glob.glob(os.path.join(plots_folder, '*_qck_vldr_0355?vp*.png'))

    vldr_images_532 = glob.glob(os.path.join(plots_folder, '*_qck_vldr_0532?va*.png'))
    if len(vldr_images_532) == 0:
        vldr_images_532 = glob.glob(os.path.join(plots_folder, '*_qck_vldr_0532?vp*.png'))
    
    vldr_images_1064 = glob.glob(os.path.join(plots_folder, '*_qck_vldr_1064?va*.png'))
    if len(vldr_images_1064) == 0:
        vldr_images_1064 = glob.glob(os.path.join(plots_folder, '*_qck_vldr_1064?vp*.png'))

    return vldr_images_355 + vldr_images_532 + vldr_images_1064

def _collect_report_data(plots_folder: str):
    """
    Collect all file paths, image metadata and base64 payloads needed to write the report.

    This keeps data collection separate from HTML writing (so that summary tables can
    be written at the beginning of the HTML file).
    """
    data = {}

    # Quicklooks ray
    for qa in ['ray', 'tlc_qua', 'tlc_rin', 'pcb']:
        qck_images = _select_quicklook_images(plots_folder, qa)
        qck_metas = {}
        for im in qck_images:
            meta = Image.open(im).text
            atlas_channel_id = meta['atlas_channel_id']
            qck_metas[atlas_channel_id] = {
                **meta,
                'path':im,
                "data_uri": _encode_png_as_data_uri(im),
            }
        data[f"qck_{qa}"] = qck_metas

    # Rayleigh-Fit plots
    ray_images = np.sort(glob.glob(os.path.join(plots_folder, '*_ray_*.png')))
    ray_images = [item for item in ray_images if '_qck_' not in item]
    ray_images = [item for item in ray_images if '_mask_' not in item]

    ray_metas = {}
    for im in ray_images:
        meta = Image.open(im).text
        atlas_channel_id = meta['atlas_channel_id']
        ray_metas[atlas_channel_id] = {
            **meta,
            'path':im,
            "data_uri": _encode_png_as_data_uri(im),
        }
    data["ray"] = ray_metas

    # Telecover plots
    tlc_qua_images = np.sort(glob.glob(os.path.join(plots_folder, '*_tlc_qua_*.png')))
    tlc_qua_images = [item for item in tlc_qua_images if '_qck_' not in item]

    tlc_rin_images = np.sort(glob.glob(os.path.join(plots_folder, '*_tlc_rin_*.png')))
    tlc_rin_images = [item for item in tlc_rin_images if '_qck_' not in item]

    tlc_qua_metas = {}
    for im in tlc_qua_images:
        meta = Image.open(im).text
        atlas_channel_id = meta['atlas_channel_id']
        tlc_qua_metas[atlas_channel_id] = {
            **meta,
            'path':im,
            "data_uri": _encode_png_as_data_uri(im),
        }
        
    data["tlc_qua"] = tlc_qua_metas

    tlc_rin_metas = {}
    for im in tlc_rin_images:
        meta = Image.open(im).text
        atlas_channel_id = meta['atlas_channel_id']
        tlc_rin_metas[atlas_channel_id] = {
            **meta,
            'path':im,
            "data_uri": _encode_png_as_data_uri(im),
        }
        
    data["tlc_rin"] = tlc_rin_metas

    # Polarization Calibration plots
    pcb_images = np.sort(glob.glob(os.path.join(plots_folder, '*_pcb_*.png')))
    pcb_images = [item for item in pcb_images if '_qck_' not in item and '_ray_' not in item]

    pcb_metas = {}
    for im in pcb_images:
        meta = Image.open(im).text
        atlas_vldr_id = meta['vldr_id']
        pcb_metas[atlas_vldr_id] = {
            **meta,
            'path':im,
            "data_uri": _encode_png_as_data_uri(im),
        }
        
    data["pcb"] = pcb_metas

    # VLDR plots
    vldr_images = _select_vldr_images(plots_folder)

    vldr_metas = {}
    for im in vldr_images:
        meta = Image.open(im).text
        atlas_vldr_id = meta['vldr_id']
        vldr_metas[atlas_vldr_id] = {
            **meta,
            'path':im,
            "data_uri": _encode_png_as_data_uri(im),
        }
        
    data["vldr"] = vldr_metas
    
    return data
 
def channel_entry(f, meta, plot_width):
    
    data_uri = meta["data_uri"]
    f.write(f'<img  width="{plot_width}" src="data:image/png;base64,{data_uri}">')
    f.write('\n')

def QA_report(plots_folder, html_filename, photon_only=False, export_all=False):

    # -----------------------------------------------------------------------------
    # Export to HTML
    # -----------------------------------------------------------------------------
    # Plot landscape or portrait?
    plot_orientation = 'portrait'

    if plot_orientation == 'landscape':
        plot_width = 2000
    elif plot_orientation == 'portrait':
        plot_width = 1420

    # 1) Collect everything first (metadata + embedded images)
    data = _collect_report_data(plots_folder)

    # 2) Write HTML (tables can now be placed at the beginning)
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

        f.write('<br clear="all" style="page-break-before:always;">\n')

        f.write('<h1>Test report summary</h1>')

        # ---------------------------------------------------------------------
        # Summary tables (moved to the beginning; content unchanged)
        # ---------------------------------------------------------------------
        channel_limit_table(f, data, photon_only, export_all)        
        channel_background_table(f, data, photon_only, export_all)
        channel_polarization(f, data, photon_only, export_all)

        f.write('<br clear="all" style="page-break-before:always;">\n')

        # ---------------------------------------------------------------------
        # Plots (HTML content unchanged; only order moved below the tables)
        # ---------------------------------------------------------------------
        # Quicklooks
        f.write('<h1>Quicklooks</h1>')
        f.write('\n')
        for ch, meta in data["qck_ray"].items():
            f.write(f'<h2>{ch}</h2>')
            channel_entry(f, meta, plot_width)
            for qa in ['tlc_qua', 'tlc_rin', 'pcb']:
                if ch in data[f"qck_{qa}"]:
                    f.write('<br>\n')
                    channel_entry(f, data[f"qck_{qa}"][ch], plot_width)
                    
        # Quicklooks VLDR
        f.write('<h1>VLDR Quicklooks</h1>')
        f.write('\n')
        for ch, meta in data["vldr"].items():
            f.write(f'<h2>{ch}</h2>')
            channel_entry(f, meta, plot_width)
        
        # Rayleigh-Fit plots
        f.write('<h1>Rayleigh Fit</h1>')
        f.write('\n')
        for ch, meta in data["ray"].items():
            if ch[6] == 'p' or (ch[6] == 'a' and float(ch[:4]) > 900) or export_all == True:
                f.write(f'<h2>{ch}</h2>')
                channel_entry(f, meta, plot_width)
                
        f.write('<br clear="all" style="page-break-before:always;">\n')

        # Telecover Plots
        f.write('<h1>Telecover</h1>')
        f.write('\n')

        for ch in data["tlc_qua"].keys() | data["tlc_rin"].keys():
            if ch[6] == 'a' or photon_only or export_all == True:
                f.write(f'<h2>{ch}</h2>')
                if ch in data["tlc_qua"]:
                    channel_entry(f, data["tlc_qua"][ch], plot_width)                    
                if ch in data["tlc_qua"] and ch in data["tlc_rin"]:
                    f.write('<br>\n')
                if ch in data["tlc_rin"]:
                    channel_entry(f, data["tlc_rin"][ch], plot_width)

        f.write('<br clear="all" style="page-break-before:always;">\n')

        # Polarization Calibration Plots
        f.write('<h1>Polarization Calibration</h1>')
        f.write('\n')

        for ch, meta in data["pcb"].items():
            ch_r = meta['atlas_channel_id_r']
            ch_t = meta['atlas_channel_id_t']
            if ch_r[6] == 'p' or (ch_r[6] == 'a' and float(ch_r[:4]) > 900) or export_all == True:
                f.write(f'<h2>{ch_r} to {ch_t}</h2>')
                channel_entry(f, meta, plot_width)

    return()


def atlas_to_scc_triggering(meta):

    background_low_bin = int(meta['background_low_bin'])
    background_high_bin = int(meta['background_high_bin'])
    zero_bin = int(meta['zero_bin'])
    range_resolution = float(meta['range_resolution'])

    if zero_bin < -50:
        background_mode = 'Pre-Trigger'
        background_low = background_low_bin
        background_high = background_high_bin
        first_signal_rangebin = -zero_bin
        trigger_delay = -999.
    else:
        background_mode = 'Far Field'
        background_mode = range_resolution * background_low_bin
        background_high = range_resolution * background_high_bin
        if zero_bin < 0:
            first_signal_rangebin = -zero_bin 
            trigger_delay = -999.
        else:
            first_signal_rangebin = -999.
            trigger_delay = zero_bin * range_resolution * 20. / 3.

    meta['background_mode'] = background_mode
    meta['background_low'] = background_low
    meta['background_high'] = background_high
    meta['first_signal_rangebin'] = first_signal_rangebin
    meta['trigger_delay'] = trigger_delay

# # Unblock to convert the html-file in an odt-file (pandoc must be installed on the computer)

#  # Path to pandoc executable
# pandoc_executable = r'c:\util\pandoc-3.1.8\pandoc.exe'

#  # Path to pandoc reference file for the odt output format
# pandoc_reference = r'c:\Users\volker\AppData\Roaming\Pandoc\templates\template.odt'

# pandoc_command = pandoc_executable + ' -f html -t odt --dpi=135 ' +  html_file + ' -o ' + odt_file + ' --reference-doc ' + pandoc_reference
# os.system(pandoc_command)
