#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 13:51:14 2025

@author: nikos
"""


import os
import html
import time
import shutil
import platform
import subprocess
import numpy as np
from PIL import Image
import datetime, os, glob, base64
import zipfile
import tempfile
import re
from pathlib import Path

qck_text_map = {
    'qck_ray': 'Rayleigh measurement', 
    'qck_ray_pcb': 'Rayleigh measurement for polarization calibration',  
    'qck_tlc_qua': 'Quadrant telecover', 
    'qck_tlc_rin': 'Ring telecover', 
    'qck_pcb': 'Polarization calibration', 
    'qck_drk': 'Long dark',
    }


# Keep the table and figure widths tied to the same value.
# Pandoc is used only for the editable DOCX export.
REPORT_CONTENT_WIDTH_PX = 1420
DOCX_PAGE_BREAK_MARKER = "__ATLAS_REPORT_PAGE_BREAK__"

REPORT_CSS = f'''
@page {{
    size: A4 landscape;
    margin: 12mm;
}}
body {{
    font-family: Liberation Sans, sans-serif;
}}
table.report-table {{
    width: {REPORT_CONTENT_WIDTH_PX}px;
    max-width: 100%;
    table-layout: fixed;
    border-collapse: collapse;
}}
table.report-table td {{
    word-wrap: break-word;
    overflow-wrap: break-word;
    vertical-align: middle;
}}
img.report-plot {{
    width: {REPORT_CONTENT_WIDTH_PX}px;
    max-width: 100%;
}}
.report-page-break {{
    display: block;
    clear: both;
    height: 0;
    line-height: 0;
    page-break-before: always;
    break-before: page;
}}
'''


def write_page_break(f):
    """Write a page-break marker that is converted to a real DOCX break."""
    f.write(
        f'<p class="report-page-break" '
        f'style="page-break-before:always; break-before:page; '
        f'font-size:1pt; color:#FFFFFF; margin:0; padding:0;">'
        f'{DOCX_PAGE_BREAK_MARKER}</p>\n'
    )


def write_report_table_start(f):
    """Start a full-width report table."""
    f.write(
        f'<table class="report-table" border="1" '
        f'width="{REPORT_CONTENT_WIDTH_PX}" '
        f'style="width:{REPORT_CONTENT_WIDTH_PX}px; border-collapse:collapse; table-layout:fixed;">\n'
    )

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
    title_format = 'font-family: Liberation Sans; font-size: 20pt; font-weight: bold;'
    header_format = 'font-family: Liberation Sans; font-size: 14pt; text-align: center; font-weight: bold; background-color: #C5DBF0; padding: 3pt;'
    units_format = 'font-family: Liberation Sans; font-size: 13pt; text-align: center; font-weight: bold; background-color: #E0E0E0; padding: 3pt;'
    table_format = 'font-family: Liberation Sans; font-size: 13pt; text-align: center; padding: 3pt;'
    footnote_format = 'font-family: Liberation Sans; font-size: 11pt; text-align: center; font-style: italic; padding: 3pt;'
    
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

    write_report_table_start(f)
    
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
    write_page_break(f)

def channel_background_table(f, data, photon_only, export_all):
    
    # Print a list of channels
    title_format = 'font-family: Liberation Sans; font-size: 20pt; font-weight: bold;'
    header_format = 'font-family: Liberation Sans; font-size: 14pt; text-align: center; font-weight: bold; background-color: #C5DBF0; padding: 3pt;'
    units_format = 'font-family: Liberation Sans; font-size: 13pt; text-align: center; font-weight: bold; background-color: #E0E0E0; padding: 3pt;'
    table_format = 'font-family: Liberation Sans; font-size: 13pt; text-align: center; padding: 3pt;'
    footnote_format = 'font-family: Liberation Sans; font-size: 11pt; text-align: center; font-style: italic; padding: 3pt;'

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

    write_report_table_start(f)
    
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
    write_page_break(f)
   
    
def channel_polarization(f, data, photon_only, export_all):
    
    # Print a list of channels
    title_format = 'font-family: Liberation Sans; font-size: 20pt; font-weight: bold;'
    header_format = 'font-family: Liberation Sans; font-size: 14pt; text-align: center; font-weight: bold; background-color: #C5DBF0; padding: 3pt;'
    units_format = 'font-family: Liberation Sans; font-size: 13pt; text-align: center; font-weight: bold; background-color: #E0E0E0; padding: 3pt;'
    table_format = 'font-family: Liberation Sans; font-size: 13pt; text-align: center; padding: 3pt;'
    footnote_format = 'font-family: Liberation Sans; font-size: 11pt; text-align: center; font-style: italic; padding: 3pt;'


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

    write_report_table_start(f)
    
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
    write_page_break(f)

    
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


def _encode_png_as_data_uri(path_png: str) -> str:
    """Return a full embedded PNG data URI for the self-contained HTML report."""
    with open(path_png, "rb") as image_file:
        encoded = base64.b64encode(image_file.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"

def _collect_report_data(plots_folder: str):
    """
    Collect all file paths and image metadata needed to write the report.

    This keeps data collection separate from HTML writing (so that summary tables can
    be written at the beginning of the HTML file).
    """
    data = {}

    # Quicklooks ray
    for qa in ['ray', 'ray_pcb', 'tlc_qua', 'tlc_rin', 'pcb', 'drk']:
        qck_images = _select_quicklook_images(plots_folder, qa)
        qck_metas = {}
        for im in qck_images:
            meta = Image.open(im).text
            atlas_channel_id = meta['atlas_channel_id']
            qck_metas[atlas_channel_id] = {
                **meta,
                'path':im,
                'data_uri': _encode_png_as_data_uri(im),
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
            'data_uri': _encode_png_as_data_uri(im),
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
            'data_uri': _encode_png_as_data_uri(im),
        }
        
    data["tlc_qua"] = tlc_qua_metas

    tlc_rin_metas = {}
    for im in tlc_rin_images:
        meta = Image.open(im).text
        atlas_channel_id = meta['atlas_channel_id']
        tlc_rin_metas[atlas_channel_id] = {
            **meta,
            'path':im,
            'data_uri': _encode_png_as_data_uri(im),
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
            'data_uri': _encode_png_as_data_uri(im),
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
            'data_uri': _encode_png_as_data_uri(im),
        }
        
    data["vldr"] = vldr_metas
    
    return data
 
def channel_entry(f, meta, plot_width):
    """Write an embedded image tag.

    The public HTML report remains fully self-contained.  The DOCX converter
    creates a temporary Pandoc-only HTML copy with extracted image files, but
    this function must always write embedded images to the real HTML report.
    """

    data_uri = meta["data_uri"]

    # _encode_png_as_data_uri() already returns a full data URI:
    #     data:image/png;base64,...
    # Do not prepend "data:image/png;base64," again, otherwise browsers see
    #     data:image/png;base64,data:image/png;base64,...
    # and show a broken-image placeholder.
    if not str(data_uri).startswith("data:image/"):
        data_uri = f"data:image/png;base64,{data_uri}"

    f.write(
        f'<img class="report-plot" width="{plot_width}" '
        f'style="width:{plot_width}px; max-width:100%;" '
        f'src="{data_uri}">'
    )
    f.write('\n')

def QA_report(
    plots_folder,
    html_filepath,
    photon_only=False,
    export_all=False,
    export_docx=True,
    docx_filepath=None,
    ):

    # -----------------------------------------------------------------------------
    # Export to HTML
    # -----------------------------------------------------------------------------
    # Plot landscape or portrait?
    plot_orientation = 'portrait'

    if plot_orientation == 'landscape':
        plot_width = 2000
    elif plot_orientation == 'portrait':
        plot_width = REPORT_CONTENT_WIDTH_PX

    # 1) Collect everything first (metadata + image file paths)
    data = _collect_report_data(plots_folder)

    # 2) Write HTML (tables can now be placed at the beginning)
    with open(html_filepath, 'w') as f:
        # write html header
        f.write('<!DOCTYPE html><html><head>')
        f.write(f'<style>{REPORT_CSS} #t1 {{ -moz-tab-size: 4; tab-size: 4; }} </style>')

        f.write('</head><body>')
        f.write('\n')
        # write the body of the file
        # Add a Title to the document
        # f.write('<p style="font-family: Helvetica,Arial,sans-serif; font-size: xx-large; font-weight: normal;">')
        f.write(f'<h1>{os.path.basename(html_filepath)}</h1>')

        write_page_break(f)

        f.write('<h1>Test report summary</h1>')

        # ---------------------------------------------------------------------
        # Summary tables (moved to the beginning; content unchanged)
        # ---------------------------------------------------------------------
        channel_limit_table(f, data, photon_only, export_all)        
        channel_background_table(f, data, photon_only, export_all)
        channel_polarization(f, data, photon_only, export_all)

        # ---------------------------------------------------------------------
        # Plots (HTML content unchanged; only order moved below the tables)
        # ---------------------------------------------------------------------
        # Quicklooks
        if data:
            f.write('<h1>Quicklooks</h1>')
            f.write('\n')
            key_first = next((k for k in data if k.startswith("qck_")), None)
            qck_list = ['qck_ray', 'qck_ray_pcb', 'qck_tlc_qua', 'qck_tlc_rin', 'qck_pcb', 'qck_drk']
            for ch, meta in data[key_first].items():
                f.write(f'<h2>{ch}</h2>')
                for key in qck_list:
                    if ch in data[key]:
                        f.write(f'<h3>{qck_text_map[key]}</h3>')
                        f.write('<br>\n')
                        channel_entry(f, data[key][ch], plot_width)
                    
        # Quicklooks VLDR
        if data["vldr"]:
            f.write('<h1>VLDR Quicklooks</h1>')
            f.write('\n')
            for ch, meta in data["vldr"].items():
                f.write(f'<h2>{ch}</h2>')
                channel_entry(f, meta, plot_width)
            
        # Rayleigh-Fit plots
        if data["ray"]:
            f.write('<h1>Rayleigh Fit</h1>')
            f.write('\n')
            for ch, meta in data["ray"].items():
                if ch[6] == 'p' or (ch[6] == 'a' and float(ch[:4]) > 900) or export_all == True:
                    f.write(f'<h2>{ch}</h2>')
                    channel_entry(f, meta, plot_width)
                    
            write_page_break(f)
    
        if data["tlc_qua"] or data["tlc_rin"]:
            # Telecover Plots
            f.write('<h1>Telecover</h1>')
            f.write('\n')

            for ch in  sorted(data["tlc_qua"].keys() | data["tlc_rin"].keys()):
                if ch[6] == 'a' or photon_only or export_all == True:
                    f.write(f'<h2>{ch}</h2>')
                    if ch in data["tlc_qua"]:
                        channel_entry(f, data["tlc_qua"][ch], plot_width)                    
                    if ch in data["tlc_qua"] and ch in data["tlc_rin"]:
                        f.write('<br>\n')
                    if ch in data["tlc_rin"]:
                        channel_entry(f, data["tlc_rin"][ch], plot_width)
    
            write_page_break(f)

        if data["pcb"]:
            # Polarization Calibration Plots
            f.write('<h1>Polarization Calibration</h1>')
            f.write('\n')
    
            for ch, meta in data["pcb"].items():
                ch_r = meta['atlas_channel_id_r']
                ch_t = meta['atlas_channel_id_t']
                if ch_r[6] == 'p' or (ch_r[6] == 'a' and float(ch_r[:4]) > 900) or export_all == True:
                    f.write(f'<h2>{ch_r} to {ch_t}</h2>')
                    channel_entry(f, meta, plot_width)

    # Optional editable document export. Conversion must never break HTML report creation.
    if export_docx:
        try:
            convert_report_data_to_docx(
                data=data,
                docx_filepath=docx_filepath,
                html_filepath=html_filepath,
                photon_only=photon_only,
                export_all=export_all,
                plot_width_px=plot_width,
                verbose=True,
            )
        except Exception as exc:
            print(f"Warning: DOCX export failed. Skipping DOCX export. Exception:\n{exc}")


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





def _safe_text(value):
    """Return a safe text value for HTML/DOCX tables."""
    if value is None:
        return ""
    text = str(value)
    if text.lower() == "nan":
        return ""
    return text


def _copy_meta_with_triggering(meta):
    """Copy metadata and apply SCC triggering fields without mutating the report data."""
    out = dict(meta)
    try:
        atlas_to_scc_triggering(out)
    except Exception:
        pass
    return out


def _int_text(value):
    try:
        return str(int(float(value)))
    except Exception:
        return _safe_text(value)


def _float_round_text(value, ndigits):
    try:
        text = str(round(float(value), ndigits))
        return "" if text.lower() == "nan" else text
    except Exception:
        return _safe_text(value)


def _limit_table_definition():
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
        "H":"H",
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
        "H":"CARS or manufacturer based on PI input & GHK script",
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
        "H":"",
    }
    return header, sub_header, units


def _background_table_definition():
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
    return header, sub_header, units


def _polarization_table_definition():
    header = {
        "atlas_channel_id_r":"Channel R ATLAS ID",
        "scc_channel_id_r":"Channel R SCC ID",
        "atlas_channel_id_t":"Channel T ATLAS ID",
        "scc_channel_id_t":"Channel T SCC ID",
        "vldr_residual":"Systematic VLDR error",
        "min_bsc_ratio":"Min (backscatter) ratio for PLDR(a)",
        "K":"K",
        "R_to_T_transmission_ratio":"ND filter transmission ratio",
        "R_to_T_transmission_ratio_error":"ND filter transmission ratio abs. uncertainty",
        "K_R_to_T_transmission_ratio":"K * ND filter transmission ratio",
        "K_R_to_T_transmission_ratio_error":"K * ND filter transmission ratio abs. uncertainty",
    }
    sub_header = {
        "atlas_channel_id_r":"",
        "scc_channel_id_r":"",
        "atlas_channel_id_t":"",
        "scc_channel_id_t":"",
        "vldr_residual":"estimated by CARS from the VLDR residual based on the pol. cal. analysis",
        "min_bsc_ratio":"minimum backscatter ratio for calculating particle depolarization with an error less than 0.050",
        "K":"calculated by CARS or manufacturer based on the GHK script and PI input",
        "R_to_T_transmission_ratio":"Pol. Cal. ND filter transmission ratio (R/T) determined by the PI or CARS",
        "R_to_T_transmission_ratio_error":"",
        "K_R_to_T_transmission_ratio":"",
        "K_R_to_T_transmission_ratio_error":"",
    }
    units = {key: "" for key in header}
    return header, sub_header, units


def _summary_table_metas(data):
    if data.get('ray'):
        return data['ray']
    if data.get('tlc_qua'):
        return data['tlc_qua']
    if data.get('tlc_rin'):
        return data['tlc_rin']
    return {}


def _limit_table_rows(data, photon_only, export_all):
    header, _, _ = _limit_table_definition()
    rows = []
    for ch, meta0 in _summary_table_metas(data).items():
        meta = _copy_meta_with_triggering(meta0)

        minimum_channel_height_qua = data.get('tlc_qua', {}).get(ch, {}).get('minimum_channel_height', "")
        minimum_channel_height_rin = data.get('tlc_rin', {}).get(ch, {}).get('minimum_channel_height', "")
        meta["minimum_channel_height"] = bigger_numeric_string(minimum_channel_height_qua, minimum_channel_height_rin)

        is_1064 = (ch[6] == 'a' and float(ch[:4]) > 900)
        if ch[6] == 'a' and not is_1064:
            meta["maximum_channel_height"] = ""

        normalization_flag = meta.get("normalization_flag", "")
        if normalization_flag in ["", "external", "default"]:
            meta["maximum_channel_height"] = ""

        if ch[6] == 'p' and not photon_only and not export_all:
            meta["minimum_channel_height"] = ""

        if ch[6] == 'a':
            meta["dead_time"] = ""

        rows.append([_safe_text(meta.get(key, "")) for key in header])
    return rows


def _background_table_rows(data, photon_only, export_all):
    header, _, _ = _background_table_definition()
    rows = []
    for ch, meta0 in _summary_table_metas(data).items():
        meta = _copy_meta_with_triggering(meta0)
        for key in ["background_low_bin", "background_high_bin", "background_low", "background_high"]:
            if key in meta:
                meta[key] = _int_text(meta[key])
        rows.append([_safe_text(meta.get(key, "")) for key in header])
    return rows


def _polarization_table_rows(data, photon_only, export_all):
    header, _, _ = _polarization_table_definition()
    rows = []
    for pair, meta0 in data.get('pcb', {}).items():
        meta = dict(meta0)
        meta["vldr_residual"] = _float_round_text(meta.get("vldr_residual", ""), 4)
        meta["min_bsc_ratio"] = _float_round_text(meta.get("sr_limit", ""), 3)
        meta["R_to_T_transmission_ratio"] = _float_round_text(meta.get("R_to_T_transmission_ratio", ""), 4)
        try:
            meta["K_R_to_T_transmission_ratio"] = str(round(float(meta.get("K", "nan")) * float(meta.get("R_to_T_transmission_ratio", "nan")), 4))
        except Exception:
            meta["K_R_to_T_transmission_ratio"] = ""

        ch_r = meta.get("atlas_channel_id_r", "")
        ch_t = meta.get("atlas_channel_id_t", "")
        if len(ch_r) < 7 or len(ch_t) < 7:
            continue

        both_analogue = (ch_r[6] == 'a' and ch_t[6] == 'a')
        is_1064 = both_analogue and float(ch_t[:4]) > 900

        if not both_analogue or is_1064:
            rows.append([_safe_text(meta.get(key, "")) for key in header])
    return rows


def _set_cell_background(cell, fill_hex):
    """Set a DOCX table-cell background color.

    Parameters
    ----------
    cell : docx.table._Cell
        Target table cell.
    fill_hex : str or None
        RGB hex color, with or without leading '#'.  If empty, no shading is
        applied.
    """
    if not fill_hex:
        return

    fill_hex = str(fill_hex).strip().lstrip('#').upper()
    if not fill_hex:
        return

    try:
        from docx.oxml import OxmlElement
        from docx.oxml.ns import qn

        tc_pr = cell._tc.get_or_add_tcPr()

        # Remove an existing shading element first, so repeated calls do not
        # accumulate duplicate w:shd nodes.
        for shd in tc_pr.findall(qn('w:shd')):
            tc_pr.remove(shd)

        shd = OxmlElement('w:shd')
        shd.set(qn('w:fill'), fill_hex)
        tc_pr.append(shd)
    except Exception:
        pass


def _set_cell_text(cell, text, bold=False, italic=False, font_size_pt=10, fill_hex=None):
    cell.text = ""

    _set_cell_background(cell, fill_hex)

    # Center the cell contents vertically inside the DOCX table cell.
    try:
        from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
    except Exception:
        pass

    paragraph = cell.paragraphs[0]

    # Center the cell contents horizontally.  Keep the integer fallback because
    # python-docx stores WD_ALIGN_PARAGRAPH.CENTER as 1.
    try:
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    except Exception:
        paragraph.alignment = 1

    run = paragraph.add_run(_safe_text(text))
    run.bold = bold
    run.italic = italic
    try:
        from docx.shared import Pt
        run.font.name = "Liberation Sans"
        run.font.size = Pt(font_size_pt)
        for p in cell.paragraphs:
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.space_before = Pt(0)
    except Exception:
        pass


def _add_docx_summary_table(document, title, header, sub_header, units, rows):
    document.add_heading(title, level=2)
    keys = list(header.keys())
    table = document.add_table(rows=3 + len(rows), cols=len(keys))
    table.style = "Table Grid"
    table.autofit = True

    for col, key in enumerate(keys):
        _set_cell_text(table.cell(0, col), header[key], bold=True, font_size_pt=9, fill_hex="C5DBF0")
        _set_cell_text(table.cell(1, col), sub_header[key], italic=True, font_size_pt=7)
        _set_cell_text(table.cell(2, col), units[key], bold=True, font_size_pt=8, fill_hex="E0E0E0")

    for row_idx, row_values in enumerate(rows, start=3):
        for col, value in enumerate(row_values):
            _set_cell_text(table.cell(row_idx, col), value, font_size_pt=8)


def _add_docx_picture(document, image_path, max_width_inches=10.8):
    if not image_path or not os.path.exists(image_path):
        return
    try:
        from docx.shared import Inches
        document.add_picture(image_path, width=Inches(max_width_inches))
    except Exception as exc:
        document.add_paragraph(f"[Could not insert image: {image_path}. {exc}]")


def convert_report_data_to_docx(
    data,
    docx_filepath=None,
    html_filepath=None,
    photon_only=False,
    export_all=False,
    plot_width_px=REPORT_CONTENT_WIDTH_PX,
    verbose=True,
):
    """Create an editable DOCX directly from the collected report data.

    This avoids Pandoc and avoids parsing the large embedded-image HTML.  The
    HTML report remains fully self-contained, while the DOCX embeds the same
    PNG files directly from the plot folder.
    """
    try:
        from docx import Document
        from docx.enum.section import WD_ORIENT
        from docx.shared import Inches, Pt
    except Exception as exc:
        if verbose:
            print(f"Warning: python-docx is not available. Skipping DOCX export. Exception:\n{exc}")
        return None

    if docx_filepath is None:
        if html_filepath is None:
            docx_filepath = os.path.abspath("QA_report.docx")
        else:
            docx_filepath = os.path.splitext(os.path.abspath(html_filepath))[0] + ".docx"
    else:
        docx_filepath = os.path.abspath(docx_filepath)

    outdir = os.path.abspath(os.path.dirname(docx_filepath) or ".")
    os.makedirs(outdir, exist_ok=True)

    document = Document()
    section = document.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width = Inches(11.69)
    section.page_height = Inches(8.27)
    section.top_margin = Inches(0.45)
    section.bottom_margin = Inches(0.45)
    section.left_margin = Inches(0.45)
    section.right_margin = Inches(0.45)

    styles = document.styles
    for style_name in ["Normal", "Title", "Heading 1", "Heading 2", "Heading 3"]:
        try:
            styles[style_name].font.name = "Liberation Sans"
        except Exception:
            pass
    styles["Normal"].font.size = Pt(10)

    title = os.path.basename(html_filepath) if html_filepath else os.path.basename(docx_filepath)
    document.add_heading(title, level=1)
    document.add_page_break()
    document.add_heading("Test report summary", level=1)

    header, sub_header, units = _limit_table_definition()
    _add_docx_summary_table(
        document,
        "Channel Range Limits - Related Optical Parameters",
        header,
        sub_header,
        units,
        _limit_table_rows(data, photon_only, export_all),
    )
    document.add_page_break()

    header, sub_header, units = _background_table_definition()
    _add_docx_summary_table(
        document,
        "Background Related Signal Parameters",
        header,
        sub_header,
        units,
        _background_table_rows(data, photon_only, export_all),
    )
    document.add_page_break()

    header, sub_header, units = _polarization_table_definition()
    _add_docx_summary_table(
        document,
        "Polarization Correction Parameters and Associated Product Uncertainties",
        header,
        sub_header,
        units,
        _polarization_table_rows(data, photon_only, export_all),
    )
    document.add_page_break()

    if data:
        document.add_heading("Quicklooks", level=1)
        key_first = next((k for k in data if k.startswith("qck_")), None)
        qck_list = ['qck_ray', 'qck_ray_pcb', 'qck_tlc_qua', 'qck_tlc_rin', 'qck_pcb', 'qck_drk']
        if key_first is not None:
            for ch, meta in data[key_first].items():
                document.add_heading(ch, level=2)
                for key in qck_list:
                    if ch in data.get(key, {}):
                        document.add_heading(qck_text_map[key], level=3)
                        _add_docx_picture(document, data[key][ch].get('path'))

    if data.get("vldr"):
        document.add_heading("VLDR Quicklooks", level=1)
        for ch, meta in data["vldr"].items():
            document.add_heading(ch, level=2)
            _add_docx_picture(document, meta.get('path'))

    if data.get("ray"):
        document.add_heading("Rayleigh Fit", level=1)
        for ch, meta in data["ray"].items():
            if ch[6] == 'p' or (ch[6] == 'a' and float(ch[:4]) > 900) or export_all:
                document.add_heading(ch, level=2)
                _add_docx_picture(document, meta.get('path'))
        document.add_page_break()

    if data.get("tlc_qua") or data.get("tlc_rin"):
        document.add_heading("Telecover", level=1)
        for ch in sorted(data.get("tlc_qua", {}).keys() | data.get("tlc_rin", {}).keys()):
            if ch[6] == 'a' or photon_only or export_all:
                document.add_heading(ch, level=2)
                if ch in data.get("tlc_qua", {}):
                    _add_docx_picture(document, data["tlc_qua"][ch].get('path'))
                if ch in data.get("tlc_rin", {}):
                    _add_docx_picture(document, data["tlc_rin"][ch].get('path'))
        document.add_page_break()

    if data.get("pcb"):
        document.add_heading("Polarization Calibration", level=1)
        for ch, meta in data["pcb"].items():
            ch_r = meta.get('atlas_channel_id_r', '')
            ch_t = meta.get('atlas_channel_id_t', '')
            if len(ch_r) >= 7 and (ch_r[6] == 'p' or (ch_r[6] == 'a' and float(ch_r[:4]) > 900) or export_all):
                document.add_heading(f"{ch_r} to {ch_t}", level=2)
                _add_docx_picture(document, meta.get('path'))

    document.save(docx_filepath)

    if verbose:
        print(f"DOCX report exported: {docx_filepath}")

    return docx_filepath
