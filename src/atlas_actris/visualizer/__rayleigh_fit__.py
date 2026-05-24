#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings
import numpy as np
from version import __version__
from visualizer import curve_fit
from visualizer import export_ascii
from collections import defaultdict
from utils.printouts import print_header
from visualizer.check import check_channels
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder, slice_by_vertical_scale, smoothing, collect_dict,
    multiply_y_values, convert_m_to_km, perform_color_reduction, add_plot_metadata
    )
from visualizer import plot_rayleigh, plot_rayleigh_mask

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def get_max_channel_height(norm_region, norm_region_flag):
    
    max_channel_height = np.mean(norm_region)
    max_channel_height = str(int(np.round(1E3*max_channel_height, decimals=-2)))
    
    return(max_channel_height)
    
def generate_rayleigh_fit(data_pack, caller_info, settings_info):
    
    process = caller_info['process']
    
    qa_test_info = defaultdict(dict)
    
    for key in ['ray','ray_pcb']:
        if key in process and key in data_pack:
    
            print_header(f'Initializing the Rayleigh fit test ({key})')
            
            # Prepare folders
            prepare_folder(caller_info, pattern = "_ray_")

            # Load arrays
            profiles = data_pack[key]['profile']
            atten_bsc = data_pack[key]['molecular'].sel({'opto_parameters':'atten_bsc'})
            vertical_scale = data_pack[key][caller_info['vertical_scale']]
                                       
            # Load settings
            settings = settings_info['ray'].copy()
            
            # Convert the range/height units to km 
            vertical_scale = convert_m_to_km(vertical_scale)
                        
            # Check if the parsed channels exist
            channels = check_channels(
                all_channels = profiles.channel.values,
                settings = settings
                )
            
            # Iterate over the channels
            for ch in channels:
                print(f"-- channel: {ch}")
        
                ch_d = dict(channel = ch)
                
                channel_settings = settings.copy()
                                        
                sig_ch = profiles.mean('time').sel(ch_d)
                atb_ch = atten_bsc.sel(ch_d)
                vertical_scale_ch = vertical_scale.sel(ch_d)
                                                
                # Trim the x and y1 using the x axis limits                
                sig_ch, vertical_scale_ch, sl_mask  = slice_by_vertical_scale(
                    da = sig_ch, 
                    vertical_scale = vertical_scale_ch, 
                    x_lims = channel_settings['smoothing_range'], 
                    )
                
                y1_vals = sig_ch.values
                x_vals = vertical_scale_ch.values
                
                # Trim the x and y2 using the x axis limits                
                y2_vals = atb_ch.where(sl_mask, drop=True).values
        
                # Smoothing of the y1 array - generates also the corresponding standard deviation
                y1_vals_sm, y1_errs = smoothing(
                    args = channel_settings, 
                    x_vals = x_vals, 
                    y_vals = y1_vals, 
                    err_type = "std"
                    )
                    
                # Smoothing of the y2 array 
                y2_vals_sm, y1_sems = smoothing(
                    args = channel_settings, 
                    x_vals = x_vals, 
                    y_vals = y2_vals, 
                    err_type = "sem"
                    )
                
                # Check for a Rayleigh fit range
                stats, masks = \
                    curve_fit.statistics(
                        y1 = y1_vals,
                        y2 = y2_vals, 
                        y1_err = y1_errs,
                        y1_avg = y1_vals_sm,
                        x = x_vals,
                        keyw_args = channel_settings
                        )
        
        
                # Identify the uppermost range, bin index and related stats there
                norm_region, norm_region_flag, idx = \
                    curve_fit.scan(
                        masks, 
                        user_norm_region = channel_settings['normalization_region'],
                        auto_fit = True,
                        prefered_range = 'far'
                        )
                    
        
                # Isolate the values of the statistics and masks on the normalization region
                stats_norm_region, masks_norm_region = \
                    curve_fit.metrics_norm_region(
                        idx = idx, 
                        stats = stats, 
                        masks = masks
                        )
                
                # Get the maximum channel height
                maximum_channel_height = get_max_channel_height(
                    norm_region= norm_region, 
                    norm_region_flag = norm_region_flag
                    )
        
                # Normalize y1_vals and y1_errs with the normalization factor from the Rayleigh fit test
                y1_vals_nrm, y1_errs_nrm = multiply_y_values(
                    sig = y1_vals_sm, 
                    sig_err = y1_errs, 
                    coef = stats_norm_region['normalization_factor']
                    )
                      
                # Gather the metadata that are common for all QA tests in a dictonary
                metadata = collect_metadata(data_pack[key], atlas_channel_id = ch)

                # Collect metadata to be returned by the QA test
                qa_test_info[key][ch] = collect_dict(
                    data_list = [
                        stats_norm_region, 
                        masks_norm_region, 
                        maximum_channel_height, 
                        norm_region, 
                        norm_region_flag,
                        ],
                    data_keys = [
                        'stats_norm_region',
                        'masks_norm_region', 
                        'maximum_channel_height', 
                        'norm_region', 
                        'norm_flag',
                        ]
                    )
                
                # Collect metadata to be added in the plot files
                plot_metadata = collect_dict(
                    data_list = [
                        maximum_channel_height, 
                        norm_region, 
                        norm_region_flag,
                        __version__,
                        'ray'
                        ],
                    data_keys = [
                        'maximum_channel_height', 
                        'norm_region', 
                        'norm_flag',
                        'ATLAS_version',
                        'QA_test_ID'
                        ],
                    add_dicts = [channel_settings, metadata]
                    )

#------------------------------------------------------------------------------
# Rayleigh Fit
#------------------------------------------------------------------------------  
        
#------------------------------------------------------------------------------  
# Text
                # Load libraris
                lib = Libraries(
                    caller_info = caller_info,
                    metadata = metadata,
                    extra_metadata = {},
                    settings = channel_settings,
                    qa_test_info = qa_test_info[key][ch]
                    )
                
                # Call GenerateText class
                text_generator = GenerateText(lib = lib)
                
                # Make titles
                qa_test_info[key][ch]['title'] = \
                    text_generator.make_rayleigh_fit_title()
                
                # Make filenames
                qa_test_info[key][ch]['filename'] = text_generator.make_filename(
                    qa_test = 'ray'
                    )
                
                # Make ascii file header
                ascii_header = text_generator.make_header_rayleigh_fit()
            
#------------------------------------------------------------------------------  
# Plot
                # Generate the Rayleigh fit plot
                qa_test_info[key][ch]['ray_plot_path'] = plot_rayleigh.generate_plot(
                    X = x_vals, 
                    Y1 = y1_vals_nrm,
                    Y2 = y2_vals_sm,
                    Y1E = y1_errs_nrm,
                    args = metadata | channel_settings | qa_test_info[key][ch] | caller_info
                    ) 
        
                # Perform color reduction        
                perform_color_reduction(
                    color_reduction = caller_info['color_reduction'], 
                    plot_path = qa_test_info[key][ch]['ray_plot_path']
                    )
                
                # Add the metadata to the Rayleigh fit plot  
                add_plot_metadata(
                    plot_path = qa_test_info[key][ch]['ray_plot_path'], 
                    plot_metadata = plot_metadata
                    )
        
                # Export to ascii (Volker's format)        
                export_ascii.rayleigh(
                    dir_out = caller_info['output_folder'], 
                    fname = f"{qa_test_info[key][ch]['filename']}.txt", 
                    alt = x_vals, 
                    atb = y2_vals, 
                    rcs = y1_vals, 
                    header = ascii_header
                    )
                
#------------------------------------------------------------------------------
# Masks
#------------------------------------------------------------------------------  

#------------------------------------------------------------------------------  
# Text
                
                mask_metadata = qa_test_info[key][ch].copy()
                
                mask_metadata['title'] = text_generator.make_rayleigh_fit_mask_title()
                
                mask_metadata['filename'] = text_generator.make_filename(
                    qa_test = 'ray', 
                    extra_type = 'mask'
                    ) 

#------------------------------------------------------------------------------  
# Plot
                # Generate the molecular mask plot
                qa_test_info[key][ch]['ray_mask_plot_path'] = \
                    plot_rayleigh_mask.generate_plot(
                        masks = masks,
                        args = channel_settings | mask_metadata | caller_info
                        )
        
                # Perform color reduction        
                perform_color_reduction(
                    color_reduction = caller_info['color_reduction'], 
                    plot_path = qa_test_info[key][ch]['ray_mask_plot_path']
                    )
        
                # Add the metadata to the molecular mask plot 
                add_plot_metadata(
                    plot_path = qa_test_info[key][ch]['ray_mask_plot_path'], 
                    plot_metadata = plot_metadata
                    )
                                
    return qa_test_info
                
