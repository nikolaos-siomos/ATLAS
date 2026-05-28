#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import warnings
from version import __version__
from collections import defaultdict
from utils.printouts import print_header
from visualizer.check import check_channels
from processor.packaging import collect_metadata
from visualizer.plot_quicklook import generate_plot
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder, slice_by_vertical_scale, smoothing_2D, collect_dict,
    convert_m_to_km, perform_color_reduction, 
    add_plot_metadata, insert_nan_time_gaps, slice_time
    )

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')
            
def generate_quicklooks(data_pack, caller_info, settings_info):
    
    process_qck = caller_info['process_qck']
    
    qa_test_info = defaultdict(dict)

    for key in process_qck:
        if key in data_pack.keys():
    
            print_header(f"Start generating quicklooks ({key})")
            
            # Prepare folders
            prepare_folder(caller_info, pattern = f"qck_{key}")

            # Load arrays
            profiles = data_pack[key]['profile']
            vertical_scale = data_pack[key][caller_info['vertical_scale']]

            # Load settings
            settings = settings_info.copy()
                        
            # Slice time
            profiles, time_sliced = slice_time(profiles, t_lims = settings['t_lims'])
            
            # Add nan values in time gaps
            profiles, _, has_time_gap = insert_nan_time_gaps(profiles, gap_factor = 1.5)

            # Load time after slicing and icluding time gaps
            time = profiles.time.values

            # Convert the range/height units to km 
            vertical_scale = convert_m_to_km(vertical_scale)

            # Check if the parsed channels exist
            channels = check_channels(
                all_channels = profiles.channel.values,
                settings = settings
                )
            
            if len(channels) > 0:
                qa_test_info[key] = {}
            
            # iterate over the channels
            for ch in channels:
                
                print(f"-- channel: {ch}")
        
                ch_d = dict(channel = ch)
                
                sig_ch = profiles.sel(ch_d)
                vertical_scale_ch = vertical_scale.sel(ch_d)
                
                # # Trim the x and y using the x axis limits                
                # sig_ch, vertical_scale_ch, _  = slice_by_vertical_scale(
                #     da = sig_ch,
                #     vertical_scale = vertical_scale_ch,
                #     x_lims = settings['x_lims'], 
                #     )
                y_vals = sig_ch.values
                x_vals = vertical_scale_ch.values

                # Smoothing
                y_vals_sm, _ = smoothing_2D(
                    args = settings, 
                    x_vals = x_vals, 
                    y_vals = y_vals, 
                    err_type = "std"
                    )

                # Gather the metadata that are common for all QA tests in a dictonary
                metadata = collect_metadata(data_pack[key], atlas_channel_id = ch)
                
                # Initialise qa_test_info dictionary
                qa_test_info[key][ch] = {
                    'time_sliced': time_sliced,
                    'has_time_gap': has_time_gap
                    }
                
                # Collect metadata to be added in the plot files
                plot_metadata = collect_dict(
                    data_list = [
                        time_sliced,
                        has_time_gap,
                        __version__,
                        f'qck_{key}'
                        ],
                    data_keys = [
                        'time_sliced',
                        'has_time_gap',
                        'ATLAS_version',
                        'QA_test_ID'
                        ],
                    add_dicts = [settings, metadata]
                    )

                # Load libraris
                lib = Libraries(
                    caller_info = caller_info,
                    metadata = metadata,
                    extra_metadata = {},
                    settings = settings,
                    qa_test_info = {}
                    )
                
                
                # Call GenerateText class
                text_generator = GenerateText(lib = lib)
                
                # Make titles
                qa_test_info[key][ch]['title'] = \
                    text_generator.make_quicklook_title()
                
                # Make filenames
                qa_test_info[key][ch]['filename'] = text_generator.make_filename(
                    qa_test = f'qck_{key}'
                    )
            
                # Make the plot
                qa_test_info[key][ch]['qck_plot_path'] = generate_plot(
                    T = time,
                    X = x_vals,
                    Y = y_vals,
                    args = settings | qa_test_info[key][ch] | caller_info
                    )  
            
                # Perform color reduction        
                perform_color_reduction(
                    color_reduction = True, 
                    plot_path = qa_test_info[key][ch]['qck_plot_path']
                    )

                # Add the metadata to the plot 
                add_plot_metadata(
                    plot_path = qa_test_info[key][ch]['qck_plot_path'], 
                    plot_metadata = plot_metadata
                    )

            print('-----------------------------------------')
            print(' ')

    