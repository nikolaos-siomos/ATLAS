#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import warnings
import numpy as np
from version import __version__
from collections import defaultdict
from utils.printouts import print_header
from visualizer.check import check_channels
from processor.packaging import collect_metadata
from visualizer.plot_background import generate_plot
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder, perform_color_reduction, 
    add_plot_metadata, insert_nan_time_gaps, slice_time
    )

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')


# QA tests for which background time-series plots are generated.
# Update this set to enable or disable background plots for additional tests.
# BACKGROUND_PLOT_QA_TESTS = {"drk", "ray", "ray_pcb"}

def _to_numpy_selected(da):
    """
    Materialize only the already-selected background time series.

    Matplotlib cannot draw a Dask-backed array directly, so a computation is
    still necessary.  The important point is that channel/time selection has
    already happened before this function is called.
    """

    if hasattr(da, "compute"):
        da = da.compute()

    return np.asarray(da.values)
            
def generate_background(data_pack, caller_info, settings_info):
    
    process_qck = caller_info['process_bgd']
    
    qa_test_info = defaultdict(dict)

    for key in process_qck:
        # if key not in BACKGROUND_PLOT_QA_TESTS:
        #     continue

        if key in data_pack.keys():
    
            print_header(f"Start generating background plots ({key})")
            
            # Prepare folders
            prepare_folder(caller_info, pattern = f"bgd_{key}")

            # Load arrays
            background = data_pack[key]['background']
            background_error = data_pack[key]['background_error']
            system_info = data_pack[key]["system_info"]
            channel_info = data_pack[key]["channel_info"]

            # Load settings
            settings = settings_info.copy()
                        
            # Slice time
            background, time_sliced = slice_time(background, t_lims = settings['t_lims'])
            background_error, _ = slice_time(background_error, t_lims = settings['t_lims'])
            
            # Add NaN profiles in time gaps. This keeps time coordinates finite
            # and only inserts NaNs into the lazy profile values.
            background, _, has_time_gap = insert_nan_time_gaps(background, gap_factor = 1.5)
            background_error, _, _ = insert_nan_time_gaps(background_error, gap_factor = 1.5)

            # Load time after slicing and including time gaps. The time coordinate
            # is eager; only the profile values are lazy.
            time = background.time.values

            # Check if the parsed channels exist and apply exclusion options
            channels = check_channels(
                all_channels = background.channel.values,
                settings = settings
                )
            
            if len(channels) > 0:
                qa_test_info[key] = {}
            
            sys_info = dict(zip(system_info.parameters.values, system_info.values))

            # iterate over the channels
            for ch in channels:
                
                print(f"-- channel: {ch}")
        
                ch_d = dict(channel = ch)
                
                ch_info = channel_info.sel({'channel':ch})  
                ch_info_d = dict(zip(ch_info.parameters.values, ch_info.values))
                
                sig_ch = background.sel(ch_d)
                sig_err_ch = background_error.sel(ch_d)

                # The background is one-dimensional after channel selection.
                # Compute only this selected channel time series.
                y_vals = _to_numpy_selected(sig_ch)
                y_err_vals = _to_numpy_selected(sig_err_ch)

                # Gather the metadata that are common for all QA tests in a dictonary
                metadata = collect_metadata(data_pack[key], atlas_channel_id = ch)
                
                # Initialise qa_test_info dictionary
                qa_test_info[key][ch] = {
                    'time_sliced': time_sliced,
                    'has_time_gap': has_time_gap,
                    "atlas_channel_id": ch,
                    "input_qa_test": key
                    }

                plot_metadata = (
                    {
                        **sys_info,
                        **ch_info_d,
                        **settings,
                        "atlas_channel_id": ch,
                        "ATLAS_version": __version__,
                        "QA_test_ID": f"bgd_{key}",
                    }
                )
                
                plot_metadata = dict(sorted(plot_metadata.items()))

                # Load libraris
                lib = Libraries(
                    caller_info = caller_info,
                    metadata = metadata,
                    extra_metadata = {},
                    settings = settings,
                    qa_test_info = {'qa_test': key}
                    )
                
                # Call GenerateText class
                text_generator = GenerateText(lib = lib)
                
                # Make titles
                qa_test_info[key][ch]['title'] = \
                    text_generator.make_background_title()
                
                # Make filenames
                qa_test_info[key][ch]['filename'] = text_generator.make_filename(
                    qa_test = f'bgd_{key}'
                    )
            
                # Make the plot
                qa_test_info[key][ch]['bgd_plot_path'] = generate_plot(
                    T = time,
                    Y = y_vals,
                    Y_E = y_err_vals,
                    args = settings | qa_test_info[key][ch] | caller_info
                    )  
            
                # Perform color reduction        
                perform_color_reduction(
                    color_reduction = True, 
                    plot_path = qa_test_info[key][ch]['bgd_plot_path']
                    )

                # Add the metadata to the plot 
                add_plot_metadata(
                    plot_path = qa_test_info[key][ch]['bgd_plot_path'], 
                    plot_metadata = plot_metadata
                    )

            print('-----------------------------------------')
            print(' ')
