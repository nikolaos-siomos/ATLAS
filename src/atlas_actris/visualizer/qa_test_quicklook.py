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
from visualizer.plot_quicklook import generate_plot
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder, smoothing_2D, collect_dict,
    convert_m_to_km, perform_color_reduction, 
    add_plot_metadata, insert_nan_time_gaps, slice_time
    )

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def _get_vertical_bin_dim(vertical_scale):
    """Return the vertical/bin dimension name from a 1D vertical scale."""

    if len(vertical_scale.dims) != 1:
        raise ValueError(
            "The selected vertical scale must be 1D after selecting one channel. "
            f"Found dimensions: {vertical_scale.dims}"
        )

    return vertical_scale.dims[0]


def _slice_vertical_scale_only(sig_ch, vertical_scale_ch, x_lims):
    """
    Slice one channel lazily using only the eager vertical scale.

    This intentionally does not inspect the profile values.  Inspecting profile
    validity with operations such as da.notnull().any(dim='time') would trigger
    a calculation over the lazy profile.  For quicklooks, filtering finite
    vertical coordinates within x_lims is enough and keeps pcolormesh happy.
    """

    bin_dim = _get_vertical_bin_dim(vertical_scale_ch)
    x_vals_all = np.asarray(vertical_scale_ch.values)

    if x_lims is None or len(x_lims) == 0:
        mask = np.isfinite(x_vals_all)
    else:
        mask = (
            np.isfinite(x_vals_all)
            & (x_vals_all >= x_lims[0])
            & (x_vals_all <= x_lims[1])
        )

    if not np.any(mask):
        selected_id = None
        for coord_name in ["channel", "pair"]:
            if coord_name in sig_ch.coords:
                try:
                    selected_id = sig_ch.coords[coord_name].values
                except Exception:
                    selected_id = sig_ch.coords[coord_name]
                break

        raise ValueError(
            "No finite vertical-scale bins were found inside x_lims "
            f"({x_lims}) for selection {selected_id}."
        )

    sig_ch = sig_ch.isel({bin_dim: mask})
    vertical_scale_ch = vertical_scale_ch.isel({bin_dim: mask})

    return sig_ch, vertical_scale_ch


def _to_numpy_selected(da):
    """
    Materialize only the already-selected quicklook slice.

    Matplotlib cannot draw a Dask-backed array directly, so a computation is
    still necessary.  The important point is that channel/time/bin selection has
    already happened before this function is called.
    """

    if hasattr(da, "compute"):
        da = da.compute()

    return np.asarray(da.values)
            
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
            system_info = data_pack[key]["system_info"]
            channel_info = data_pack[key]["channel_info"]

            # Load settings
            settings = settings_info.copy()
                        
            # Slice time
            profiles, time_sliced = slice_time(profiles, t_lims = settings['t_lims'])
            
            # Add NaN profiles in time gaps. This keeps time coordinates finite
            # and only inserts NaNs into the lazy profile values.
            profiles, _, has_time_gap = insert_nan_time_gaps(profiles, gap_factor = 1.5)

            # Load time after slicing and including time gaps. The time coordinate
            # is eager; only the profile values are lazy.
            time = profiles.time.values

            # Convert the range/height units to km 
            vertical_scale = convert_m_to_km(vertical_scale)

            # Check if the parsed channels exist and apply exclusion options
            channels = check_channels(
                all_channels = profiles.channel.values,
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
                
                sig_ch = profiles.sel(ch_d)
                vertical_scale_ch = vertical_scale.sel(ch_d)

                # Trim bins before materializing the lazy profile.  This uses
                # only the eager vertical scale, so it does not compute over the
                # full lazy profile array.
                sig_ch, vertical_scale_ch = _slice_vertical_scale_only(
                    sig_ch = sig_ch,
                    vertical_scale_ch = vertical_scale_ch,
                    x_lims = settings['x_lims'],
                    )

                # Matplotlib and the current smoothing functions need NumPy
                # arrays.  Compute only this selected 2D quicklook slice.
                y_vals = _to_numpy_selected(sig_ch)
                x_vals = np.asarray(vertical_scale_ch.values)

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

                plot_metadata = (
                    {
                        **sys_info,
                        **ch_info_d,
                        **settings,
                        "atlas_channel_id": ch,
                        "ATLAS_version": __version__,
                        "QA_test_ID": f"qck_{key}",
                    }
                )
                
                plot_metadata = dict(sorted(plot_metadata.items()))

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
                    Y = y_vals_sm,
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
