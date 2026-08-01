#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings
import numpy as np
import xarray as xr
from version import __version__
from visualizer import plot_dark
from visualizer import export_ascii
from collections import defaultdict
from utils.printouts import print_header
from visualizer.check import check_channels
from scipy.stats import linregress, shapiro
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder, 
    smoothing,
    smoothing_2D,
    collect_dict,
    convert_m_to_km, 
    slice_by_vertical_scale,
    perform_color_reduction, 
    add_plot_metadata,
    )

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

key_translation = {
    "molecular_mask_window": "fit_mask_window",
    "molecular_mask_window_step": "fit_mask_window_step",
    "molecular_mask_region": "fit_mask_region",
    "rsem_threshold": "rsem_threshold",
    "first_derivative_threshold": "first_derivative_threshold",
    "second_derivative_threshold": "second_derivative_threshold",
    "shapiro_wilk_threshold": "shapiro_wilk_threshold",
    "cross_criterion_threshold": "cross_criterion_threshold",
    "durbin_watson_threshold": "durbin_watson_threshold",
}

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

def extract_arrays(ch_d, key, data_pack, data_pack_bc, data_pack_rc, 
                   caller_info, settings):
    
    vertical_scale = data_pack[key][caller_info['vertical_scale']].sel(ch_d)
    vertical_scale_bc = data_pack[key][caller_info['vertical_scale_bc']].sel(ch_d)

    time_info = data_pack[key]['time_info']
    time_info_bc = data_pack_bc[key]['time_info']
    
    if settings['averaging_rate'] == 'raw':
        profiles = _to_numpy_selected(data_pack[key]['profile'].sel(ch_d))
        profiles_bc = _to_numpy_selected(data_pack_bc[key]['profile'].sel(ch_d))
        
    elif settings['averaging_rate'] == 'low_res':
        profiles = _to_numpy_selected(data_pack[key]['profile_low_res'].sel(ch_d))
        profiles_bc = _to_numpy_selected(data_pack_bc[key]['profile_low_res'].sel(ch_d))
    
    elif settings['averaging_rate'] == 'high_res':
        profiles = _to_numpy_selected(data_pack[key]['profile_high_res'].sel(ch_d))
        profiles_bc = _to_numpy_selected(data_pack_bc[key]['profile_high_res'].sel(ch_d))
    
    profiles_rc = data_pack_rc[key]['profile_mean'].sel(ch_d)
    profile_error_rc = data_pack_rc[key]['profile_error_mean'].sel(ch_d)
    vertical_scale_rc = data_pack_rc[key][caller_info['vertical_scale']].sel(ch_d)
    time_info_rc = data_pack_rc[key]['time_info']
    
    if 'ray' in data_pack:
        profiles_ray_rc = data_pack_rc['ray']['profile_mean'].sel(ch_d)
        profile_error_ray_rc = data_pack_rc['ray']['profile_error_mean'].sel(ch_d)
        vertical_scale_ray_rc = data_pack_rc['ray'][caller_info['vertical_scale']].sel(ch_d)
        time_info_ray_rc = data_pack_rc['ray']['time_info']
        
    else:
        profiles_ray_rc = np.nan * xr.zeros_like(profiles_rc)
        profile_error_ray_rc = np.nan * xr.zeros_like(profiles_rc)
        vertical_scale_ray_rc = np.nan * xr.zeros_like(vertical_scale_rc)
        time_info_ray_rc = np.nan
        
    if 'drk_ray' in data_pack:
        profiles_drk_ray_rc = data_pack_rc['drk_ray']['profile_mean'].sel(ch_d)
        profile_error_drk_ray_rc = data_pack_rc['drk_ray']['profile_error_mean'].sel(ch_d)
        vertical_scale_drk_ray_rc = data_pack_rc['drk_ray'][caller_info['vertical_scale']].sel(ch_d)
        time_info_drk_ray_rc = data_pack_rc['drk_ray']['time_info']
    
    elif 'drk_ray' in caller_info['loading_map']:
        alias_key = caller_info['loading_map']['drk_ray'].sel(ch_d)
        
        profiles_drk_ray_rc = data_pack_rc[alias_key]['profile_mean'].sel(ch_d)
        profile_error_drk_ray_rc = data_pack_rc[alias_key]['profile_error_mean'].sel(ch_d)
        vertical_scale_drk_ray_rc = data_pack_rc[alias_key][caller_info['vertical_scale']].sel(ch_d)
        time_info_drk_ray_rc = data_pack_rc[alias_key]['time_info']
   
    else:
        profiles_drk_ray_rc = np.nan * xr.zeros_like(profiles_rc)
        profile_error_drk_ray_rc = np.nan * xr.zeros_like(profiles_rc)
        vertical_scale_drk_ray_rc = np.nan * xr.zeros_like(vertical_scale_rc)    
        time_info_drk_ray_rc = np.nan
    
    # Convert the range/height units to km 
    vertical_scale = convert_m_to_km(vertical_scale)
    vertical_scale_bc = convert_m_to_km(vertical_scale_bc)
    vertical_scale_rc = convert_m_to_km(vertical_scale_rc)
    vertical_scale_ray_rc = convert_m_to_km(vertical_scale_ray_rc)
    vertical_scale_drk_ray_rc = convert_m_to_km(vertical_scale_drk_ray_rc)

    
    y_dict = {
        "av" : np.asarray(profiles.values),
        "bc" : np.asarray(profiles_bc.values),
        "sm" : np.asarray(profiles_bc.values),
        "rc" : np.asarray(profiles_rc.values),
        "rc_er" : np.asarray(profile_error_rc.values),
        "rc_ray" : np.asarray(profiles_ray_rc.values),
        "rc_drk_ray" : np.asarray(profiles_drk_ray_rc.values),
        "rc_ray_er" : np.asarray(profile_error_ray_rc.values),
        }
    
    y_err_dict = {
        "av" : np.nan * np.asarray(profiles.values),
        "bc" : np.nan * np.asarray(profiles_bc.values),
        "sm" : np.nan * np.asarray(profiles_bc.values),
        "rc" : np.asarray(profile_error_rc.values),
        "rc_ray" : np.asarray(profile_error_ray_rc.values),
        "rc_drk_ray" : np.asarray(profile_error_drk_ray_rc.values),
        }
    
    x_dict = {
        "av" : np.asarray(vertical_scale.values),
        "bc" : np.asarray(vertical_scale_bc.values),
        "sm" : np.asarray(vertical_scale_bc.values),
        "rc" : np.asarray(vertical_scale_rc.values),
        "rc_ray" : np.asarray(vertical_scale_ray_rc.values),
        "rc_drk_ray" : np.asarray(vertical_scale_drk_ray_rc.values),
        }
    
    bins_dict = {
        "av" : np.asarray(vertical_scale.bins.values),
        "bc" : np.asarray(vertical_scale_bc.bins.values),
        "sm" : np.asarray(vertical_scale_bc.bins.values),
        "rc" : np.asarray(vertical_scale_rc.bins.values),
        "rc_ray" : np.asarray(vertical_scale_ray_rc.bins.values),
        "rc_drk_ray" : np.asarray(vertical_scale_drk_ray_rc.bins.values),
        }
    
    time_dict = {
        "av" : time_info,
        "bc" : time_info_bc,
        "sm" : time_info_bc,
        "rc" : time_info_rc,
        "rc_ray" : time_info_ray_rc,
        "rc_drk_ray" : time_info_drk_ray_rc,
        
        }
    
    return y_dict, y_err_dict, x_dict, bins_dict, time_dict
    
def smooth_arrays(x_dict, y_dict, settings):
    
    y_dict['sm'], _ = smoothing_2D(
        args = settings, 
        x_vals = x_dict['sm'], 
        y_vals = y_dict['sm'], 
        err_type = "std"
        )

    y_dict['sm'], _ = smoothing_2D(
        args = settings, 
        x_vals = x_dict['sm'], 
        y_vals = y_dict['sm'], 
        err_type = "std"
        )
    
    y_dict['rc'], _ = smoothing(
        args = settings, 
        x_vals = x_dict['rc'], 
        y_vals = y_dict['rc'], 
        err_type = "std"
        )   
    
    y_dict['rc_ray'], _ = smoothing(
        args = settings, 
        x_vals = x_dict['rc_ray'], 
        y_vals = y_dict['rc_ray'], 
        err_type = "std"
        )   
    
    y_dict['rc_drk_ray'], _ = smoothing(
        args = settings, 
        x_vals = x_dict['rc_drk_ray'], 
        y_vals = y_dict['rc_drk_ray'], 
        err_type = "std"
        ) 
    
    return y_dict
    
    
def generate_rayleigh_fit(
        data_pack, data_pack_bc, data_pack_rc, caller_info, settings_info
        ):

    qa_test_info = defaultdict(dict)
    
    if 'drk' not in caller_info['process']:
        return
    
    for key in data_pack:
        
        if key.startswith('drk'):
            
            print_header(f'Initializing the Dark test ({key})')
            
            # Prepare folders
            prepare_folder(
                caller_info, pattern = key, 
                exclude_patterns = [f'_qck_{key}_']
                )
    
            # Load settings
            settings = settings_info.copy()
            
            # Load arrays
            system_info = data_pack[key]["system_info"]
            channel_info = data_pack[key]["channel_info"]
            
            shots = data_pack[key]["system_info"]
        
            profiles = data_pack[key]['profile']
            
            # Check if the parsed channels exist and apply exclusion options
            channels = check_channels(
                all_channels = profiles.channel.values,
                settings = settings
                )
                    
            # iterate over the channels
            for ch in channels:
                print(f"-- channel: {ch}")
                
                ch_d = dict(channel = ch)

                ch_info = channel_info.sel(ch_d)  
                ch_info_d = dict(zip(ch_info.parameters.values, ch_info.values))
                
                channel_settings = settings.copy()
                
                y_dict, y_err_dict, x_dict, bins_dict, time_dict = \
                    extract_arrays(
                        ch_d = ch_d, 
                        key = key, 
                        data_pack = data_pack, 
                        data_pack_bc = data_pack_bc, 
                        data_pack_rc = data_pack_rc, 
                        caller_info = caller_info, 
                        settings = channel_settings
                        )
        
                y_dict = smooth_arrays(
                    x_dict = x_dict, 
                    y_dict = y_dict, 
                    settings = channel_settings, 
                    )
                
                qa_test_info['shots'] = shots.loc[ch_d].mean().values
                         
                zero_bin = get_zero_bin(ch_info)
            
                if channel_settings['far_range'][1] is None:
                    channel_settings['far_range'][1] = x_dict['av'][-1]
                
                # Slice range to perform statistics
                y_region, x_region = slice_by_vertical_scale(
                    da = y_dict['bc'], 
                    vertical_scale = x_dict['bc'], 
                    x_lims = channel_settings['smoothing_range'], 
                    )

                qa_test_info = calculate_statistics(
                    sig = y_region, 
                    ranges = x_region,
                    time = time_dict['bc'],
                    stats = qa_test_info
                    )

                # Gather the metadata that are common for all QA tests in a dictonary
                metadata = collect_metadata(data_pack[key], atlas_channel_id = ch)
                            
                plot_metadata = (
                    {
                        **system_info,
                        **ch_info_d,
                        **settings,
                        "atlas_channel_id": ch,
                        "ATLAS_version": __version__,
                        "QA_test_ID": key,
                    }
                )
                
                plot_metadata = dict(sorted(plot_metadata.items()))

#------------------------------------------------------------------------------
# Dark
#------------------------------------------------------------------------------  
        
#------------------------------------------------------------------------------  
# Text
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
                text_generator.make_quicklook_title()
            
            # Make filenames
            qa_test_info[key][ch]['filename'] = text_generator.make_filename(
                qa_test = f'qck_{key}'
                )

            # Ascii header            
            ascii_header = text_generator.make_header_rayleigh_fit()

        
            # Raw signal plot x and y axis limits
            xlims_av, xlims_range_av, ylims_av = \
                raw_lims(
                    y_dict['av'], 
                    bins = bins_dict['av'],
                    ranges = x_dict['av'],
                    region = ch_info_d['background_range']
                    )     
                
            # Pretrig zoomed
            xlims_bc, xlims_range_bc, ylims_bc = \
                pretrig_lims(
                    y_dict['bc'], 
                    bins = bins_dict['bc'],
                    ranges = x_dict['bc'],
                    zero_bin = zero_bin,
                    region = ch_info_d['background_range']
                    )    
                
            # Zero bin zoomed
            xlims_zb, xlims_range_zb, ylims_zb = \
                zero_bin_lims(
                    y_dict['bc'], 
                    bins = bins_dict['bc'],
                    ranges = x_dict['bc'],
                    zero_bin = zero_bin,
                    )    
                
            # Smoothed zoomed
            xlims_sm, xlims_range_sm, ylims_sm = \
                smoothed_lims(
                    y_dict['sm'], 
                    bins = bins_dict['sm'],
                    ranges = x_dict['sm'],
                    zero_bin = zero_bin,
                    )  

            # RC smoothed
            xlims_rc, ylims_rc = \
                rc_smoothed_lims(
                    sig = y_dict['rc'], 
                    sig_ray = y_dict['rc_drk'], 
                    sig_drk_ray = y_dict['rc_drk_ray'], 
                    ranges = x_dict['rc'],
                    )  

            # Pass all generated scalar or list parameters relevant to the plots to the args dictionary
            qa_test_info[key][ch] = pass_to_args(
                args = qa_test_info[key][ch], 
                data_list = [
                    xlims_av,
                    xlims_range_av,
                    ylims_av,
                    xlims_bc,
                    xlims_range_bc,
                    ylims_bc,
                    xlims_zb,
                    xlims_range_zb,
                    ylims_zb,
                    xlims_sm,
                    xlims_range_sm,
                    ylims_sm,
                    xlims_rc,
                    ylims_rc,
                    ch[6],
                    ],
                data_keys = [
                    'xlims_av',
                    'xlims_range_av',
                    'ylims_av',
                    'xlims_bc',
                    'xlims_range_bc',
                    'ylims_bc',
                    'xlims_zb',
                    'xlims_range_zb',
                    'ylims_zb',
                    'xlims_sm',
                    'xlims_range_sm',
                    'ylims_sm',
                    'xlims_rc',
                    'ylims_rc',
                    'channel_mode',
                    ]
                )
            
            # Make the plot
            qa_test_info[key][ch]['drk_plot_path'] = plot_dark.generate_plot(
                time_dict = time_dict,
                bins_dict = bins_dict,
                x_dict = x_dict,
                y_dict = y_dict,
                y_err_dict = y_err_dict,
                args = channel_settings | qa_test_info[key][ch] | caller_info
                )  
        
            # Perform color reduction        
            perform_color_reduction(
                color_reduction = True, 
                plot_path = qa_test_info[key][ch]['drk_plot_path']
                )

            # Add the metadata to the plot 
            add_plot_metadata(
                plot_path = qa_test_info[key][ch]['drk_plot_path'], 
                plot_metadata = plot_metadata
                )

            # Export to ascii (Volker's format)        
            export_ascii.dark(
                dir_out = caller_info['ascii_folder'], 
                fname = f"{qa_test_info[key][ch]['filename']}.txt", 
                x_dict = x_dict,
                y_dict = y_dict,
                header = ascii_header,
                )

    return qa_test_info
                        

def slice_signal_by_range(sig, ranges, region):
    """
    Slice a 2D signal array by a range interval.

    Parameters
    ----------
    sig : ndarray (time, bins)
        Signal array.
    ranges : ndarray (bins,)
        Range coordinate corresponding to bins.
    region : list or tuple [min_range, max_range]
        Range limits (same units as ranges).

    Returns
    -------
    sig_slice : ndarray
        Sliced signal array.
    ranges_slice : ndarray
        Corresponding sliced ranges.
    """

    mask = (ranges >= region[0]) & (ranges <= region[1])

    sig_slice = sig.where(mask, drop=True)
    ranges_slice = ranges.where(mask, drop=True)

    return sig_slice, ranges_slice

def calculate_statistics(sig, ranges, time_info, stats= None):
    
    time = (time_info.time - time_info.time[0]).dt.seconds.values + \
        1e-6 * (time_info.time - time_info.time[0]).dt.microseconds
        
    if stats is None:
       stats = {}
           
    if ranges.size > 5:
                
        sig_m_t = np.mean(sig, axis = 1)
        sig_m_b = np.mean(sig, axis = 0)
        sig_m = np.mean(sig)
        sig_m_b_c = sig_m_b - sig_m
        
        vert_fit = linregress(
            x = ranges, 
            y = sig_m_b
            )
        
        temp_fit = linregress(x = time.values, y = sig_m_t)
        
        stats['profiles'] = time.size
        stats['bins'] = ranges.size
        stats['sample'] = sig.size
        
        stats['baseline_offset'] = sig_m
        stats['baseline_offset_sdev'] = np.std(sig)
        stats['baseline_offset_sem'] = stats['baseline_offset_sdev'] / np.sqrt(stats['sample'])

        stats['noise_per_bin'] = stats['baseline_offset_sdev'] * np.sqrt(stats['shots'] * stats['block_samples'])

        stats['vert_slope'] = vert_fit[0]
        stats['vert_slope_sign'] = vert_fit[3] <= 0.05
        stats['temp_slope'] = temp_fit[0]
        stats['temp_slope_sign'] = temp_fit[3] <= 0.05
        stats['gaussian_noise'] = shapiro(sig_m_b_c)[1] > 0.05
        
        if stats['gaussian_noise']: stats['gaussian_noise_flag'] = 'Yes'
        else: stats['gaussian_noise_flag'] = 'No'
        
        if stats['vert_slope_sign']: stats['vert_slope_flag'] = 'Significant'
        else: stats['vert_slope_flag'] = 'Insignificant'

        if stats['temp_slope_sign']: stats['temp_slope_flag'] = 'Significant'
        else: stats['temp_slope_flag'] = 'Insignificant'
    
    else:
        print(f"--Warning: Insufficient number of points for the dark test statistics. Please check the provided stats_range parameter: {stats['stats_range']}")
                    
    return(stats)

def get_zero_bin(channel_info):
    
    if channel_info.loc['DAQ_Trigger_Offset'] <= 0:
        zero_bin = -channel_info.loc['DAQ_Trigger_Offset']
    else:
        zero_bin = 0
        
    return zero_bin

def find_region_ind(arr, llim = None, ulim = None):
    
    zone = np.where((arr >= llim) & (arr <= llim))[0]
    
    if llim is None:
       llim = arr[0]

    if ulim is None:
       ulim = arr[-1]
      
    if not (isinstance(llim,int) or isinstance(llim,float)):
        raise Exception("llim must be integer or float")

    if not (isinstance(ulim,int) or isinstance(ulim,float)):
        raise Exception("ulim must be integer or float")

    if not llim < ulim:
        raise Exception("llim must be smaller than ulim")
        
    zone = np.where((arr >= llim) & (arr <= ulim))[0]
    
    if len(zone) > 3:
        lind = zone[0]
        uind = zone[-1]
    else:
        lind = 0
        uind = -1

    return(lind, uind)
    
def get_span(arr, lind = 0, uind = -1, edge_fraction = 0.03):

    first = arr[lind]
    last = arr[uind]
    
    span = (last - first)
    
    edge = np.ceil(edge_fraction * span)        

    return (first, last, span, edge)
    
def raw_lims(sig, bins, ranges, region):
    
    first_bin, last_bin, bin_span, bin_edge = get_span(bins.values)
    
    first_range, last_range, range_span, range_edge = get_span(ranges.values)
    
    xlims_bins = [first_bin - bin_edge, last_bin + bin_edge]
    
    xlims_range = [
        1E-3 * (first_range - range_edge), 
        1E-3 * (last_range + range_edge)
        ]
    
    min_y, max_y, mean_y, y_edge = \
        region_extrema(
            sig = sig,
            bins = bins,
            region = region
            )
    
    ylims = [mean_y - 3.3 * y_edge, mean_y + 3.3 * y_edge]
    
    return(xlims_bins, xlims_range, ylims)

def pretrig_lims(sig, bins, ranges, zero_bin, region):
    
    if region[0] > zero_bin:
        
        last_range = ranges.values[-1]
        
        lind, uind = find_region_ind(ranges.values, llim = last_range - 4000.)
        
        first_bin, last_bin, bin_span, bin_edge = \
            get_span(bins.values, lind, uind)
        
        first_range, last_range, range_span, range_edge = \
            get_span(ranges.values, lind, uind)
        
    else:
        
        lind, uind = find_region_ind(bins.values, llim = 0, ulim = zero_bin)
        
        first_bin, last_bin, bin_span, bin_edge = \
            get_span(bins.values, lind, uind)
        
        first_range, last_range, range_span, range_edge = \
            get_span(ranges.values, lind, uind)
        
    xlims_bins = [first_bin - bin_edge, last_bin + bin_edge]
    
    xlims_range = [
        1E-3 * (first_range - range_edge), 
        1E-3 * (last_range + range_edge)
        ]
    
    min_y, max_y, mean_y, y_edge = \
        region_extrema(
            sig = sig,
            bins = bins,
            region = region
            )
    
    ylims = [- 3.3 * y_edge, 3.3 * y_edge]
    
    return(xlims_bins, xlims_range, ylims)

def zero_bin_lims(sig, bins, ranges, zero_bin):
    
    r_buffer = 300
    l_buffer = 100
    
    if zero_bin < 200:
        l_buffer = 0
        
    lind, uind = find_region_ind(
        bins.values, 
        llim = zero_bin - l_buffer, 
        ulim = zero_bin + r_buffer
        )
    
    first_bin, last_bin, bin_span, bin_edge = \
        get_span(bins.values, lind, uind)
    
    first_range, last_range, range_span, range_edge = \
        get_span(ranges.values, lind, uind)
    
    xlims_bins = [first_bin - bin_edge, last_bin + bin_edge]
    
    xlims_range = [
        1E-3 * (first_range - range_edge), 
        1E-3 * (last_range + range_edge)
        ]
    
    min_y, max_y, mean_y, y_edge = \
        region_extrema(
            sig = sig,
            bins = bins,
            region = xlims_bins
            )
    
    ylims = [- 1.1 * y_edge, 1.1 * y_edge]
            
    return(xlims_bins, xlims_range, ylims)

def smoothed_lims(sig, bins, ranges, zero_bin):
    
    first_bin, last_bin, bin_span, bin_edge = \
        get_span(bins.values)
    
    first_range, last_range, range_span, range_edge = \
        get_span(ranges.values)
        
    xlims_bins = [first_bin - bin_edge, last_bin + bin_edge]
    
    xlims_range = [
        1E-3 * (first_range - range_edge), 
        1E-3 * (last_range + range_edge)
        ]
    
    min_y, max_y, mean_y, y_edge = \
        region_extrema(
            sig = sig,
            bins = bins,
            region = [zero_bin + 500, last_bin]
            )
    
    ylims = [mean_y - 3.3 * y_edge, mean_y + 3.3 * y_edge]
    
    return(xlims_bins, xlims_range, ylims)

def rc_smoothed_lims(sig, sig_ray, sig_drk_ray, ranges):
    
    mask = (sig == sig)
    
    first_range, last_range, range_span, range_edge = \
        get_span(ranges.copy().where(mask, drop=True).values)
            
    xlims_range = [
        1E-3 * (first_range - range_edge), 
        1E-3 * (last_range + range_edge)
        ]
    
    sig_sl = slice_signal_by_range(sig, ranges, region = [5000., ranges[-1].values])
    sig_ray_sl = slice_signal_by_range(sig_ray, ranges, region = [5000., ranges[-1].values])
    sig_drk_ray_sl = slice_signal_by_range(sig_drk_ray, ranges, region = [5000., ranges[-1].values])
    
    max_val = np.nanmax([np.nanmax(sig_sl), np.nanmax(sig_ray_sl), np.nanmax(sig_drk_ray_sl)])
    ylims = [-3. * max_val, 3. * max_val]
    
    return(xlims_range, ylims)
    
def region_extrema(sig, bins, region):
        
    mask_bins = (bins >= region[0]) &\
        (bins <= region[1])
        
    min_y = sig.min('time').where(mask_bins).min().values
    max_y = sig.max('time').where(mask_bins).max().values
    mean_y = sig.where(mask_bins).mean().values
    
    edge = np.max([mean_y - min_y, max_y - mean_y])

    return(min_y, max_y, mean_y, edge)
    
    
def pass_to_args(args, data_list, data_keys):
    
    for i in range(len(data_keys)):
        args[data_keys[i]] = data_list[i]
        
    return(args)

def add_extra_plot_metadata(plot_metadata, norm_region_flag, 
                            stats_norm_region, maximum_channel_height):
    
    plot_metadata['norm_region_flag'] = f"{norm_region_flag}"
    for key in stats_norm_region.keys():
        plot_metadata[f"stats_{key}"] = f"{stats_norm_region[key]}"
    for key in stats_norm_region.keys():
        plot_metadata[f"masks_{key}"] = f"{stats_norm_region[key]}"
        
    plot_metadata['maximum_channel_height'] = f"{maximum_channel_height}"
    
    return(plot_metadata)
