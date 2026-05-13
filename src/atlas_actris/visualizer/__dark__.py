#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings
import numpy as np
import pandas as pd
import xarray as xr
from .readers.parse_drk_args import call_parser, check_parser
from .readers.check import check_channels
from .readers.read_prepro import unpack
from .readers.read_raw import unpack_raw
from .plotting import make_title, plot_utils, plot_dark, plot_dark_mask
from ..visualizer.tools.smoothing import sliding_average_2D_fast
from .tools.smoothing import sliding_average_1D_fast as smooth_1D
from .writters import make_header, export_ascii
from .tools import curve_fit
from scipy.stats import linregress, shapiro
from bokeh.palettes import Category10, turbo
from .tools import curve_fit
from matplotlib import pyplot as plt


# from .processor.lidar_processing.signal import dark_correction

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
    # Check the command line argument information
    args = check_parser(args)

    print('-----------------------------------------')
    print('Initializing the Dark test...')
    print('-----------------------------------------')
    
    profiles, metadata = unpack(args['input_file'])

    if 'ray_file' in args.keys():
        profiles_ray, metadata_ray = unpack(args['ray_file'])
    else:
        profiles_ray = None
        metadata_ray = None
    
    
    args['input_raw_file']
    

    system_info, channel_info, time_info,\
           sig_t, ranges_t, shots = unpack_raw(args['input_raw_file'])
           
    # Check if the parsed channels exist
    channels = \
        check_channels(sel_channels = args['channels'], 
                       all_channels = metadata['atlas_channel_id'],
                       exclude_telescope_type = args['exclude_telescope_type'], 
                       exclude_channel_type = args['exclude_channel_type'], 
                       exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                       exclude_channel_subtype = args['exclude_channel_subtype'])

    meas_duration = int(time_info.loc[:,'Raw_Data_Start_Time'].iloc[-1] / 60.)
    
    print()
    # iterate over the channels
    for ch in channels:
        print(f"-- channel: {ch}")

        args_ch = args.copy()

        ch_d = dict(channel = ch)
        
        ranges_ch, sig_ch, sig_ray_ch = \
            extract_xy_arrays(profiles, profiles_ray, ch_d)
        
        sig_t_ch = sig_t.loc[ch_d]
        ranges_t_ch = ranges_t.loc[ch_d]
        bins_t_ch = sig_t_ch.bins

        # sig_t_up_ch = sig_t_ch.max(dim='time')
        # sig_t_dn_ch = sig_t_ch.min(dim='time')
        
        args_ch['shots'] = shots.loc[ch_d].mean().values
        
        sig_t_av_ch, args_ch['block_samples'] = block_mean_by_duration(
            sig_t_ch, 
            duration = args_ch['averaging_rate']
            )
        
        sig_av_ch = sig_t_ch.mean('time')
                        
        if ch[6] == 'a':
            
            set_background_range(
                args = args_ch, 
                channel_info = channel_info.loc[ch,:]
                )
            
            sig_t_bc_ch, bgr_t_ch = background_correction(
                sig = sig_t_av_ch, 
                background_range = args_ch["background_range"],
                )
            
            # er_per_bin_ch = error_per_bin(sig_t_ch, args_ch, channel_info)
            
            sig_t_sm_ch, sig_t_er_ch = smoothing_2D(
                sig = sig_t_bc_ch, 
                ranges = ranges_t_ch,
                args = args_ch,
                channel_info = channel_info.loc[ch,:],
                mode = 'bin'
                )
            
            sig_t_sm_raw_ch, sig_t_raw_er_ch = smoothing_2D(
                sig = sig_t_ch, 
                ranges = ranges_t_ch,
                args = args_ch,
                channel_info = channel_info.loc[ch,:],
                mode = 'bin'
                )
            
            # sig_t_rc_ch = range_correction(sig_t_sm_ch, ranges_t_ch)
            
            sig_sm_ch, sig_er_ch = smoothing_1D(
                sig = sig_ch,
                ranges = ranges_ch,
                args = args_ch
                )
            
            sig_ray_sm_ch, sig_ray_er_ch = smoothing_1D(
                sig = sig_ray_ch,
                ranges = ranges_ch,
                args = args_ch
                )
        
            last_index = np.where(sig_av_ch == sig_av_ch)[0][-1]
            max_bin = ranges_t_ch.bins.values[last_index] - 1
            
            if args_ch['fit_mask_region'] is None:
                # args['fit_mask_region'] = [
                #     np.ceil(args['fit_mask_window'][0]/2.),
                #     max_bin - np.ceil(args['fit_mask_window'][0]/2.)]
                args_ch['fit_mask_region'] = [0, max_bin]
                
            # Check for a Rayleigh fit range
            stats, masks = \
                curve_fit.statistics(
                    y1 = sig_av_ch.values,
                    y2 = np.ones_like(sig_av_ch), 
                    x  = ranges_t_ch.bins.values,
                    sm_win =  1,
                    keyw_args = args_ch,
                    cancel_stats = ['pos', 'ccr', 'ext', 'psn']
                    )
                
            zero_bin = get_zero_bin(channel_info.loc[ch])
            
            middle_point = masks['total'].middle_point
            pretrig_middle_points = np.where(middle_point <= zero_bin)[0]
            
            if len(pretrig_middle_points) > 0 and zero_bin > 200:
                masks['pretrig'] = masks['total'][:,:pretrig_middle_points[-1]]
            
            if args_ch['far_range'][1] is None:
                args_ch['far_range'][1] = ranges_t_ch[-1].values
            
            if args_ch['far_range'][1] - args_ch['far_range'][0] > 2.:
                far_range_llim = np.where(ranges_t_ch >= 1E3 * args_ch['far_range'][0])[0][0]
                far_range_ulim = np.where(ranges_t_ch <= 1E3 * args_ch['far_range'][1])[0][-1]
    
                masks['far_range'] = masks['total'].loc[:,far_range_llim:far_range_ulim]
            
            # Slice range to perform statistics
            sig_t_sl_ch, ranges_t_sl_ch = slice_signal_by_range(
                sig = sig_t_bc_ch,
                ranges = ranges_t_ch * 1E-3,
                region = args['stats_range']
                )
            
            args_ch = calculate_statistics(
                sig_t_sl_ch, 
                ranges = ranges_t_sl_ch,
                stats = args_ch
                )
        
            # Make title
            title = make_title.dark(
                channel = ch, 
                metadata = metadata, 
                metadata_r = metadata_ray, 
                args = args_ch
                )
            
            # Make plot filename
            fname = plot_utils.make_filename(
                metadata = metadata, 
                channel = ch, 
                meas_type = 'drk', 
                version = __version__
                )
            
            # Raw signal plot x and y axis limits
            xlims_av, xlims_range_av, ylims_av = \
                raw_lims(
                    sig_t_av_ch, 
                    bins = bins_t_ch,
                    ranges = ranges_t_ch,
                    region = args_ch['background_range']
                    )     
                
            # Pretrig zoomed
            xlims_bc, xlims_range_bc, ylims_bc = \
                pretrig_lims(
                    sig_t_bc_ch, 
                    bins = bins_t_ch,
                    ranges = ranges_t_ch,
                    zero_bin = zero_bin,
                    region = args_ch['background_range']
                    )    
                
            # Zero bin zoomed
            xlims_zb, xlims_range_zb, ylims_zb = \
                zero_bin_lims(
                    sig_t_bc_ch, 
                    bins = bins_t_ch,
                    ranges = ranges_t_ch,
                    zero_bin = zero_bin,
                    )    
                
            # Smoothed zoomed
            xlims_sm, xlims_range_sm, ylims_sm = \
                smoothed_lims(
                    sig_t_sm_ch, 
                    bins = bins_t_ch,
                    ranges = ranges_t_ch,
                    zero_bin = zero_bin,
                    )  

            # RC smoothed
            xlims_rc, ylims_rc = \
                rc_smoothed_lims(
                    sig = sig_sm_ch, 
                    sig_er = sig_er_ch, 
                    sig_ray = sig_ray_sm_ch, 
                    sig_ray_er = sig_ray_er_ch, 
                    ranges = ranges_ch,
                    )  

            # Pass all generated scalar or list parameters relevant to the plots to the args dictionary
            args_ch = pass_to_args(
                args = args_ch, 
                data_list = [
                    title,
                    fname,
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
                    meas_duration,
                    ch[6],
                    ],
                data_keys = [
                    'title',
                    'fname',
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
                    'meas_duration',
                    'channel_mode',
                    ]
                )
            
            data_pack = {
                "av" : sig_t_av_ch,
                "bc" : sig_t_bc_ch,
                "sm" : sig_t_sm_ch,
                "rc" : sig_sm_ch,
                "rc_er" : sig_er_ch,
                "rc_ray" : sig_ray_sm_ch,
                "rc_ray_er" : sig_ray_er_ch,
                "ranges" : ranges_t_ch,
                "ranges_rc" : ranges_ch,
                "bins" : bins_t_ch,
                }
            
            drk_plot_path = \
                plot_dark.generate_plot(
                    profiles = data_pack,
                    args = args_ch)
                
            # Make mask plot filename
            fname_mask = plot_utils.make_filename(metadata = metadata, 
                                                  channel = ch, 
                                                  meas_type = 'drk', 
                                                  extra_type = 'mask',
                                                  version = __version__)
        
            # Make title
            title_mask = make_title.dark(
                channel = ch, 
                metadata = metadata, 
                metadata_r = metadata_ray, 
                args = args_ch,
                is_mask = True
                )
            
            # Pass all additional scalar or list parameters relevant to the mask plots to the args dictionary
            args_mask_ch = args_ch.copy()
            
            args_mask_ch = pass_to_args(
                args = args_mask_ch,
                data_list = [
                    title_mask,
                    fname_mask
                    ],
                data_keys = [
                    'title',
                    'fname'
                    ]
                )
            
            # Generate the molecular mask plot
            drk_mask_plot_path = \
                plot_dark_mask.generate_plot(masks = masks,
                                             args = args_mask_ch)

            # Perform color reduction        
            plot_utils.perform_color_reduction(color_reduction = args_ch['color_reduction'], 
                                               plot_path = drk_plot_path)
            
            # Perform color reduction        
            plot_utils.perform_color_reduction(color_reduction = args_ch['color_reduction'], 
                                               plot_path = drk_mask_plot_path)
            

                        

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

def calculate_statistics(sig, ranges, stats= None):
    
    time = (sig.copy().time - sig.copy().time[0]).dt.seconds.values + \
        1e-6 * (sig.copy().time - sig.copy().time[0]).dt.microseconds
        
    if stats is None:
       stats = {}
           
    if ranges.size > 5:
                
        sig_m_t = sig.mean(dim = 'bins')
        sig_m_b = sig.mean(dim = 'time')
        sig_m = sig.mean()
        sig_m_b_c = sig_m_b - sig_m
        
        vert_fit = linregress(x = ranges.values, 
                              y = sig_m_b.values)        
        temp_fit = linregress(x = time.values, y = sig_m_t.values)
        
        stats['profiles'] = time.size
        stats['bins'] = ranges.size
        stats['sample'] = sig.size
        
        stats['baseline_offset'] = sig_m.values
        stats['baseline_offset_sdev'] = sig.std().values
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
        print(f"--Warning: Insufficient number of points for the dark test statistics. Please check the provided stats_range parameter: {args['stats_range']}")
                    
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

def rc_smoothed_lims(sig, sig_er, sig_ray, sig_ray_er, ranges):
    
    mask = (sig == sig)
    
    first_range, last_range, range_span, range_edge = \
        get_span(ranges.copy().where(mask, drop=True).values)
            
    xlims_range = [
        1E-3 * (first_range - range_edge), 
        1E-3 * (last_range + range_edge)
        ]
    
    sig_sl = slice_signal_by_range(sig, ranges, region = [5000., ranges[-1].values])
    sig_ray_sl = slice_signal_by_range(sig_ray, ranges, region = [5000., ranges[-1].values])
    
    max_val = np.nanmax([np.nanmax(sig_sl), np.nanmax(sig_ray_sl)])
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
    
    
def _default_processor(y_tw, x_w, t):
    """
    Default stats processor for one (channel, bin_center) window.

    Parameters
    ----------
    y_tw : np.ndarray, shape (time, window)
        Signal values for this rolling window.
    x_w : np.ndarray, shape (window,)
        Range values (for this channel) corresponding to the same window bins.
    t : np.ndarray, shape (time,)
        Time axis in seconds since start.

    Returns
    -------
    tuple of scalars (mean, sdev, sem, vert_slope, temp_slope, gaussian_noise,
                      profiles, bins, points)
    """
    # If rolling window not fully available, xarray will feed NaNs -> return NaNs
    if np.isnan(y_tw).any() or np.isnan(x_w).any():
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                np.nan, np.nan, np.nan)

    profiles = float(t.size)
    bins = float(x_w.size)
    points = float(np.isfinite(y_tw).sum())

    mean = float(np.nanmean(y_tw))
    sdev = float(np.nanstd(y_tw, ddof=0))
    sem = float(sdev / np.sqrt(points)) if points > 0 else np.nan

    # Vertical trend: mean over time -> profile over window bins
    y_profile = np.nanmean(y_tw, axis=0)  # (window,)
    vert_p = linregress(x=x_w, y=y_profile).pvalue
    vert_slope = float(vert_p <= 0.05)

    # Temporal trend: mean over bins -> time series over time
    y_time = np.nanmean(y_tw, axis=1)  # (time,)
    temp_p = linregress(x=t, y=y_time).pvalue
    temp_slope = float(temp_p <= 0.05)

    # Gaussian noise test on demeaned profile residuals (matches your idea)
    # Shapiro requires at least 3 samples
    if y_profile.size >= 3:
        gaussian_noise = float(shapiro(y_profile - np.mean(y_profile)).pvalue > 0.05)
    else:
        gaussian_noise = np.nan

    return (mean, sdev, sem, vert_slope, temp_slope, gaussian_noise,
            profiles, bins, points)


def calculate_statistics_rolling(sig, ranges, window, region=None, processor=None, center=True):
    """
    Rolling-bin statistics with an easy-to-modify processor hook.

    Parameters
    ----------
    sig : xr.DataArray
        dims: (time, channel, bins)
    ranges : xr.DataArray
        dims: (channel, bins)
    window : int
        Rolling window size in *bins*.
    region : None or [lower_bin, upper_bin]
        Bin boundaries (inclusive) where rolling stats are computed.
        Outside region -> NaNs.
    processor : callable or None
        Function(y_tw, x_w, t) -> tuple of scalar stats.
        If None, uses _default_processor above.
    center : bool
        Center the rolling window on each bin (recommended).

    Returns
    -------
    xr.Dataset
        dims: (channel, bins)
        vars: mean, sdev, sem, vert_slope, temp_slope, gaussian_noise,
              profiles, bins_in_window, points
    """
    if processor is None:
        processor = _default_processor

    # time in seconds since start (same intent as your current code)
    dt = sig["time"] - sig["time"].isel(time=0)
    t = (dt.dt.seconds + 1e-6 * dt.dt.microseconds).astype(float).values  # (time,)

    # Select region in bin-index space (inclusive)
    n_bins = sig.sizes["bins"]
    if region is None:
        lo, hi = 0, n_bins - 1
    else:
        lo, hi = int(region[0]), int(region[1])
        lo = max(lo, 0)
        hi = min(hi, n_bins - 1)
        if hi < lo:
            raise ValueError("region must satisfy upper_bin >= lower_bin")

    sig_reg = sig.isel(bins=slice(lo, hi + 1))
    ranges_reg = ranges.isel(bins=slice(lo, hi + 1))

    # Build rolling windows (NaN edges when not fully covered)
    sig_win = (
        sig_reg.rolling(bins=window, center=center, min_periods=window)
        .construct("window")
    )  # dims: (time, channel, bins, window)

    ranges_win = (
        ranges_reg.rolling(bins=window, center=center, min_periods=window)
        .construct("window")
    )  # dims: (channel, bins, window)

def range_correction(sig, ranges):
    
    sig_rc = sig * np.power(ranges, 2)
    
    mask = (ranges > 0.).broadcast_like(sig_rc)

    sig_rc = sig_rc.where(mask,sig)
    
    return(sig_rc)
    
    
def smoothing_1D(sig, ranges, args):

    sig_sm = sig.copy()
    sig_er = np.nan * sig.copy()
    ranges_sm = 1E-3 * ranges.copy()
       
    window = args['smoothing_window_rc']
    half_window = np.ceil(window / 2)
    
    if args['smoothing_range_rc'] is None:
        args['smoothing_range_rc'] = [0.5, ranges_sm[-1].values]
            
    if args['smooth_rc']:

        sig_sm.values, sig_er.values = \
            smooth_1D(y_vals = sig_sm.values, 
                      x_vals = ranges_sm.values,
                      x_sm_lims = args['smoothing_range_rc'],
                      x_sm_win = args['smoothing_window_rc'],
                      expo = False,
                      err_type = 'std',
                      mode = 'range')
            
        range_mask = (ranges <= ranges[-1].values - half_window) 
        
        sig_sm = sig_sm.where(range_mask)
        sig_er = sig_er.where(range_mask)
    
    return(sig_sm, sig_er)

def smoothing_2D(sig, ranges, args, channel_info, mode):
    
    zero_bin = channel_info.loc['DAQ_Trigger_Offset']
    
    # if args['smoothing_range'] is None:
    #     if zero_bin <= 100:
    #         args['smoothing_range'] = [
    #             -int(channel_info.loc['DAQ_Trigger_Offset']) + 100,
    #             sig.bins.size - 1
    #             ]
    #     else:
    #         args['smoothing_range'] = [
    #             0,
    #             sig.bins.size - 1
    #             ]  
    
    # resol = channel_info.loc['Raw_Data_Range_Resolution']
    # smoothing_window = args['smoothing_window'] * resol
    window = args['smoothing_window']
    half_window = np.ceil(args['smoothing_window']/2)
    
    if args['smoothing_range'] is None:
        args['smoothing_range'] = [0,
            sig.bins.size - 1
            ]    
    
    sig_sm = sig.copy()
    sig_er = np.nan * sig.copy()
    ranges_sm = 1E-3 * ranges.copy()
    bins = sig.bins.copy()
    
    if args['smooth']:
        sig_sm.values, sig_er.values = sliding_average_2D_fast(
            z_vals = sig_sm.values, 
            y_vals = ranges_sm.values,
            y_sm_lims = args['smoothing_range'],
            y_sm_win = window,
            err_type = 'std',
            mode = mode
            )
        
        
        bin_mask = (bins >= half_window) & \
            (bins <= sig.bins.size - 1 - half_window) 
            
        sig_sm = sig_sm.where(bin_mask)
        sig_er = sig_er.where(bin_mask)
        
    return(sig_sm, sig_er)

def error_per_bin(sig, args, channel_info):

    sig_av = sig.mean('time')

    n_time = sig.time.size
    
    if args['background_range'] is None:
        args['background_range'] = [
            int(channel_info.loc['Background_Low_Bin']),
            int(channel_info.loc['Background_High_Bin'])
            ]

    if args['background_range'] is not None:
        
        background_llim = args['background_range'][0]#np.where(ranges.values >= args['background_range'][0])[0][0]
        background_ulim = args['background_range'][1]#np.where(ranges.values <= args['background_range'][1])[0][-1]
        
        bin_d = dict(bins = slice(background_llim,background_ulim))

        er_region = sig_av.loc[bin_d].std(dim = 'bins').values
        er_bin = er_region * np.sqrt(n_time)
    
    return(er_bin)

def set_background_range(args, channel_info):
    
    if args['background_range'] is None:
        args['background_range'] = [
            int(channel_info.loc['Background_Low_Bin']),
            int(channel_info.loc['Background_High_Bin'])
            ]

    if args['background_range'] is not None:
        
        background_llim = args['background_range'][0]#np.where(ranges.values >= args['background_range'][0])[0][0]
        background_ulim = args['background_range'][1]#np.where(ranges.values <= args['background_range'][1])[0][-1]
        args['background_range'] = [background_llim, background_ulim]
        
    
def background_correction(sig, background_range):
    
    sig_bgc = sig.copy()

    bin_d = dict(bins = slice(background_range[0], background_range[1]))

    sig_bgr = sig.loc[bin_d].mean(dim = 'bins')

    sig_bgc = sig_bgc - sig_bgr.copy()

    return(sig_bgc, sig_bgr)

def extract_xy_arrays(profiles, profiles_ray, ch_d):
    
    ch = ch_d['channel']
    
    sig_ch = profiles['sig'].loc[ch_d].copy()
    ranges_ch = profiles['ranges'].loc[ch_d].copy()
    
    if profiles_ray is not None:
        sig_ray_ch = profiles_ray['sig'].loc[ch_d].copy()
        ranges_ray_ch = profiles_ray['ranges'].loc[ch_d].copy()
        if not coords_are_subset(ranges_ray_ch, ranges_ch):
            raise Exception(f"--Error: The range is different for the dark and normal files for channel {ch} ")
        sig_ray_ch = sig_ray_ch + sig_ch
            
        index_ne = np.where((sig_ray_ch == sig_ray_ch) & (sig_ch == sig_ch))[0]
        sig_ch = sig_ch[index_ne]
        sig_ray_ch = sig_ray_ch[index_ne]
        ranges_ch = ranges_ch[index_ne]
        
    else:
        sig_ray_ch = np.nan * xr.zeros_like(sig_ch)
        ranges_ray_ch = ranges_ch.copy()
        
    return(ranges_ch, sig_ch, sig_ray_ch)

def coords_are_subset(da_small, da_big):
    for dim in da_small.dims:
        if dim not in da_big.dims:
            return False
        
        small_vals = da_small.coords[dim].values
        big_vals = da_big.coords[dim].values
        
        if not np.isin(small_vals, big_vals).all():
            return False

    return True

def resampling(sig_raw, averaging_rate, averaging_threshold):

    # Averaging the Rayleigh measurement for quicklooks
    if averaging_rate == None:
        mask_incomplete = xr.full_like(sig_raw, 
                                       fill_value = False, 
                                       dtype=bool)
    else:
        sig_avg_ray_qck, mask_incomplete = temporal_averaging(
            sig = sig_raw, 
            averaging_rate = averaging_rate, 
            averaging_threshold = averaging_threshold
            )
        
    return()

def temporal_averaging(sig: xr.DataArray, averaging_rate: str, 
                       averaging_threshold: float):
    
    delta_t_min = np.min(sig.time[1:].values-sig.time[:-1].values)
    
    expected_profiles = expected_profiles_per_window(averaging_rate, delta_t_min)
    
    if expected_profiles < 1:
        print(f"The provided averaging_rate ({averaging_rate}) is smaller than the temporal resolution. Averaging is not possible, the raw temporal resolution will be used")
        sig_avg = sig
    
    else:    
        
        origin = pd.Timestamp(sig.time.values[0])

        actual_profiles = sig.resample({"time" : averaging_rate},
                                       label = "left", 
                                       origin = origin).count()
        
        mask_incomplete = actual_profiles / expected_profiles < averaging_threshold
        
        sig_avg = sig.resample({"time" : averaging_rate},
                               label = "left", 
                               origin = origin).mean()
    
        sig_avg = sig_avg.where(~mask_incomplete, np.nan)
        
    return(sig_avg, mask_incomplete)

def expected_profiles_per_window(freq: str,
                                 sample_dt: np.timedelta64,
                                 include_end: bool = False) -> int:
    """
    Return the expected number of profiles within one resampling window.

    freq        : pandas offset alias (e.g. '3min', '1H', '90s', '1H30min')
    sample_dt   : np.timedelta64 (e.g. np.timedelta64(15_000_000_000, 'ns'))
    include_end : if True, count includes the right edge when it lands exactly on a sample
                  (i.e., closed='both' / right-closed windows)

    Assumes windows are left-closed, right-open by default (include_end=False),
    matching xarray/pandas resample defaults.
    """
    # Convert both to integer nanoseconds to avoid unit mismatches
    W_ns = pd.to_timedelta(freq).value                  # window length in ns (int)
    dt_ns = pd.to_timedelta(sample_dt).value            # sample spacing in ns (int)

    if W_ns <= 0 or dt_ns <= 0:
        raise ValueError("freq and sample_dt must be positive durations.")

    count = W_ns // dt_ns  # floor(W/dt) for [left-closed, right-open) bins

    # If you want to include the right edge when exactly aligned:
    if include_end and (W_ns % dt_ns == 0):
        count += 1

    return int(count)

def pass_to_args(args, data_list, data_keys):
    
    for i in range(len(data_keys)):
        args[data_keys[i]] = data_list[i]
        
    return(args)
    
def get_max_channel_height(norm_region, norm_region_flag):
    
    max_channel_height = np.mean(norm_region)
    max_channel_height = str(int(np.round(1E3*max_channel_height, decimals=-2)))
    
    return(max_channel_height)
    

def get_y_limits(y1_vals, y2_vals, y_lims, wavelength, use_lin_scale):
    
    # Get the max signal bin and value       
    y_max = np.nanmax(y1_vals)
    y_min = np.nanmin(y2_vals)

    scale_f = wavelength / 355.
    scat_ratio_f = 2.5
    
    # Get the signal axis upper limit
    if use_lin_scale == False:
        if y_lims[-1] == None:
            y_ulim = scat_ratio_f * scale_f * y_max
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_ulim = 1
            else:
                y_ulim =  y_lims[-1]
    else:
        if y_lims[-1] == None:
            y_ulim = scat_ratio_f * scale_f * y_max
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_ulim = 1
            else:
                y_ulim =  y_lims[-1]
        
    # Get the signal axis lower limit
    if use_lin_scale == False:
        if y_lims[0] == None:
            y_llim = y_min / 2.
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_llim = 0.
            else:
                y_llim =  y_lims[0]
    else:
        if y_lims[0] == None:
            y_llim = y_min / 2.
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_llim = 0.
            else:
                y_llim =  y_lims[0]
    
    y_lims = [y_llim, y_ulim]
    
    return(y_lims)

def x_unit_conversions(ranges):
   
    # Convert meters to kilometers and select ranges or heights for the x axis depending on the use_dis value 
    x_vals = 1E-3 * ranges

    return(x_vals)

def get_x_label(use_range):
   
    # Convert meters to kilometers and select ranges or heights for the x axis depending on the use_dis value 
    if use_range:
        x_label = "Range from the lidar [km]"
        
    else:
        x_label = "Height above the lidar [km]"

    return(x_label)

def y_unit_conversions(sig, sig_err, norm_coef):
    
    # Multiply the   
    y_vals  = norm_coef * sig.copy()
    
    y_errs = norm_coef * sig_err.copy()
    
    return(y_vals, y_errs)

def slice_arrays(x_lims, x_vals, y1_vals, y2_vals):
    
    x_mask = (x_vals >= x_lims[0]) & (x_vals <= x_lims[1])
    
    X = x_vals[x_mask]

    Y1  = y1_vals[x_mask]
    Y2  = y2_vals[x_mask]
    
    return(X, Y1, Y2)

def slice_2D_array(x_lims, x_vals, y_vals):
    
    x_mask = (x_vals >= x_lims[0]) & (x_vals <= x_lims[1])
    
    X = x_vals[x_mask]

    Y  = y_vals[:,x_mask]
    
    return(X, Y)

def add_extra_plot_metadata(plot_metadata, norm_region_flag, 
                            stats_norm_region, maximum_channel_height):
    
    plot_metadata['norm_region_flag'] = f"{norm_region_flag}"
    for key in stats_norm_region.keys():
        plot_metadata[f"stats_{key}"] = f"{stats_norm_region[key]}"
    for key in stats_norm_region.keys():
        plot_metadata[f"masks_{key}"] = f"{stats_norm_region[key]}"
        
    plot_metadata['maximum_channel_height'] = f"{maximum_channel_height}"
    
    return(plot_metadata)

if __name__ == '__main__':
    # Get the command line argument information
    args = call_parser()

    # Call main
    main(args)

def _parse_duration_to_timedelta64(freq: str) -> np.timedelta64:
    """
    Parse strings like '2min', '1h', '30s', '500ms', '1D' into np.timedelta64.
    Supported units: ns, us, ms, s, min, h, D
    """
    s = freq.strip()
    # split into leading integer + unit
    i = 0
    while i < len(s) and s[i].isdigit():
        i += 1
    if i == 0:
        raise ValueError(f"Duration must start with an integer, got: {freq!r}")

    n = int(s[:i])
    unit = s[i:].strip()

    unit_map = {
        "s": "s",
        "sec": "s",
        "secs": "s",
        "min": "m",
        "mins": "m",
        "m": "m",
        "h": "h",
        "hr": "h",
        "hrs": "h",
    }
    if unit not in unit_map:
        raise ValueError(
            f"Unsupported unit {unit!r} in {freq!r}. "
            "Use one of: s, min, h"
        )

    return np.timedelta64(n, unit_map[unit])

def infer_dt(time_values: np.ndarray, method: str = "median") -> np.timedelta64:
    """
    Infer a representative sampling interval from a 1D np.datetime64 array.
    Uses median (robust) by default.
    """
    t = np.asarray(time_values)
    if t.size < 2:
        raise ValueError("Need at least 2 time samples to infer resolution.")
    diffs = np.diff(t).astype("timedelta64[ns]").astype(np.int64)
    if method == "mean":
        dt_ns = int(np.mean(diffs))
    elif method == "mode":
        # simple mode-ish: most common diff
        vals, counts = np.unique(diffs, return_counts=True)
        dt_ns = int(vals[np.argmax(counts)])
    else:
        raise ValueError("method must be 'mean' or 'mode'")
    return np.timedelta64(dt_ns, "ns")

def block_mean_fixed_samples(
    da: xr.DataArray,
    N: int,
    dim: str = "time",
    label: str = "start",  # 'start' | 'middle' | 'end'
) -> xr.DataArray:
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    if da.sizes[dim] < N:
        raise ValueError(f"Not enough samples ({da.sizes[dim]}) for N={N}.")

    # Trim to exact multiple of N so every block has exactly N samples
    nblocks = da.sizes[dim] // N
    da_trim = da.isel({dim: slice(0, nblocks * N)})

    # Coarsen does exactly what we want
    out = da_trim.coarsen({dim: N}, boundary="exact").mean()

    # Assign a representative time coordinate
    t = da_trim[dim].values
    if label == "start":
        t_rep = t[::N]
    elif label == "middle":
        t_rep = t[(N // 2)::N]
    elif label == "end":
        t_rep = t[(N - 1)::N]
    else:
        raise ValueError("label must be 'start', 'middle', or 'end'")

    out = out.assign_coords({dim: t_rep[: out.sizes[dim]]})
    return out

def block_mean_by_duration(
    da: xr.DataArray,
    duration: str,               # e.g. "2min", "1h"
    dim: str = "time",
    dt_method: str = "mean",   # 'mean' or 'mode'
    label: str = "start",
) -> xr.DataArray:
    target = _parse_duration_to_timedelta64(duration)
    dt = infer_dt(da[dim].values, method=dt_method)

    # Compute N = target / dt (rounded to nearest)
    target_ns = target.astype("timedelta64[ns]").astype(np.int64)
    dt_ns = dt.astype("timedelta64[ns]").astype(np.int64)
    if dt_ns <= 0:
        raise ValueError("Inferred dt is non-positive; check time coordinate ordering.")

    N = int(np.round(target_ns / dt_ns))
    N = max(N, 1)

    return (block_mean_fixed_samples(da, N=N, dim=dim, label=label), N)

def round_it(x, sig):
    
    if not np.isfinite(x) or np.isnan(x):
        x = -999.
        sig = 3
        
    if x != 0.:
        x_out = np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)
    else:
        x_out = 0.
        
    return x_out

def make_colorscale(colorscale, n_time):
    
    if colorscale == 'discrete':
        colors = int(np.ceil(n_time / 10)) * Category10[10]
    elif colorscale == 'sequential':
        if n_time <= 256:
            colors = turbo(n_time)
        else:
            max_colors = turbo(256)
            indexes = np.arange(0, n_time)
            color_index = (np.round(255 * indexes / n_time, decimals = 0)).astype(int)
            colors = tuple([max_colors[ind] for ind in color_index])
                
    return(colors)