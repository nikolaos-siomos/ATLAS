#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings
import numpy as np
from .readers.parse_ray_args import call_parser, check_parser
from .readers.check import check_channels
from .readers.read_prepro import unpack
from .plotting import make_title, plot_utils, plot_rayleigh, plot_rayleigh_mask
from .writters import make_header, export_ascii
from .tools import curve_fit
       
# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
    # Check the command line argument information
    args = check_parser(args)

    print('-----------------------------------------')
    print('Initializing the Rayleigh Fit...')
    print('-----------------------------------------')
    
    profiles, metadata = unpack(args['input_file'])
        
    # Check if the parsed channels exist
    channels = \
        check_channels(sel_channels = args['channels'], 
                       all_channels = metadata['atlas_channel_id'],
                       exclude_telescope_type = args['exclude_telescope_type'], 
                       exclude_channel_type = args['exclude_channel_type'], 
                       exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                       exclude_channel_subtype = args['exclude_channel_subtype'])

    # iterate over the channels
    for ch in channels:
        print(f"-- channel: {ch}")

        ch_d = dict(channel = ch)
        
        args_ch = args.copy()

        sig_ch = profiles['sig'].loc[ch_d].copy().values
        atb_ch = profiles['atb'].loc[ch_d].copy().values
        
        ranges_ch = profiles['ranges'].loc[ch_d].copy().values
        height_ch = profiles['heights'].loc[ch_d].copy().values

        wave_ch = metadata['dwl'].loc[ch_d].copy().values
        
        # Convert the range/height units to km 
        x_vals = x_unit_conversions(ranges = ranges_ch, 
                                    heights = height_ch, 
                                    use_range = args_ch['use_range'])

        # Get the x axis label depending on use_dis
        x_label = get_x_label(args_ch['use_range'])

        # Trim the x and y related arrays using the x axis limits
        x_vals, y1_vals, y2_vals  = \
            slice_arrays(x_lims = args_ch['x_lims'], 
                         x_vals = x_vals, 
                         y1_vals = sig_ch, 
                         y2_vals = atb_ch)

        # Smoothing of the y1 array - generates also the corresponding standard deviation
        y1_vals_sm, y1_errs = smoothing(args = args_ch, 
                                        x_vals = x_vals, 
                                        y_vals = y1_vals, 
                                        err_type = "std")
            
        # Smoothing of the y2 array 
        y2_vals_sm, y1_sems = smoothing(args = args_ch, 
                                        x_vals = x_vals, 
                                        y_vals = y2_vals, 
                                        err_type = "sem")
        
        # Check for a Rayleigh fit range
        stats, masks = \
            curve_fit.statistics(
                y1 = y1_vals,
                y2 = y2_vals, 
                y1_err = y1_errs,
                y1_avg = y1_vals_sm,
                x  = x_vals,
                keyw_args = args_ch
                )


        # Identify the uppermost range, bin index and related stats there
        norm_region, norm_region_flag, idx = \
            curve_fit.scan(
                masks, 
                user_norm_region = args_ch['normalization_region'],
                auto_fit = args_ch['auto_fit'],
                prefered_range = 'far')
            

        # Isolate the values of the statistics and masks on the normalization region
        stats_norm_region, masks_norm_region = \
            curve_fit.metrics_norm_region(idx = idx, 
                                          stats = stats, 
                                          masks = masks)
        
        # Get the maximum channel height
        maximum_channel_height = \
            get_max_channel_height(norm_region= norm_region, 
                                   norm_region_flag = norm_region_flag)

        # Normalize y1_vals and y1_errs with the normalization factor from the Rayleigh fit test
        y1_vals_nrm, y1_errs_nrm = y_unit_conversions(
            sig = y1_vals_sm, 
            sig_err = y1_errs, 
            norm_coef = stats_norm_region['normalization_factor']
            )

        # Get the limits of the y axis in case they are not provided by the user (default)
        y_lims = get_y_limits(y1_vals = y1_vals_nrm, 
                              y2_vals = y2_vals_sm, 
                              y_lims = args_ch['y_lims'], 
                              wavelength = wave_ch, 
                              use_lin_scale = args_ch['use_lin_scale'])
    
        # Make title
        title = make_title.rayleigh(channel = ch, 
                                    metadata = metadata, 
                                    args = args_ch)
        
        # Make plot filename
        fname = plot_utils.make_filename(metadata = metadata, 
                                        channel = ch, 
                                        meas_type = 'ray', 
                                        version = __version__)
        
        # Pass all generated scalar or list parameters relevant to the plots to the args dictionary
        args_ch = pass_to_args(args = args_ch,
                     data_list = [stats_norm_region, 
                                  masks_norm_region, 
                                  maximum_channel_height, 
                                  norm_region, 
                                  norm_region_flag,
                                  x_label,
                                  y_lims,
                                  title,
                                  fname],
                     data_keys = ['stats_norm_region''',
                                  'masks_norm_region', 
                                  'maximum_channel_height', 
                                  'normalization_region', 
                                  'normalization_flag',
                                  'x_label',
                                  'y_lims',
                                  'title',
                                  'fname'])
        
        # Generate the Rayleigh fit plot
        ray_plot_path = plot_rayleigh.generate_plot(X = x_vals, 
                                                    Y1 = y1_vals_nrm,
                                                    Y2 = y2_vals_sm,
                                                    Y1E = y1_errs_nrm,
                                                    args = args_ch) 

        # Perform color reduction        
        plot_utils.perform_color_reduction(color_reduction = args_ch['color_reduction'], 
                                           plot_path = ray_plot_path)

        # Make ascii file header
        header = \
            make_header.rayleigh(channel = ch, 
                                 metadata = metadata, 
                                 norm_region = norm_region)

        # Export to ascii (Volker's format)        
        export_ascii.rayleigh(dir_out = args_ch['output_folder'], 
                              fname = f"{fname}.txt", 
                              alt = x_vals, 
                              atb = y2_vals, 
                              rcs = y1_vals, 
                              header = header)
        
        
        # Make mask plot filename
        fname_mask = plot_utils.make_filename(metadata = metadata, 
                                              channel = ch, 
                                              meas_type = 'ray', 
                                              extra_type = 'mask',
                                              version = __version__)
    
        # Make title
        title_mask = make_title.rayleigh(channel = ch, 
                                         metadata = metadata, 
                                         args = args_ch,
                                         is_mask = True)
        # Pass all additional scalar or list parameters relevant to the mask plots to the args dictionary
        args_mask_ch = args_ch.copy()
        
        args_mask_ch = pass_to_args(args = args_mask_ch,
                     data_list = [title_mask,
                                  fname_mask],
                     data_keys = ['title',
                                  'fname'])
        
        # Generate the molecular mask plot
        ray_mask_plot_path = \
            plot_rayleigh_mask.generate_plot(masks = masks,
                                             args = args_mask_ch)

        # Perform color reduction        
        plot_utils.perform_color_reduction(
            color_reduction = args_mask_ch['color_reduction'], 
            plot_path = ray_mask_plot_path
            )

            
        # Gather the metadata that are common for all QA tests in a dictonary
        plot_metadata = plot_utils.get_plot_metadata(metadata = metadata, 
                                                     args = args_mask_ch, 
                                                     channel = ch,
                                                     meas_type = 'ray', 
                                                     version = __version__)
        
        # Add additional metada that are special to the Rayleigh fit test in the dictionary
        plot_metadata = \
            add_extra_plot_metadata(plot_metadata = plot_metadata, 
                                    norm_region_flag = norm_region_flag, 
                                    stats_norm_region = stats_norm_region, 
                                    maximum_channel_height = maximum_channel_height)
        
        # Add the metadata to the Rayleigh fit plot  
        plot_utils.add_plot_metadata(plot_path = ray_plot_path, 
                                     plot_metadata = plot_metadata)

        # Add the metadata to the molecular mask plot 
        plot_utils.add_plot_metadata(plot_path = ray_plot_path, 
                                     plot_metadata = plot_metadata)
        

def pass_to_args(args, data_list, data_keys):
    
    for i in range(len(data_keys)):
        args[data_keys[i]] = data_list[i]
        
    return(args)
    
def get_max_channel_height(norm_region, norm_region_flag):
    
    max_channel_height = np.mean(norm_region)
    max_channel_height = str(int(np.round(1E3*max_channel_height, decimals=-2)))
    
    return(max_channel_height)
    
def smoothing(args, x_vals, y_vals, err_type = "std"):
    
    if args['smooth']:
        if not isinstance(args['smoothing_window'],list):
            # from .tools.smoothing import sliding_average_1D_fast as smooth_1D
            from .tools.smoothing import sliding_average_1D_fast as smooth_1D

        else:
            # from .tools.smoothing import sliding_average_1D as smooth_1D
            from .tools.smoothing import sliding_average_1D as smooth_1D


        y_vals_sm, y_errs = \
            smooth_1D(y_vals = y_vals, 
                      x_vals = x_vals,
                      x_sm_lims = args['smoothing_range'],
                      x_sm_win = args['smoothing_window'],
                      expo = args['smooth_exponential'],
                      err_type = err_type)
    
    else:
        y_vals_sm = y_vals.copy()
        y_errs = np.nan * y_vals.copy()
        
    return(y_vals_sm, y_errs)

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

def x_unit_conversions(ranges, heights, use_range):
   
    # Convert meters to kilometers and select ranges or heights for the x axis depending on the use_dis value 
    if use_range:
        x_vals = 1E-3 * ranges
        
    else:
        x_vals = 1E-3 * heights      

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
