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
from .plotting import make_axis, make_title, make_plot
from .writters import make_header, export_ascii
from .tools import curve_fit
import os
       
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
        sig_ch = profiles['sig'].loc[ch_d].values
        atb_ch = profiles['atb'].loc[ch_d].values

        
        ranges_ch = profiles['ranges'].copy().loc[ch_d].values
        height_ch = profiles['heights'].copy().loc[ch_d].values
        
        # Create the y axis (height/range)
        x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
            make_axis.rayleigh_x(heights = height_ch, 
                                 ranges = ranges_ch,
                                 x_lims = args['x_lims'], 
                                 use_dis = args['use_range'])
    
        # Smoothing
        if args['smooth']:
            if not isinstance(args['smoothing_window'],list):
                from .tools.smoothing import sliding_average_1D_fast as smooth_1D
            else:
                from .tools.smoothing import sliding_average_1D as smooth_1D

            y_vals_sm, y_errs = \
                smooth_1D(y_vals = sig_ch.copy(), 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')

            atb_ch, y_errs_mol = \
                smooth_1D(y_vals = atb_ch.copy(), 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')
        
        else:
            y_vals_sm = sig_ch.copy()
            y_errs = np.nan * y_vals_sm   
        
        ulim = 16.
        min_win = 1.
        max_win = 4.
        rsem_lim = 0.02
        
        if args['normalization_region'][0] < args['cross_check_lim']:
            raise Exception(f"The provided cross_check_lim {args['cross_check_lim']} is lower than the lower edge of the normalization region {args['normalization_region'][0]}. Please revise the rayleight_fit section of the settings file")
        if args['normalization_region'][0] > ulim:
            ulim = args['normalization_region'][1]
        if args['normalization_region'][1] - args['normalization_region'][0] < min_win:
            min_win = args['normalization_region'][1] - args['normalization_region'][0]
        if args['normalization_region'][1] - args['normalization_region'][0] > max_win:
            max_win = args['normalization_region'][1] - args['normalization_region'][0]
        
        # Check for a fit range
        rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef = \
            curve_fit.stats(y1 = sig_ch.copy(),
                            y2 = atb_ch.copy(), 
                            x  = x_vals,
                            min_win = 1.,
                            max_win = 4.,
                            step = 0.1,
                            llim = args['cross_check_lim'],
                            ulim = 16.,
                            rsem_lim = rsem_lim)   
            
        norm_region, idx, fit = \
            curve_fit.scan(mfit = mfit,
                           dflt_region = args['normalization_region'],
                           auto_fit = args['auto_fit'])
                
        nder_c = float(nder[idx].values)
        mder_c = float(mder[idx].values)
        rsem_c = float(rsem[idx].values)
        coef_c = float(coef[idx].values)
        
        # rsem_c = coef_c[0] * sem_c / atb_c
        
        # Create the y axis (signal)
        y_llim, y_ulim, y_label = \
            make_axis.rayleigh_y(sig = coef_c * y_vals_sm[slice(x_lbin,x_ubin+1)], 
                                 atb = atb_ch.copy()[slice(x_lbin,x_ubin+1)], 
                                 y_lims = args['y_lims'],
                                 wave = metadata['dwl'].copy().loc[ch].values,
                                 use_lin = args['use_lin_scale'])
        
        # Make title
        title = make_title.rayleigh(channel = ch, 
                                    metadata = metadata, 
                                    args = args)
        
        # Make plot filename
        fname = make_plot.make_filename(metadata = metadata, 
                                        channel = ch, 
                                        meas_type = 'ray', 
                                        version = __version__)
        
        # Make the png file
        plot_path = make_plot.rayleigh(dir_out = os.path.join(args['output_folder'],'plots'), 
                                       fname = f"{fname}.png", title = title,
                                       dpi_val = args['dpi'],
                                       color_reduction = args['color_reduction'],
                                       use_lin = args['use_lin_scale'],
                                       norm_region = norm_region,
                                       x_vals = x_vals, 
                                       y1_vals = y_vals_sm,
                                       y2_vals = atb_ch,
                                       y1_errs = y_errs,
                                       y2_errs = np.nan * np.zeros(y_errs.shape),
                                       coef = coef_c,
                                       rsem = rsem_c,
                                       rslope = nder_c,
                                       pval = mder_c,
                                       rsem_lim = rsem_lim,
                                       fit = fit,
                                       auto_fit = args['auto_fit'],
                                       x_lbin = x_lbin, x_ubin = x_ubin,
                                       x_llim = x_llim, x_ulim = x_ulim, 
                                       y_llim = y_llim, y_ulim = y_ulim, 
                                       x_label = x_label, y_label = y_label,
                                       x_tick = args['x_tick'],
                                       label_1 = 'measured',
                                       label_2 = 'molecular') 
        
        # Make ascii file header
        header = \
            make_header.rayleigh(channel = ch, 
                                 metadata = metadata, 
                                 norm_region = norm_region)

        # Export to ascii (Volker's format)        
        export_ascii.rayleigh(dir_out = args['output_folder'], 
                              fname = f"{fname}.txt", 
                              alt = x_vals, 
                              atb = atb_ch, 
                              rcs = sig_ch, 
                              header = header)
        
        if args['auto_fit']:
            
            fname_mask = make_plot.make_filename(metadata = metadata, 
                                                 channel = ch, 
                                                 meas_type = 'ray', 
                                                 extra_type = 'mask',
                                                 version = __version__)
        
            # Make title
            title_mask = make_title.rayleigh(channel = ch, 
                                             metadata = metadata, 
                                             args = args,
                                             is_mask = True)
            
            make_plot.rayleigh_mask(dir_out = os.path.join(args['output_folder'],'plots'), 
                                    fname = f"{fname_mask}.png", title = title_mask,
                                    dpi_val = args['dpi'],
                                    color_reduction = args['color_reduction'],
                                    mfit = mfit, mder = mder, msec = msec, 
                                    mshp = mshp, mcrc = mcrc, rsem = rsem,
                                    rsem_lim = rsem_lim)
                
        # Add metadata to the plots
        plot_metadata = make_plot.get_plot_metadata(metadata = metadata, 
                                                    args = args, 
                                                    channel = ch,
                                                    meas_type = 'ray', 
                                                    version = __version__)
        
        plot_metadata['fit'] = f"{fit}"
        plot_metadata['rsem'] = f"{rsem_c}"
        plot_metadata['rslope'] = f"{nder_c}"
        plot_metadata['maximum_channel_height'] = f"{norm_region[-1]}"
        
        make_plot.add_plot_metadata(plot_path = plot_path, 
                                    plot_metadata = plot_metadata)
        
    print('-----------------------------------------')
    print(' ')
    return()

if __name__ == '__main__':
    # Get the command line argument information
    args = call_parser()

    # Call main
    main(args)
