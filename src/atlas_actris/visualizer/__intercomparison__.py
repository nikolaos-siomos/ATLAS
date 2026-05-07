#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, sys
import numpy as np
from .readers.parse_cmp_args import call_parser, check_parser
from .readers.check import check_channels_no_exclude as check_channels
from .readers.read_prepro import unpack
from .plotting import make_axis, make_title, make_plot
from .tools import curve_fit

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print('Initializing the Lidar Intercomparison...')
    print('-----------------------------------------')
    
    profiles_1, metadata_1 = unpack(args['input_file'][0])
    
    channels_1 = metadata_1['channels']

    profiles_2, metadata_2 = unpack(args['input_file'][1])
    
    channels_2 = metadata_2['channels']
    
    # Relative Stnadard error of the mean upper limit for the auto_fit
    rsem_lim = 0.02
    
    # Check if the parsed channels exist
    channels_1 = check_channels(sel_channels = args['channels_1'], 
                                all_channels = channels_1)
    channels_2 = check_channels(sel_channels = args['channels_2'], 
                                all_channels = channels_2)
    
    # iterate over the channels
    for i in range(channels_1.size):
   
        ch_1 = channels_1[i]
        ch_2 = channels_2[i]

        print(f"-- ch: {ch_1} & {ch_2}")
        
        ch_1_d = dict(channel = ch_1)
        ch_2_d = dict(channel = ch_2)

        sig_1_ch = profiles_1['sig'].copy().loc[ch_1_d].values
        atb_1_ch = profiles_1['atb'].copy().loc[ch_1_d].values
        
        # ranges_1_ch = profiles_1['ranges_ray'].copy().loc[ch_1_d].values
        heights_1_ch = profiles_1['heights_ray'].copy().loc[ch_1_d].values
        
        sig_2_ch = profiles_2['sig'].copy().loc[ch_2_d].values
        atb_2_ch = profiles_2['atb'].copy().loc[ch_2_d].values
                
        # ranges_2_ch = profiles_2['ranges_cal'].copy().loc[ch_2_d].values
        heights_2_ch = profiles_2['heights_cal'].copy().loc[ch_2_d].values
        
        min_height = np.max(heights_1_ch[0],heights_2_ch[0])
        max_height = np.min(heights_1_ch[-1],heights_2_ch[-1])
        
        mask_heights_1 = (heights_1_ch >= min_height) & (heights_1_ch <= max_height)
        mask_heights_2 = (heights_2_ch >= min_height) & (heights_2_ch <= max_height)
        
        sig_1_ch = sig_1_ch[mask_heights_1]
        atb_1_ch = atb_1_ch[mask_heights_1]
        heights_1_ch = heights_1_ch[mask_heights_1]
        
        sig_2_ch = sig_2_ch[mask_heights_2]
        atb_2_ch = atb_2_ch[mask_heights_2]
        heights_2_ch = heights_2_ch[mask_heights_2]
                        
        if heights_1_ch == heights_2_ch:
            heights_ch = heights_2_ch
        else:
            com_bins, com_ind_1, com_ind_2 = np.intersect1d(heights_1_ch, heights_2_ch, return_indices = True)
            if len(com_bins) <= 0.5 * np.min([len(heights_1_ch), len(heights_2_ch)]):
                step_1 = heights_1_ch[1] - heights_1_ch[0]
                step_2 = heights_2_ch[1] - heights_2_ch[0]
                max_step = np.max([step_1, step_2])
                
                heights_ch = np.arange(min_height + max_step, max_height, 2. * max_step)
                heights_llim_ch = np.arange(min_height, max_height - max_step, 2. * max_step)
                heights_ulim_ch = np.arange(min_height + 2. * max_step, max_height + max_step, 2. * max_step)
                
                sig_1_ch_com = np.nan * np.zeros(heights_ch.size)
                sig_2_ch_com = np.nan * np.zeros(heights_ch.size)
                atb_1_ch_com = np.nan * np.zeros(heights_ch.size)
                atb_2_ch_com = np.nan * np.zeros(heights_ch.size)

                for j in range(heights_ch.size):
                    mask_1_com = (heights_1_ch >= heights_llim_ch[j]) & (heights_1_ch < heights_ulim_ch[j])
                    mask_2_com = (heights_2_ch >= heights_llim_ch[j]) & (heights_2_ch < heights_ulim_ch[j])
                    sig_1_ch_com[j] = np.mean(sig_1_ch[mask_1_com])
                    sig_2_ch_com[j] = np.mean(sig_2_ch[mask_2_com])
                    atb_1_ch_com[j] = np.mean(atb_1_ch[mask_1_com])
                    atb_2_ch_com[j] = np.mean(atb_2_ch[mask_2_com])
                
                sig_1_ch = sig_1_ch_com
                sig_2_ch = sig_2_ch_com
                atb_1_ch = atb_1_ch_com
                atb_2_ch = atb_2_ch_com
                
            else:
                heights_ch = heights_2_ch[com_ind_2]
                sig_1_ch = sig_1_ch[com_ind_1]
                sig_2_ch = sig_2_ch[com_ind_2]
                atb_1_ch = atb_1_ch[com_ind_1]
                atb_2_ch = atb_2_ch[com_ind_2]
                

        # Create the y axis (height/range)
        x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
            make_axis.rayleigh_x(heights = heights_ch, 
                                 ranges = np.nan * heights_ch,
                                 x_lims = args['x_lims'], 
                                 use_dis = False)
    
        # Smoothing
        if args['smooth']:
            if not isinstance(args['smoothing_window'],list):
                from .tools.smoothing import sliding_average_1D_fast as smooth_1D
            else:
                from .tools.smoothing import sliding_average_1D as smooth_1D

            y1_vals_sm, y1_errs = \
                smooth_1D(y_vals = sig_1_ch, 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')
                
            y2_vals_sm, y2_errs = \
                smooth_1D(y_vals = sig_2_ch.copy(), 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')      
        
        else:
            y1_vals_sm = sig_1_ch.copy()
            y2_vals_sm = sig_2_ch.copy()
            y1_errs = np.nan * y1_vals_sm
            y2_errs = np.nan * y2_vals_sm
    
        rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef = \
            curve_fit.statistics(y1 = sig_1_ch.copy(),
                                 y2 = sig_2_ch.copy(), 
                                 x  = x_vals,
                                 min_win = 1.,
                                 max_win = 4.,
                                 step = 0.1,
                                 llim = 2.,
                                 ulim = 16.,
                                 rsem_lim = rsem_lim,
                                 cross_check_all_points = False,
                                 cross_check_type = 'all')

            
        # Negate the cross-check mask for the signal comparison
        # mneg = np.ones(mder.shape, dtype = bool)
        
        norm_region, idx, ismol = \
            curve_fit.scan(mfit = mfit,
                           dflt_region = args['normalization_region'],
                           auto_fit = args['auto_fit'])
                
        coef_c = float(coef[idx].values)
        rsem_c = float(rsem[idx].values)
        nder_c = float(nder[idx].values)
        mder_c = float(mder[idx].values)
    
        # Create the y axis (signal)
        y_llim, y_ulim, y_label = \
            make_axis.intercomparison_y(sig1 = coef_c * y1_vals_sm[slice(x_lbin,x_ubin+1)], 
                                        sig2 = y2_vals_sm[slice(x_lbin,x_ubin+1)], 
                                        y_lims = args['y_lims'],
                                        use_lin = args['use_lin_scale'])
            
                
        # Make title
        
        title = make_title.intercomparison(channel_1 = ch_1, 
                                           channel_2 = ch_2, 
                                           metadata_1 = metadata_1, 
                                           metadata_2 = metadata_2, 
                                           args = args)
            
        # Make filename
        fname = make_plot.make_filename_intercomparison(metadata_1 = metadata_1, 
                                                        metadata_2 = metadata_2, 
                                                        channel_1 = ch_1, 
                                                        channel_2 = ch_2, 
                                                        version = __version__)
        # fname = f'{meas_id_1}_{meas_id_2}_{lidar_name_1}_{lidar_name_2}_cmp_{channels_1[j]}_{channels_2[j]}_ATLAS_{__version__}.png'

        # Make the plot
        plot_path = make_plot.intercomparison(dir_out = args['output_folder'], 
                                           fname = f"{fname}.png", 
                                           title = title,
                                           dpi_val = args['dpi'],
                                           color_reduction = args['color_reduction'],
                                           use_lin = args['use_lin_scale'],
                                           norm_region = norm_region,
                                           x_vals = x_vals, 
                                           y1_vals = y1_vals_sm,
                                           y2_vals = y2_vals_sm,
                                           y1_errs = y1_errs,
                                           y2_errs = y2_errs,
                                           coef = coef_c,
                                           rsem = rsem_c,
                                           rslope = nder_c,
                                           pval = mder_c,
                                           rsem_lim = rsem_lim,
                                           x_lbin = x_lbin, x_ubin = x_ubin,
                                           x_llim = x_llim, x_ulim = x_ulim, 
                                           y_llim = y_llim, y_ulim = y_ulim, 
                                           x_label = x_label, y_label = y_label,
                                           x_tick = args['x_tick'],
                                           label_1 = ch_1,
                                           label_2 = ch_2)  
        
        if args['auto_fit']:
            fname_mask = make_plot.make_filename_intercomparison(metadata_1 = metadata_1, 
                                                                 metadata_2 = metadata_2, 
                                                                 channel_1 = ch_1, 
                                                                 channel_2 = ch_2, 
                                                                 version = __version__,
                                                                 extra_type = 'mask')
    
            make_plot.rayleigh_mask(dir_out = args['output_folder'], 
                                    fname = fname_mask, title = title,
                                    dpi_val = args['dpi'],
                                    color_reduction = args['color_reduction'],
                                    mfit = mfit, mder = mder, msec = msec, 
                                    mshp = mshp, mcrc = mcrc, rsem = rsem,
                                    rsem_lim = rsem_lim)
        
        # Add metadata to the quicklook plot
        plot_metadata_1 = get_plot_metadata(metadata = metadata_1, 
                                            args = args, 
                                            channel = ch_1,
                                            meas_type = 'cmp', 
                                            version = __version__,
                                            data_source_id = '1')
       
        plot_metadata_2 = get_plot_metadata(metadata = metadata_2, 
                                            args = args, 
                                            channel = ch_2,
                                            meas_type = 'cmp', 
                                            version = __version__,
                                            data_source_id = '2')
        
        add_plot_metadata(plot_path = plot_path, 
                          plot_metadata = plot_metadata_1,
                          plot_metadata_extra = plot_metadata_2)
        

    print('-----------------------------------------')
    print(' ')
    return()

if __name__ == '__main__':
    
    sys.path.append('../')
    
    from .version import __version__
    
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args, __version__)