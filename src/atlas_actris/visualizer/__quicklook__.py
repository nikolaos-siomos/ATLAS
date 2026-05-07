#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import warnings, os, sys
from .readers.parse_qck_args import call_parser, check_parser
from .readers.check import check_channels
from .readers.read_prepro import unpack
from .plotting import make_axis, make_title, make_plot
       
# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__, meas_type):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print(f'Start generating Quicklooks ({meas_type})...')
    print('-----------------------------------------')

    profiles, metadata = unpack(args['input_file'])
    
    # # Extract signal time, channels, and bins
    time = profiles['sig'].time.values
    delta_t = metadata['Raw_Data_Stop_Time'] - metadata['Raw_Data_Start_Time']
    
    # Check if the parsed channels exist
    channels = \
        check_channels(sel_channels = args['channels'], 
                       all_channels = metadata['atlas_channel_id'],
                       exclude_telescope_type = args['exclude_telescope_type'], 
                       exclude_channel_type = args['exclude_channel_type'], 
                       exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                       exclude_channel_subtype = args['exclude_channel_subtype'])

    # Create the x axis (time)
    x_lbin, x_ubin, x_tick, t_vals, t_tick, = \
        make_axis.quicklook_x(t_lims = args['t_lims'],
                              t_tick = args['t_tick'], 
                              time = time)

    # iterate over the channels
    for ch in channels:
        
        print(f"-- channel: {ch}")

        ch_d = dict(channel = ch)
        
        sig_ch = profiles['sig'].loc[ch_d].values
        
        ranges_ch = profiles['ranges'].copy().loc[ch_d].values
        heights_ch = profiles['heights'].copy().loc[ch_d].values
            
        # Create the y axis (height/range)
        y_lbin, y_ubin, y_llim, y_ulim, y_vals, y_label = \
            make_axis.quicklook_y(heights = heights_ch, 
                                  ranges = ranges_ch,  
                                  y_lims = args['y_lims'], 
                                  use_dis = args['use_range'])

        # Smoothing
        if args['smooth']:
            if not isinstance(args['smoothing_window'],list):
                from .tools.smoothing import sliding_average_2D_fast as smooth_2D
            else:
                from .tools.smoothing import sliding_average_2D as smooth_2D

            z_vals, _ = smooth_2D(z_vals = sig_ch, 
                                  y_vals = y_vals,
                                  y_sm_lims = args['smoothing_range'],
                                  y_sm_win = args['smoothing_window'],
                                  expo = args['smooth_exponential'])
        else:
            z_vals = sig_ch

        # Create the z axis (signal)
        z_llim, z_ulim, z_vals = \
            make_axis.quicklook_z(z_vals = z_vals, 
                                  y_vals = y_vals,
                                  z_lims = args['z_lims'] , 
                                  use_log = args['use_log_scale'],
                                  z_max_zone = args['z_max_zone'],
                                  z_min_zone = args['z_min_zone'])
                
        # Make title
        title = make_title.quicklook(channel = ch, 
                                     metadata = metadata, 
                                     args = args)
  
        # Make filename
        fname = make_plot.make_filename(metadata = metadata, 
                                        channel = ch, 
                                        meas_type = f'qck_{meas_type}', 
                                        version = __version__)
        # Make the plot
        plot_path = make_plot.quicklook(dir_out = os.path.join(args['output_folder'],'plots'), 
                                        fname = f"{fname}.png", title = title,
                                        dpi_val = args['dpi'],
                                        color_reduction = args['color_reduction'],
                                        use_log = args['use_log_scale'],
                                        delta_t = delta_t,
                                        t_vals = t_vals, y_vals = y_vals, 
                                        z_vals = z_vals, 
                                        x_lbin = x_lbin, x_ubin = x_ubin, 
                                        y_lbin = y_lbin, y_ubin = y_ubin, 
                                        y_llim = y_llim, y_ulim = y_ulim, 
                                        z_llim = z_llim, z_ulim = z_ulim,
                                        y_label = y_label, 
                                        t_tick = t_tick, x_tick = x_tick,
                                        y_tick = args['y_tick'])  
    
     
        # Add metadata to the quicklook plot
        plot_metadata = make_plot.get_plot_metadata(metadata = metadata, 
                                                    args = args, 
                                                    channel = ch,
                                                    meas_type = f'qck_{meas_type}', 
                                                    version = __version__)
        
        make_plot.add_plot_metadata(plot_path = plot_path, 
                                    plot_metadata = plot_metadata)

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
    