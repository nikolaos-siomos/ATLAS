#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import warnings, os, sys, glob
import xarray as xr
import numpy as np
import netCDF4 as nc
from .readers.parse_qck_args import call_parser, check_parser
from .readers.check import check_channels
from .plotting import make_axis, make_title, make_plot

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__, meas_type):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print('Start generating Quicklooks...')
    print('-----------------------------------------')

    # Read the quicklook file
    data = xr.open_dataset(args['input_file'])
    
    # Extract signal
    sig = data.Range_Corrected_Signals
    sig = sig.copy().where(sig != nc.default_fillvals['f8'])

    # Extract signal time, channels, and bins
    time = sig.time.values
    bins = sig.bins.values

    # Check if the parsed channels exist
    channels = \
        check_channels(sel_channels = args['channels'], 
                       all_channels = data.channel.values,
                       exclude_telescope_type = args['exclude_telescope_type'], 
                       exclude_channel_type = args['exclude_channel_type'], 
                       exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                       exclude_channel_subtype = args['exclude_channel_subtype'])

    # Create the x axis (time)
    x_lbin, x_ubin, x_vals, t_vals, x_label, t_label, x_tick, t_tick, nodes, = \
        make_axis.quicklook_x(x_lims = args['x_lims'], 
                              x_tick = args['x_tick'],
                              t_tick = args['t_tick'], 
                              time = time)

    # iterate over the channels
    for ch in channels:
        
        print(f"-- channel: {ch}")

        ch_d = dict(channel = ch)
        sig_ch = sig.copy().loc[ch_d]
            
        # Create the y axis (height/range)
        y_lbin, y_ubin, y_llim, y_ulim, y_vals, y_label = \
            make_axis.quicklook_y(heights = data.Height_levels.loc[ch_d].values, 
                                  ranges = data.Range_levels.loc[ch_d].values,
                                  y_lims = args['y_lims'], 
                                  use_dis = args['use_distance'])

        # Smoothing
        if args['smooth']:
            if not isinstance(args['smoothing_window'],list):
                from .tools.smoothing import sliding_average_2D_fast as smooth_2D
            else:
                from .tools.smoothing import sliding_average_2D as smooth_2D

            z_vals = smooth_2D(z_vals = sig_ch, 
                               y_vals = y_vals,
                               y_sm_lims = args['smoothing_range'],
                               y_sm_win = args['smoothing_window'],
                               expo = args['smooth_exponential'])
        else:
            z_vals = sig_ch

        # Create the z axis (signal)
        z_llim, z_ulim, z_vals = \
            make_axis.quicklook_z(sig = z_vals, 
                                  y_vals = y_vals,
                                  z_lims = args['z_lims'] , 
                                  use_log = args['use_log_scale'],
                                  z_max_zone = args['z_max_zone'],
                                  z_min_zone = args['z_min_zone'])
                
        # Make title
        title = make_title.quicklook(start_date = data.RawData_Start_Date,
                                    start_time = data.RawData_Start_Time_UT, 
                                    end_time = data.RawData_Stop_Time_UT, 
                                    lidar = data.Lidar_Name, 
                                    channel = ch, 
                                    zan = data.Laser_Pointing_Angle,
                                    loc = data.Station_Name,
                                    smooth = args['smooth'],
                                    sm_lims = args['smoothing_range'],
                                    sm_win = args['smoothing_window'],
                                    sm_expo = args['smooth_exponential'])
      
  
        # Make filename
        fname = f'{data.Measurement_ID}_{data.Lidar_Name}_{meas_type}_qck_{ch}_ATLAS_{__version__}.png'

        # Make the plot
        fpath = make_plot.quicklook(dir_out = args['output_folder'], 
                                    fname = fname, title = title,
                                    dpi_val = args['dpi'],
                                    color_reduction = args['color_reduction'],
                                    use_log = args['use_log_scale'],
                                    x_vals = x_vals, t_vals = t_vals, 
                                    y_vals = y_vals, z_vals = z_vals, 
                                    x_lbin = x_lbin, x_ubin = x_ubin, 
                                    y_lbin = y_lbin, y_ubin = y_ubin, 
                                    y_llim = y_llim, y_ulim = y_ulim, 
                                    z_llim = z_llim, z_ulim = z_ulim,
                                    x_label = x_label, t_label = t_label, 
                                    y_label = y_label, 
                                    x_tick = x_tick, t_tick = t_tick, 
                                    y_tick = args['y_tick'], nodes = nodes)  
    
     
        # Add metadata to the quicklook plot
        from PIL import Image
        from PIL import PngImagePlugin
       
        METADATA = {"processing_software" : f"ATLAS_{data.version}",
                    "measurement_id" : f"{data.Measurement_ID}",
                    "channel" : f"{ch}",
                    "smooth" : f"{args['smooth']}",
                    "smoothing_exponential" : f"{args['smooth_exponential']}",
                    "smoothing_range" : f"{args['smoothing_range']}",
                    "smoothing_window": f"{args['smoothing_window']}",
                    "dpi" : f"{args['dpi']}",
                    "color_reduction" : f"{args['color_reduction']}",
                    "use_log_scale" : f"{args['use_log_scale']}",
                    "use_distance" : f"{args['use_distance']}",
                    "z_max_zone" : f"{args['z_max_zone']}",
                    "x_lims" : f"{args['x_lims']}",
                    "y_lims" : f"{args['y_lims']}",
                    "z_lims" : f"{args['z_lims']}",
                    "x_tick" : f"{args['x_tick']}",
                    "y_tick" : f"{args['y_tick']}"}
                
        im = Image.open(fpath)
        meta = PngImagePlugin.PngInfo()
    
        for x in METADATA.keys():
            meta.add_text(x, METADATA[x])
            
        im.save(fpath, "png", pnginfo = meta)

    print('-----------------------------------------')
    print(' ')

    return()

if __name__ == '__main__':
    
    sys.path.append('../')
    
    from version import __version__
    
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args, __version__)
    