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
from PIL import Image
from PIL import PngImagePlugin
       
# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__, meas_type):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print(f'Start generating Quicklooks ({meas_type})...')
    print('-----------------------------------------')

    # Read the quicklook file
    data = xr.open_dataset(args['input_file'])
    
    # Extract signal
    sig = data.Range_Corrected_Signals
    sig = sig.copy().where(sig != nc.default_fillvals['f8'])

    # Extract signal time, channels, and bins
    time = sig.time.values
    
    # Extract SCC info
    station_id = data.Station_ID.lower()
    lidar_id = data.Lidar_ID
    version_id = data.Version_ID
    config_id = data.Configuration_ID
    config_name = data.Configuration_Name
    scc_id = data.channel_ID

    # Extract date info
    start_date = data.RawData_Start_Date
    start_time = data.RawData_Start_Time_UT
    stop_time = data.RawData_Stop_Time_UT

    delta_t = data.Raw_Data_Stop_Time - data.Raw_Data_Start_Time

    # Check if the parsed channels exist
    channels = \
        check_channels(sel_channels = args['channels'], 
                       all_channels = data.channel.values,
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
        sig_ch = sig.copy().loc[ch_d]
        scc_id_ch = scc_id.copy().loc[ch_d].values
            
        # Create the y axis (height/range)
        y_lbin, y_ubin, y_llim, y_ulim, y_vals, y_label = \
            make_axis.quicklook_y(heights = data.Height_levels.loc[ch_d].values, 
                                  ranges = data.Range_levels.loc[ch_d].values,
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
                                    station_id = station_id,
                                    lidar_id = lidar_id,
                                    version_id = version_id,
                                    config_id = config_id,
                                    config_name = config_name,
                                    scc_id = scc_id_ch,
                                    zan = data.Laser_Pointing_Angle,
                                    loc = data.Station_Name,
                                    smooth = args['smooth'],
                                    sm_lims = args['smoothing_range'],
                                    sm_win = args['smoothing_window'],
                                    sm_expo = args['smooth_exponential'])
      
  
        # Make filename
        fname = f'{station_id}_{lidar_id}_{version_id}_{config_id}_{start_date}_{start_time}_qck_{meas_type}_{ch}_{scc_id_ch}_ATLAS_{__version__}.png'

        # Make the plot
        fpath = make_plot.quicklook(dir_out = os.path.join(args['output_folder'],'plots'), 
                                    fname = fname, title = title,
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
        METADATA = {"processing_software" : f"ATLAS_{data.version}",
                    "station_id" : f"{station_id}",
                    "lidar_id" : f"{lidar_id}",
                    "version_id" : f"{version_id}",
                    "config_id" : f"{config_id}",
                    "channel" : f"{ch}",
                    "scc_id" : f"{scc_id_ch}",
                    "smooth" : f"{args['smooth']}",
                    "smoothing_exponential" : f"{args['smooth_exponential']}",
                    "smoothing_range" : f"{args['smoothing_range']}",
                    "smoothing_window": f"{args['smoothing_window']}",
                    "dpi" : f"{args['dpi']}",
                    "color_reduction" : f"{args['color_reduction']}",
                    "use_log_scale" : f"{args['use_log_scale']}",
                    "use_range" : f"{args['use_range']}",
                    "z_max_zone" : f"{args['z_max_zone']}",
                    "y_lims" : f"{args['y_lims']}",
                    "z_lims" : f"{args['z_lims']}",
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
    