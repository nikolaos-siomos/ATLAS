#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, os, sys, glob
import xarray as xr
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from readers.parse_config import rayleigh_parser, check_channels
from plotting import make_axis, make_title, make_plot
from lidar_processing.smoothing import sliding_average_1D
from lidar_processing import normalize

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the command line argument information
args = rayleigh_parser()

# Read the quicklook file
data = xr.open_dataset(args['input_file'])

# Delete all existing png files within
pngs = glob.glob(os.path.join(args['output_folder'],'*.png'))
for file in pngs:
    os.remove(file)

# Extract signal
sig = data.Range_Corrected_Signals
sig = sig.copy().where(sig != nc.default_fillvals['f8'])

# Extract Attenuated Backscatter
atb = data.Attenuated_Backscatter

# Extract signal time, channels, and bins
time = sig.time.values
bins = sig.bins.values

# Check if the parsed channels exist
channels = check_channels(sel_channels = args['channels'], 
                          all_channels = data.channel.values)

# iterate over the channels
for ch in channels:
    
    ch_d = dict(channel = ch)
    sig_ch = sig.loc[ch_d].values
    atb_ch = atb.loc[ch_d].values

    # Create the y axis (height/range)
    x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
        make_axis.rayleigh_x(heights = data.Height_levels.loc[ch_d].values, 
                             ranges = data.Range_levels.loc[ch_d].values,
                             x_lims = args['x_lims'], 
                             use_dis = args['use_distance'])

    # Smoothing
    if args['smooth']:
        y_vals_sm, y_vals_err = \
            sliding_average_1D(y_vals = sig_ch.copy(), 
                               x_vals = x_vals,
                               x_sm_lims = args['smoothing_range'],
                               x_sm_hwin = args['half_window'],
                               expo = args['smooth_exponential'])
    
        # Normalization for smoothed signals
        y_vals_sm = normalize.to_a_point(sig = y_vals_sm, 
                                         sig_b = atb_ch.copy(), 
                                         x_vals = x_vals,
                                         norm = args['reference_height'],
                                         hwin = args['half_reference_window'],
                                         axis = 0)

    
    # Normalization for unsmoothed signals
    y_vals = normalize.to_a_point(sig = sig_ch.copy(), 
                                  sig_b = atb_ch.copy(), 
                                  x_vals = x_vals,
                                  norm = args['reference_height'],
                                  hwin = args['half_reference_window'],
                                  axis = 0)

    # Create the y axis (signal)
    y_llim, y_ulim, y_label = \
        make_axis.rayleigh_y(sig = y_vals[slice(x_lbin,x_ubin+1)], 
                             atb = atb_ch.copy(), 
                             y_lims = args['y_lims'] , 
                             use_lin = args['use_lin_scale'])
    
            
    # Make title
    title = make_title.rayleigh(start_date = data.RawData_Start_Date,
                                start_time = data.RawData_Start_Time_UT, 
                                end_time = data.RawData_Stop_Time_UT, 
                                lidar = data.Lidar_Name, 
                                channel = ch, 
                                zan = data.Laser_Pointing_Angle,
                                lat = data.Latitude_degrees_north, 
                                lon = data.Longitude_degrees_east, 
                                elv = data.Altitude_meter_asl)

    # Make filename
    fname = f'ray_{data.Measurement_ID}_{ch}.png'

    # Make the plot
    fpath = make_plot.rayleigh(dir_out = args['output_folder'], 
                               fname = fname, title = title,
                               dpi_val = args['dpi'],
                               use_lin = args['use_lin_scale'],
                               x_refr = args['reference_height'],
                               refr_hwin = args['half_reference_window'],
                               x_vals = x_vals, y1_vals = y_vals_sm,
                               y2_vals = atb_ch.copy(), y3_vals = y_vals,
                               x_lbin = x_lbin, x_ubin = x_ubin,
                               x_llim = x_llim, x_ulim = x_ulim, 
                               y_llim = y_llim, y_ulim = y_ulim, 
                               x_label = x_label, y_label = y_label,
                               x_tick = args['x_tick'])  

    # sys.exit()
    # # Add metadata to the quicklook plot
    # from PIL import Image
    # from PIL import PngImagePlugin
   
    # METADATA = {"processing_software" : f"ATLAS_{data.version}",
    #             "measurement_id" : f"{data.Measurement_ID}",
    #             "channel" : f"{ch}",
    #             "smooth" : f"{args['smooth']}",
    #             "smoothing_exponential" : f"{args['smooth_exponential']}",
    #             "smoothing_range (lower)" : f"{args['smoothing_range'][0]}",
    #             "smoothing_range (upper)" : f"{args['smoothing_range'][-1]}",
    #             "half_window (lower)": f"{args['half_window'][0]}",
    #             "half_window (upper)": f"{args['half_window'][-1]}",
    #             "dpi" : f"{args['dpi']}",
    #             "use_log_scale" : f"{args['use_log_scale']}",
    #             "use_distance" : f"{args['use_distance']}",
    #             "x_lims (lower)" : f"{x_llim}",
    #             "x_lims (upper)" : f"{x_ulim}",
    #             "y_lims (lower)" : f"{y_vals[y_llim]}",
    #             "y_lims (upper)" : f"{y_vals[y_ulim]}",
    #             "z_lims (lower)" : f"{z_llim}",
    #             "z_lims (upper)" : f"{z_ulim}",
    #             "x_tick" : f"{x_tick}",
    #             "y_tick" : f"{args['y_tick']}"}
            
    # im = Image.open(fpath)
    # meta = PngImagePlugin.PngInfo()

    # for x in METADATA.keys():
    #     meta.add_text(x, METADATA[x])
        
    # im.save(fpath, "png", pnginfo = meta)
