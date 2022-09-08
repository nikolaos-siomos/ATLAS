#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, os, sys
import xarray as xr
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from readers.parse_config import polarization_calibration_parser, check_channels
from plotting import make_axis, make_title, make_plot
from lidar_processing.smoothing import sliding_average_1D
from lidar_processing import average

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the command line argument information
args = polarization_calibration_parser()

# Read the quicklook file
data = xr.open_dataset(args['input_file'])

# Extract signal
sig_ray = data.Range_Corrected_Signals_Rayleigh
sig_ray = sig_ray.copy().where(sig_ray != nc.default_fillvals['f8'])

sig_m45 = data.Range_Corrected_Signals_minus_45
sig_m45 = sig_m45.copy().where(sig_m45 != nc.default_fillvals['f8'])\
    .mean(dim = 'time_m45')

sig_p45 = data.Range_Corrected_Signals_plus_45
sig_p45 = sig_p45.copy().where(sig_p45 != nc.default_fillvals['f8'])\
    .mean(dim = 'time_p45')

# Extract signal time, channels, and bins
channels = data.channel.values
channels_r = args['ch_r']
channels_t = args['ch_t']

# Extract the K value
if args['K'] == None:
    K = len(channels_r) * [1.]
else:
    K = args['K']

# Check if the parsed channels exist
channels_r = check_channels(sel_channels = channels_r, all_channels = channels)
channels_t = check_channels(sel_channels = channels_t, all_channels = channels)

# Extract Molecular Depolarization Ratio and Calucalte the Atm. Parameter alpha
mldr = data.Molecular_Linear_Depolarization_Ratio
a_m = (1. - mldr) / (1. + mldr)

# Iterate over the channels
for ch_r, ch_t, K_ch in zip(channels_r, channels_t, K):
    
    ch_r_d = dict(channel = ch_r)
    ch_t_d = dict(channel = ch_t)
    
    sig_r_p45_ch = sig_p45.loc[ch_r_d].values
    sig_t_p45_ch = sig_p45.loc[ch_t_d].values
    sig_r_m45_ch = sig_m45.loc[ch_r_d].values
    sig_t_m45_ch = sig_m45.loc[ch_t_d].values
    sig_r_ray_ch = sig_ray.loc[ch_r_d].values
    sig_t_ray_ch = sig_ray.loc[ch_t_d].values
    
    a_m_ch = a_m.loc[ch_r_d].values

    # Create the y axis (height/range)
    x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
        make_axis.rayleigh_x(heights = data.Height_levels_Calibration.loc[ch_r_d].values, 
                             ranges = data.Range_levels_Calibration.loc[ch_r_d].values,
                             x_lims = args['x_lims'], 
                             use_dis = args['use_distance'])

    x_lbin_ray, x_ubin_ray, x_llim_ray, x_ulim_ray, x_vals_ray, x_label_ray = \
        make_axis.rayleigh_x(heights = data.Height_levels_Rayleigh.loc[ch_r_d].values, 
                             ranges = data.Range_levels_Rayleigh.loc[ch_r_d].values,
                             x_lims = args['x_lims'], 
                             use_dis = args['use_distance'])

    # Smoothing
    if args['smooth']:
        # Smoothing averaged sectors
        y_r_m45, _ = \
            sliding_average_1D(y_vals = sig_r_m45_ch, 
                               x_vals = x_vals,
                               x_sm_lims = args['smoothing_range'],
                               x_sm_hwin = args['half_window'],
                               expo = args['smooth_exponential'])

        y_t_m45, _ = \
            sliding_average_1D(y_vals = sig_t_m45_ch, 
                               x_vals = x_vals,
                               x_sm_lims = args['smoothing_range'],
                               x_sm_hwin = args['half_window'],
                               expo = args['smooth_exponential'])
            
        y_r_p45, _ = \
            sliding_average_1D(y_vals = sig_r_p45_ch, 
                               x_vals = x_vals,
                               x_sm_lims = args['smoothing_range'],
                               x_sm_hwin = args['half_window'],
                               expo = args['smooth_exponential'])

        y_t_p45, _ = \
            sliding_average_1D(y_vals = sig_t_p45_ch, 
                               x_vals = x_vals,
                               x_sm_lims = args['smoothing_range'],
                               x_sm_hwin = args['half_window'],
                               expo = args['smooth_exponential'])

        y_r_ray, _ = \
            sliding_average_1D(y_vals = sig_r_ray_ch, 
                               x_vals = x_vals,
                               x_sm_lims = args['smoothing_range'],
                               x_sm_hwin = args['half_window'],
                               expo = args['smooth_exponential'])    
            
        y_t_ray, _ = \
            sliding_average_1D(y_vals = sig_t_ray_ch, 
                               x_vals = x_vals,
                               x_sm_lims = args['smoothing_range'],
                               x_sm_hwin = args['half_window'],
                               expo = args['smooth_exponential'])            

    eta_ch = np.sqrt((sig_r_p45_ch * sig_r_m45_ch) / \
                     (sig_t_p45_ch * sig_t_m45_ch)) / K_ch

    eta_m45_ch = (sig_r_m45_ch / sig_t_m45_ch) / K_ch

    eta_p45_ch = (sig_r_p45_ch / sig_t_p45_ch) / K_ch
    
    sys.exit(0)
    
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
