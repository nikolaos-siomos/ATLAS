#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import warnings, os, sys
import xarray as xr
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from readers.parse_config import quicklook_parser
from plotting import make_axis, color_lib, make_colormap, make_title
import matplotlib.dates as mdates

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the command line argument information
args = quicklook_parser()

# Read the quicklook file
data = xr.open_dataset(args['input_file'])

# Extract variables
sig = data.Range_Corrected_Signals
sig = sig.copy().where(sig != nc.default_fillvals['f8'])
zan = np.round(data.Laser_Pointing_Angle, decimals=1)
lat = np.round(data.Latitude_degrees_north, decimals=1)
lon = np.round(data.Longitude_degrees_east)
elv = np.round(data.Altitude_meter_asl, decimals=0)
meas_id = data.Measurement_ID
lidar = data.Lidar_Name

channels = data.channel.values
sel_channels = args['channels'] 
if sel_channels != None:
    if any([ch not in channels for ch in sel_channels]):
        sys.exit("-- Error: The following provided channels do not exist: "+\
                 "{sel_channels[[ch not in channels for ch in sel_channels]]}")

time = sig.time.values
bins = sig.bins.values

dpi_val = args['dpi']    
use_log = args['use_log_scale'] 
use_dis = args['use_distance']
z_lvls = args['z_lvls'] 
x_lims = args['x_lims'] 
x_tick = args['x_tick'] 
y_lims = args['y_lims'] 
y_tick = args['y_tick']

# Get x axis info
x_llim, x_ulim, x_vals, t_vals, x_label, t_label = \
    make_axis.quicklook_x(x_lims, time)

if x_tick == None:
    x_tick = np.round(time.size / 15., decimals = -1)
    mins = (time[-1]-time[0]).astype('timedelta64[m]') / np.timedelta64(5,'m')
    t_tick = 5. * np.round(mins / 15., decimals = 0)

# iterate over the channels
for ch in channels:
    
    ch_d = dict(channel = ch)
    sig_ch = sig.copy().loc[ch_d]
    
    # Use altitude or distance for the y axis  
    if use_dis:
        yvals = data.Distance_levels.loc[ch_d].values
        
    else:
        yvals = data.Height_levels.loc[ch_d].values

    y_llim, y_ulim, y_vals, y_label = \
        make_axis.quicklook_y(y_vals = yvals, 
                              y_lims = y_lims, 
                              use_dis = use_dis)
    
    sel = dict(time = slice(x_llim, x_ulim), bins = slice(y_llim, y_ulim))
            
    z_llim, z_ulim, z_label, z_vals = \
        make_axis.quicklook_z(sig = sig_ch, use_log = use_log)
        
    rgb = color_lib.volkers_rgb()
    my_cmap = make_colormap.custom_rgb(rgb, name = 'volkers')

    # Make the quicklook figures
    [X, Y] = np.meshgrid(t_vals[slice(x_llim, x_ulim)], 
                         y_vals[slice(y_llim, y_ulim)])
    Z = z_vals[slice(y_llim, y_ulim), slice(x_llim, x_ulim)]
    
    title = make_title.quicklook(start_time = t_vals[x_llim], 
                                 end_time = t_vals[x_ulim-1], 
                                 lidar = lidar, zan = zan, 
                                 lat = lat, lon = lon, elv = elv)
    #+ np.round(time.size/100, decimals = 1)
    fig = plt.figure(figsize=(9. , 4.5))
    ax = fig.subplots()
    plt.title(title, pad = 10)
    qck = ax.pcolormesh(X, Y, Z, vmin = z_llim, vmax = z_ulim, cmap = my_cmap)
    cbar = fig.colorbar(qck, label = z_label, extend = 'both')
    ax.set_xlabel(t_label)
    ax.set_ylabel(y_label)
    ax.set_yticks(np.arange(0, y_vals[y_ulim] + y_tick, y_tick))
    ax.set_ylim([y_vals[y_llim],y_vals[y_ulim]])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval = int(t_tick)))
    plt.xticks(rotation = 35)
    ax1 = ax.twiny()
    ax1.set_xlabel(x_label)
    ax1.set_xticks(np.arange(x_llim, x_ulim + x_tick, x_tick))

    # if y_tick != None:
    #     y_ticks = np.arange(y_llim, y_ulim, y_tick)
    #     plt.xticks(y_ticks,
    #                labels = np.round(y_ticks, decimals = 3).astype(str)
    fname = f'qck_{meas_id}_{ch}.png'
    plt.tight_layout()
    plt.savefig(os.path.join(args['output_folder'],fname),
                dpi = args['dpi'])
    plt.close()
    