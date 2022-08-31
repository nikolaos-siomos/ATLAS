#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 13:01:51 2022

@author: nick
"""

import numpy as np

def quicklook_x(x_lims, time):
    
    # Get the x lower limit
    if x_lims[0] == None or x_lims[0] < 0:
        x_llim = 0
    else:
        x_llim = x_lims[0] - 1

    # Get the x upper limit
    if x_lims[-1] == None or x_lims[-1] > time.size:
        x_ulim = time.size
    else:
        x_ulim = x_lims[-1]
    
    # Get the timeframe levels
    x_vals = np.arange(0, time.size, 1)
    t_vals = time

    # Get the x axis labels
    x_label = 'Number of Timeframes'
    t_label = 'Time UTC'
    
    return(x_llim, x_ulim, x_vals, t_vals, x_label, t_label)

def quicklook_y(y_vals, y_lims, use_dis):
    
    y_vals = 1E-3 * y_vals

    # Get the vertical lower limit
    if y_lims[0] == None or y_lims[0] < y_vals[0]:
        y_llim = 0
    else:
        y_llim = np.where(y_vals >= y_lims[0])[0][0]

    # Get the vertical upper limit
    if y_lims[-1] == None or y_lims[-1] > y_vals[-1]:
        y_ulim = 0
    else:
        y_ulim = np.where(y_vals <= y_lims[-1])[0][-1] 

    # Get the y axis labels
    if use_dis:
        y_label = 'Distance from the Lidar [km]'
    else:
        y_label = 'Altitude [km]'
    
    return(y_llim, y_ulim, y_vals, y_label)

def quicklook_z(sig, use_log):
    
    # Get the max signal bin and value
    sig_m = sig.mean(dim = 'time').rolling(bins = 100, center =True).mean()
       
    z_maxz = sig_m.max().values
    z_vals = sig.transpose('bins','time').values / z_maxz

    # Get the signal upper limit
    z_ulim = 1.
    
    # Get the vertical lower limit
    if use_log:
        z_llim = z_ulim * 1E-5
    else:
        z_llim = 0.
    
    # Get the y axis labels
    z_label = 'Range-corrected Signal [A.U.]'
    
    return(z_llim, z_ulim, z_label, z_vals)
    