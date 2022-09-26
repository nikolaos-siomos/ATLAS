#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 13:01:51 2022

@author: nick
"""

import numpy as np

def quicklook_x(x_lims, x_tick, t_tick, time):
    
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
    x_label = ' '
    t_label = 'Time UTC'

    # Calculate the x_tick (number of timeframes) if not provided
    if x_tick == None:
        if time.size / 15. > 10.:
            x_tick = np.round(time.size / 15., decimals = -1)
        else:
            x_tick = np.round(time.size / 15., decimals = 0)
        if x_tick == 0:
            x_tick = 1.
        
    # Calculate the t_tick (in minutes) if not provided
    if t_tick == None:
        mins = \
            (time[-1]-time[0]).astype('timedelta64[m]') / np.timedelta64(5,'m') 
        t_tick = 5. * np.round(mins / 15., decimals = 0)
    
        if t_tick == 0:
            t_tick = 1
    
    # Identify bins where temporal gaps are encountered (10% acceptance)
    nodes = np.where(time[1:]-time[:-1] > 1.15 * np.min(time[1:]-time[:-1]))[0]

    
    return(x_llim, x_ulim, x_vals, t_vals, x_label, t_label, 
           x_tick, t_tick, nodes)

def quicklook_y(heights, ranges, y_lims, use_dis):

    # Use altitude or distance for the y axis  
    if use_dis:
        y_vals = 1E-3 * ranges
        
        y_label = 'Distance from the Lidar [km]'
        
    else:
        y_vals = 1E-3 * heights       
        y_label = 'Altitude [km]'

    # Get the altitude/distance lower limit and bin
    if y_lims[0] == None or y_lims[0] < y_vals[0]:
        y_lbin = 0
        y_llim = np.round(y_vals[y_lbin], decimals = 2)

    else:
        y_lbin = np.where(y_vals >= y_lims[0])[0][0]
        
        if y_lbin > 0:
            y_lbin = y_lbin - 1
        
        y_llim = y_lims[0]


    # Get the altitude/distance upper limit and bin
    if y_lims[-1] == None or y_lims[-1] > y_vals[-1]:
        y_ubin = y_vals.size - 1
        y_ulim = np.round(y_vals[y_ubin], decimals = 2)

    else:
        y_ubin = np.where(y_vals <= y_lims[-1])[0][-1] 
        
        if y_ubin < y_vals.size:
            y_ubin = y_ubin + 1

        y_ulim = y_lims[-1]

    return(y_lbin, y_ubin, y_llim, y_ulim, y_vals, y_label)

def quicklook_z(sig, y_vals, z_lims, use_log, z_max_zone):
    
    # Get the max signal bin and value
    sig_m = sig.rolling(bins = 100, center =True).mean().values
    z_vals = sig.transpose('bins','time').values
       
    if z_max_zone[0] == None: 
        z_max = np.nanmax(sig_m)
    else:
        mask = (y_vals >= z_max_zone[0]) & (y_vals <= z_max_zone[1])
        z_max = np.nanmax(sig_m[:,mask])
    # from matplotlib import pyplot as plt
    # plt.plot(sig_m[mask])
    # plt.show()

    z_vals = z_vals / z_max

    # Get the signal upper limit
    if z_lims[-1] == None:
        z_ulim = 1.
        
    else:
       z_ulim = z_lims[-1]
    
    # Get the vertical lower limit
    if use_log and z_lims[0] == None:
        z_llim = z_ulim * 1E-5
        
    elif use_log == False and z_lims[0] == None:
        z_llim = 0.
    
    else:
       z_llim = z_lims[0]
    
    # Get the y axis labels
    z_label = 'Range-corrected Signal [A.U.]'
    
    return(z_llim, z_ulim, z_label, z_vals)


def rayleigh_x(heights, ranges, x_lims, use_dis):

    # Use altitude or distance for the y axis  
    if use_dis:
        x_vals = 1E-3 * ranges
        
        x_label = 'Distance from the Lidar [km]'
        
    else:
        x_vals = 1E-3 * heights       
        x_label = 'Altitude [km]'

    # Get the altitude/distance lower limit and bin
    if x_lims[0] == None or x_lims[0] < x_vals[0]:
        x_lbin = 0
        x_llim = np.round(x_vals[x_lbin], decimals = 2)

    else:
        x_lbin = np.where(x_vals >= x_lims[0])[0][0]
        
        if x_lbin > 0:
            x_lbin = x_lbin - 1
        
        x_llim = x_lims[0]

    # Get the altitude/distance upper limit and bin
    if x_lims[-1] == None or x_lims[-1] > x_vals[-1]:
        x_ubin = x_vals.size - 1
        x_ulim = np.round(x_vals[x_ubin], decimals = 2)

    else:
        x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
        
        if x_ubin < x_vals.size:
            x_ubin = x_ubin + 1

        x_ulim = x_lims[-1]

    return(x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label)

def rayleigh_y(sig, atb, y_lims, use_lin):
    
    # Get the max signal bin and value       
    y_max = np.nanmax([np.nanmax(sig),np.nanmax(atb)])
    y_min = np.nanmin(np.abs(atb))

    # Get the signal upper limit
    if use_lin == False:
        if y_lims[-1] == None:
            if np.isnan(y_max) or np.isinf(y_max) or y_max <= 0:
                y_ulim = 1
            else:
                y_ulim = 2. * y_max
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_ulim = 1
            else:
                y_ulim =  y_lims[-1]
        
    # Get the vertical lower limit
    if use_lin == False:
        if y_lims[0] == None:
            if np.isnan(y_min) or np.isinf(y_min) or y_min <= 0:
                y_llim = 1E-3
            else:
                y_llim = 0.5 * y_min
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_llim = 0.
            else:
                y_llim =  y_lims[0]
    
    # Get the y axis labels
    y_label = 'Attenuated Backscatter [$m^{-1} sr^{-1}$]'
    
    return(y_llim, y_ulim, y_label)

def telecover_x(heights, ranges, x_lims, use_dis):

    # Use altitude or distance for the y axis  
    if use_dis:
        x_vals = 1E-3 * ranges
        
        x_label = 'Distance from the Lidar [km]'
        
    else:
        x_vals = 1E-3 * heights       
        x_label = 'Altitude [km]'

    # Get the altitude/distance lower limit and bin
    if x_lims[0] == None or x_lims[0] < x_vals[0]:
        x_lbin = 0
        x_llim = x_vals[x_lbin]
        x_llim = np.round(x_vals[x_lbin], decimals = 2)

    else:
        x_lbin = np.where(x_vals >= x_lims[0])[0][0]
        
        if x_lbin > 0:
            x_lbin = x_lbin - 1
        
        x_llim = x_lims[0]

    # Get the altitude/distance upper limit and bin
    if x_lims[-1] == None or x_lims[-1] > x_vals[-1]:
        x_ubin = x_vals.size - 1
        x_ulim = np.round(x_vals[x_ubin], decimals = 2)


    else:
        x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
        
        if x_ubin < x_vals.size:
            x_ubin = x_ubin + 1

        x_ulim = x_lims[-1]
        

    return(x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label)

def telecover_y(sig, sig_nr, y_lims):
    
    # Get the max signal bin and value       
    y_max = np.nanmax(sig)

    y_max_nr = np.nanmax(sig_nr)
    
    coef = y_max / y_max_nr 
    
    # Get the signal upper limit
    if y_lims[-1] == None:
        y_ulim = 1.2 * y_max_nr * coef
        y_ulim_nr = 1.2 * y_max_nr
        
    else:
       y_ulim = y_lims[-1] * coef
       y_ulim_nr = y_lims[-1]
      
    if np.isnan(y_ulim) or y_ulim <= 0.:
        y_ulim = 1.

    if np.isnan(y_ulim_nr) or y_ulim_nr <= 0.:
        y_ulim_nr = 1.
    
    # Get the vertical lower limit
    if y_lims[0] == None:
        y_llim = 0.
        y_llim_nr = 0.
        
    else:
       y_llim = y_lims[0] * coef
       y_llim_nr = y_lims[0] 
    
    return(y_llim, y_ulim, y_llim_nr, y_ulim_nr)


def polarization_calibration_cal_x(ratio_m, ratio_p, x_lims_cal):
    
    # Get the max signal bin and value       
    x_max_cal = np.nanmax([ratio_m, ratio_p])
    x_min_cal = np.nanmin([ratio_m, ratio_p])

    # Get the eta upper limit
    if x_lims_cal[-1] == None:
        x_ulim_cal = 1.2 * x_max_cal
        
    else:
       x_ulim_cal = x_lims_cal[-1]
    
    # Get the vertical lower limit
    if x_lims_cal[0] == None:
        x_llim_cal = 0.8 * x_min_cal
        
    else:
       x_llim_cal = x_lims_cal[0]
    
    return(x_llim_cal, x_ulim_cal)

def polarization_calibration_ray_x(ratio, x_lims_ray):
    
    # Get the max signal bin and value       
    x_max_ray = np.nanmax(ratio)

    # Get the delta upper limit
    if x_lims_ray[-1] == None:
        x_ulim_ray = 1.2 * x_max_ray
        
    else:
       x_ulim_ray = x_lims_ray[-1]
    
    # Get the vertical lower limit
    if x_lims_ray[0] == None:
        x_llim_ray = 0.
        
    else:
       x_llim_ray = x_lims_ray[0]
    
    return(x_llim_ray, x_ulim_ray)


def polarization_calibration_y(heights, ranges, y_lims, use_dis):

    # Use altitude or distance for the y axis  
    if use_dis:
        y_vals = 1E-3 * ranges
        
        y_label = 'Distance from the Lidar [km]'
        
    else:
        y_vals = 1E-3 * heights       
        y_label = 'Altitude [km]'

    # Get the altitude/distance lower limit and bin
    if y_lims[0] == None or y_lims[0] < y_vals[0]:
        y_lbin = 0
        y_llim = np.round(y_vals[y_lbin], decimals = 2)


    else:
        y_lbin = np.where(y_vals >= y_lims[0])[0][0]
        
        if y_lbin > 0:
            y_lbin = y_lbin - 1
        
        y_llim = y_lims[0]


    # Get the altitude/distance upper limit and bin
    if y_lims[-1] == None or y_lims[-1] > y_vals[-1]:
        y_ubin = y_vals.size - 1
        y_ulim = np.round(y_vals[y_ubin], decimals = 2)

    else:
        y_ubin = np.where(y_vals <= y_lims[-1])[0][-1] 
        
        if y_ubin < y_vals.size:
            y_ubin = y_ubin + 1

        y_ulim = y_lims[-1]
        

    return(y_lbin, y_ubin, y_llim, y_ulim, y_vals, y_label)


def intercomparison_x(heights, ranges, x_lims, use_dis):

    # Use altitude or distance for the y axis  
    if use_dis:
        x_vals = 1E-3 * ranges
        
        x_label = 'Distance from the Lidar [km]'
        
    else:
        x_vals = 1E-3 * heights       
        x_label = 'Altitude [km]'

    # Get the altitude/distance lower limit and bin
    if x_lims[0] == None or x_lims[0] < x_vals[0]:
        x_lbin = 0
        x_llim = np.round(x_vals[x_lbin], decimals = 2)

    else:
        x_lbin = np.where(x_vals >= x_lims[0])[0][0]
        
        if x_lbin > 0:
            x_lbin = x_lbin - 1
        
        x_llim = x_lims[0]


    # Get the altitude/distance upper limit and bin
    if x_lims[-1] == None or x_lims[-1] > x_vals[-1]:
        x_ubin = x_vals.size - 1
        x_ulim = np.round(x_vals[x_ubin], decimals = 2)

    else:
        x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
        
        if x_ubin < x_vals.size:
            x_ubin = x_ubin + 1

        x_ulim = x_lims[-1]
        

    return(x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label)

def intercomparison_y(sig1, sig2, y_lims, use_lin):
    
    # Get the max signal bin and value       
    y_max = np.nanmax([np.nanmax(sig1[:int(sig1.size/2)]),
                       np.nanmax(sig2[:int(sig2.size/2)])])
    y_min = np.nanmin([np.nanmin(sig1[int(sig1.size/2):]),
                       np.nanmin(sig2[int(sig2.size/2):])])

    # Get the signal upper limit
    if use_lin == False:
        if y_lims[-1] == None:
            if np.isnan(y_max) or np.isinf(y_max) or y_max <= 0:
                y_ulim = 1
            else:
                y_ulim = 2. * y_max
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_ulim = 1
            else:
                y_ulim =  y_lims[-1]
        
    # Get the vertical lower limit
    if use_lin == False:
        if y_lims[0] == None:
            if np.isnan(y_min) or np.isinf(y_min) or y_min <= 0:
                y_llim = 1E-3
            else:
                y_llim = 0.5 * y_min
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_llim = 0.
            else:
                y_llim =  y_lims[0]
    
    # Get the y axis labels
    y_label = 'Attenuated Backscatter [$m^{-1} sr^{-1}$]'
    
    return(y_llim, y_ulim, y_label)