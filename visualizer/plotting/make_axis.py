#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 13:01:51 2022

@author: nick
"""

import numpy as np

def quicklook_x(t_lims, t_tick, time):

    # # Identify bins where temporal gaps are encountered (10% acceptance)
    # nodes = np.where(time[1:]-time[:-1] > 1.50 * np.nanmin(time[1:]-time[:-1]))[0]
    # print(nodes)
    # if nodes > 0:
    #     print('-- Warning: The quicklook will contain gaps as the dataset is not continuous. The ascending number of timeframes will not be display as a secondary x_axis')
    
    # Get the lowest time based on the t-lims
    if t_lims[0] == None:
        ltime = time[0]
    else:
        lmins = int(str(t_lims[0])[:2]) * 60 + int(str(t_lims[0])[2:])
        ltime = time.astype('datetime64[D]')[0] + np.timedelta64(lmins, 'm')

    if t_lims[-1] == None:
        utime = time[-1]
    else:
        umins = int(str(t_lims[-1])[:2]) * 60 + int(str(t_lims[-1])[2:])
        utime = time.astype('datetime64[D]')[-1] + np.timedelta64(umins, 'm')

    t_vals = time[(time >= ltime) & (time <= utime)]

    # Get the x lower limit
    x_lbin = np.where((time >= ltime))[0][0]
    
    # Get the x upper limit
    x_ubin = np.where((time >= utime))[0][0]
    
    # Calculate the x_tick (number of timeframes) if not provided
    if time.size / 15. > 10.:
        x_tick = np.round(time.size / 15., decimals = -1)
    else:
        x_tick = np.round(time.size / 15., decimals = 0)
    if x_tick == 0:
        x_tick = 1.
            
    # Get the timeframe levels
    t_vals = time

    # Calculate the t_tick (in minutes) if not provided
    if t_tick == None:
        mins = \
            (t_vals[x_ubin]-t_vals[x_lbin]).astype('timedelta64[m]') / np.timedelta64(1,'m') 
        if mins < 5:
            t_tick = 1.
        elif mins >= 5 and mins < 20:
            t_tick = 2.
        elif mins >= 20 and mins < 40.:
            t_tick = 4.
        elif mins >= 40 and mins < 120.:
            t_tick = 10.
        elif mins >= 120 and mins < 240.:
            t_tick = 20.
        elif mins >= 240 and mins < 480.:
            t_tick = 30.
        else:
            t_tick = 60.
    
    return(x_lbin, x_ubin, x_tick, t_vals, t_tick)

def quicklook_y(heights, ranges, y_lims, use_dis):

    # Use Height or range above the lidar for the y axis  
    if use_dis:
        y_vals = 1E-3 * ranges
        
        y_label = 'Range above the lidar [km]'
        
    else:
        y_vals = 1E-3 * heights       
        y_label = 'Height above the lidar [km]'

    # # Get the altitude/distance lower limit and bin
    # if y_lims[0] == None or y_lims[0] < y_vals[0]:
    #     y_lbin = 0
    #     y_llim = np.round(y_vals[y_lbin], decimals = 2)

    # else:
    #     y_lbin = np.where(y_vals >= y_lims[0])[0][0]
        
    #     if y_lbin > 0:
    #         y_lbin = y_lbin - 1
        
    #     y_llim = y_lims[0]


    # # Get the altitude/distance upper limit and bin
    # if y_lims[-1] == None or y_lims[-1] > y_vals[-1]:
    #     y_ubin = y_vals.size - 1
    #     y_ulim = np.round(y_vals[y_ubin], decimals = 2)

    # else:
    #     y_ubin = np.where(y_vals <= y_lims[-1])[0][-1] 
        
    #     if y_ubin < y_vals.size:
    #         y_ubin = y_ubin + 1

    #     y_ulim = y_lims[-1]
        
    # Get the altitude/distance lower limit and bin
    y_lbin = np.where(y_vals >= y_lims[0])[0][0]
    
    if y_lbin > 0:
        y_lbin = y_lbin - 1
    
    y_llim = y_lims[0]


    # Get the altitude/distance upper limit and bin
    y_ubin = np.where(y_vals <= y_lims[-1])[0][-1] 
    
    if y_ubin < y_vals.size:
        y_ubin = y_ubin + 1

    y_ulim = y_lims[-1]

    return(y_lbin, y_ubin, y_llim, y_ulim, y_vals, y_label)

def quicklook_z(z_vals, y_vals, z_lims, use_log, z_max_zone, z_min_zone):

    # Get the max signal bin and value
    mask_max_zone = (y_vals >= z_max_zone[0]) & (y_vals <= z_max_zone[1])

    z_vals_sm = np.nanmean(z_vals[:,mask_max_zone],axis = 0)
        
    z_max = round_it(np.nanmax(z_vals_sm),1)
    
    # Normalize with the max
    z_vals = z_vals / z_max

    # Get the signal upper and lower limits
    z_ulim = z_lims[-1]
    z_llim = z_lims[0]

    # Get the vertical lower limit
    if use_log and z_llim == None:

        # Get the min signal bin and value
        mask_min_zone = (y_vals >= z_min_zone[0]) & (y_vals <= z_min_zone[1])

        z_vals_sm = np.nanmean(z_vals[:,mask_min_zone],axis = 0)
        
        z_llim = round_it(np.nanmin(z_vals_sm),2)
        
    elif use_log == False and z_llim == None:
        z_llim = 0.

    return(z_llim, z_ulim, z_vals)


def rayleigh_x(heights, ranges, x_lims, use_dis):

    # Use Height or range above the lidar for the y axis  
    if use_dis:
        x_vals = 1E-3 * ranges
        
        x_label = 'Range above the lidar [km]'
        
    else:
        x_vals = 1E-3 * heights       
        x_label = 'Height above the lidar [km]'

    # # Get the altitude/distance lower limit and bin
    # if x_lims[0] == None or x_lims[0] < x_vals[0]:
    #     x_lbin = 0
    #     x_llim = np.round(x_vals[x_lbin], decimals = 2)

    # else:
    #     x_lbin = np.where(x_vals >= x_lims[0])[0][0]
        
    #     if x_lbin > 0:
    #         x_lbin = x_lbin - 1
        
    #     x_llim = x_lims[0]

    # # Get the altitude/distance upper limit and bin
    # if x_lims[-1] == None or x_lims[-1] > x_vals[-1]:
    #     x_ubin = x_vals.size - 1
    #     x_ulim = np.round(x_vals[x_ubin], decimals = 2)

    # else:
    #     x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
        
    #     if x_ubin < x_vals.size:
    #         x_ubin = x_ubin + 1

    #     x_ulim = x_lims[-1]
    
    # Get the altitude/distance lower limit and bin
    x_lbin = np.where(x_vals >= x_lims[0])[0][0]
    
    if x_lbin > 0:
        x_lbin = x_lbin - 1
    
    x_llim = x_lims[0]

    # Get the altitude/distance upper limit and bin
    x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
    
    if x_ubin < x_vals.size:
        x_ubin = x_ubin + 1

    x_ulim = x_lims[-1]

    return(x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label)

def rayleigh_y(sig, atb, y_lims, wave, use_lin):
    
    # Get the max signal bin and value       
    y_max = np.nanmax(atb)
    y_min = np.nanmin(atb)

    scale_f = wave / 355.
    scat_ratio_f = 2.5
    # Get the signal axis upper limit
    if use_lin == False:
        if y_lims[-1] == None:
            y_ulim = scat_ratio_f * scale_f * y_max
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_ulim = 1
            else:
                y_ulim =  y_lims[-1]
    else:
        if y_lims[-1] == None:
            y_ulim = scat_ratio_f * scale_f * y_max
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_ulim = 1
            else:
                y_ulim =  y_lims[-1]
        
    # Get the signal axis lower limit
    if use_lin == False:
        if y_lims[0] == None:
            y_llim = y_min / 2.
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_llim = 0.
            else:
                y_llim =  y_lims[0]
    else:
        if y_lims[0] == None:
            y_llim = y_min / 2.
        else:
            if y_lims[0] <= 0:
                print('-- Warning: rayleigh y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_llim = 0.
            else:
                y_llim =  y_lims[0]
    
    # Get the y axis labels
    y_label = 'Attn. Bsc. rel. to fit range [$m^{-1} sr^{-1}$]'
    
    return(y_llim, y_ulim, y_label)

def telecover_x(heights, ranges, x_lims, x_tick, use_dis, telescope_type):

    # Use Height or range above the lidar for the y axis  
    if use_dis:
        x_vals = 1E-3 * ranges
        
        x_label = 'Range above the lidar [km]'
        
    else:
        x_vals = 1E-3 * heights       
        x_label = 'Height above the lidar [km]'

    
    # Set the lower x limit
    if x_lims[0] == None:
        x_llim = 0.
    else:
        x_llim = x_lims[0]
        
    # Set the upper x limit depending on the telescope_type
    if x_lims[-1] == None:
        if telescope_type in ['l', 'm', 'n', 'x']:
            x_ulim = 2.5
        else:
            x_ulim = 5.
    else:
        x_ulim = x_lims[-1]

    # Set the x_tick depending on the telescope_type
    if x_tick == None and telescope_type in ['l', 'm', 'n', 'x']:
        x_tick = 0.5
    elif x_tick == None and telescope_type in ['f', 'g',' h']:
        x_tick = 1.

    # Get the altitude/distance lower limit and bin
    x_lbin = np.where(x_vals >= x_llim)[0][0]
    
    if x_lbin > 0:
        x_lbin = x_lbin - 1
    
    # Get the altitude/distance upper limit and bin
    x_ubin = np.where(x_vals <= x_ulim)[0][-1] 
    
    if x_ubin < x_vals.size:
        x_ubin = x_ubin + 1
    
    return(x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_tick, x_label)

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
      
    if not np.isfinite(y_ulim) or y_ulim <= 0.:
        y_ulim = 1.

    if not np.isfinite(y_ulim_nr) or y_ulim_nr <= 0.:
        y_ulim_nr = 1.
   
    # Get the vertical lower limit
    if y_lims[0] == None:
        y_llim = -0.1 * y_max_nr * coef
        y_llim_nr = -0.1 * y_max_nr
        
    else:
       y_llim = y_lims[0] * coef
       y_llim_nr = y_lims[0] 

      
    if not np.isfinite(y_llim) or y_llim <= 0.:
        y_llim = 0.

    if not np.isfinite(y_llim_nr) or y_llim_nr <= 0.:
        y_llim_nr = 0.
    
    return(y_llim, y_ulim, y_llim_nr, y_ulim_nr)


def polarization_calibration_cal_y(ratio_m, ratio_p, y_lims_cal):
    
    # Get the max signal bin and value       
    y_max_cal = np.nanmax([ratio_m, ratio_p])
    y_min_cal = np.nanmin([ratio_m, ratio_p])

    # Get the eta upper limit
    if np.isfinite(y_max_cal) and y_lims_cal[-1] == None:
        y_ulim_cal = y_max_cal * 2.
            
    elif np.isfinite(y_max_cal) and y_lims_cal[-1] != None:
        y_ulim_cal = y_lims_cal[-1]
        
    else:
        y_ulim_cal = 1.
    
    # Get the vertical lower limit
    if np.isfinite(y_min_cal) and y_lims_cal[0] == None:
        y_llim_cal = y_min_cal / 1.5
            
    elif np.isfinite(y_min_cal) and y_lims_cal[0] != None:
           y_llim_cal = y_lims_cal[0]
    
    else:
        y_llim_cal = 0.
          
    y_label_cal = r'Gain ratio $η^{\star}_{f}$'
        
    
    return(y_llim_cal, y_ulim_cal, y_label_cal)

def polarization_calibration_ray_y(ratio, y_lims_ray):
    
    # Get the max signal bin and value       
    y_max_ray = ratio
    
    # Get the delta upper limit
    if np.isfinite(y_max_ray) and y_lims_ray[-1] == None:
        y_ulim_ray = y_max_ray * 2.5
    
    elif np.isfinite(y_max_ray) and y_lims_ray[-1] != None:
        y_ulim_ray = y_lims_ray[-1]
            
    else:
        y_ulim_ray = 0.1
    
    # Get the vertical lower limit
    if y_lims_ray[0] == None:
        y_llim_ray = 0.
        
    else:
       y_llim_ray = y_lims_ray[0]
      
    y_label_ray = 'Linear Dep. Ratio'

    
    return(y_llim_ray, y_ulim_ray, y_label_ray)


def polarization_calibration_x(heights, ranges, x_lims, use_dis):

    # Use Height or range above the lidar for the y axis  
    if use_dis:
        x_vals = 1E-3 * ranges
        
        x_label = 'Range above the lidar [km]'
        
    else:
        x_vals = 1E-3 * heights       
        x_label = 'Height above the lidar [km]'

    # Get the altitude/distance lower limit and bin
    x_lbin = np.where(x_vals >= x_lims[0])[0][0]
    
    if x_lbin > 0:
        x_lbin = x_lbin - 1
    
    x_llim = x_lims[0]

    # Get the altitude/distance upper limit and bin
    x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
    
    if x_ubin < x_vals.size:
        x_ubin = x_ubin + 1

    x_ulim = x_lims[-1]
        

    return(x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label)


def intercomparison_x(heights, ranges, x_lims, use_dis):

    # Use Height or range above the lidar for the y axis  
    if use_dis:
        x_vals = 1E-3 * ranges
        
        x_label = 'Range above the lidar [km]'
        
    else:
        x_vals = 1E-3 * heights       
        x_label = 'Height above the lidar [km]'

    # Get the altitude/distance lower limit and bin
    x_lbin = np.where(x_vals >= x_lims[0])[0][0]
    
    if x_lbin > 0:
        x_lbin = x_lbin - 1
    
    x_llim = x_lims[0]

    # Get the altitude/distance upper limit and bin
    x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
    
    if x_ubin < x_vals.size:
        x_ubin = x_ubin + 1

    x_ulim = x_lims[-1]
                

    return(x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label)

# def intercomparison_y(sig1, sig2, y_lims, use_lin):
    
#     # Get the max signal bin and value       
#     y_max = np.nanmax([np.nanmax(sig1[:int(sig1.size/2)]),
#                        np.nanmax(sig2[:int(sig2.size/2)])])
#     y_min = np.nanmin([np.nanmin(sig1[int(sig1.size/2):]),
#                        np.nanmin(sig2[int(sig2.size/2):])])

#     # Get the signal upper limit
#     if use_lin == False:
#         if y_lims[-1] == None:
#             if not np.isfinite(y_max) or y_max <= 0:
#                 y_ulim = 1
#             else:
#                 y_ulim = 2. * y_max
#         else:
#             if y_lims[0] <= 0:
#                 print('-- Warning: rayleigh y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
#                 y_ulim = 1
#             else:
#                 y_ulim =  y_lims[-1]
        
#     # Get the vertical lower limit
#     if use_lin == False:
#         if y_lims[0] == None:
#             if not np.isfinite(y_min) or y_min <= 0:
#                 y_llim = y_ulim * 1E-3
#             else:
#                 y_llim = 0.5 * y_min
#         else:
#             if y_lims[0] <= 0:
#                 print('-- Warning: rayleigh y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
#                 y_llim = 0.
#             else:
#                 y_llim =  y_lims[0]
    
#     # Get the y axis labels
#     y_label = 'Attenuated Backscatter [$m^{-1} sr^{-1}$]'

#     return(y_llim, y_ulim, y_label)

def intercomparison_y(sig1, sig2, y_lims, use_lin):
    
    # Get the max signal bin and value       
    y_max = np.nanmax([np.nanmax(sig1[:int(sig1.size/2)]),
                       np.nanmax(sig2[:int(sig2.size/2)])])
    y_min = np.nanmin([np.nanmin(sig1[int(sig1.size/2):]),
                       np.nanmin(sig2[int(sig2.size/2):])])


    # Get the signal upper limit
    if use_lin == False:
        if y_lims[-1] == None:
            if not np.isfinite(y_max) or y_max <= 0:
                y_ulim = 1
            else:
                y_ulim = 2. * y_max
        else:
            if y_lims[0] <= 0:
                print('-- Warning: intercomparison y axis upper limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_ulim = 1
            else:
                y_ulim =  y_lims[-1]
        
    # Get the vertical lower limit
    if use_lin == False:
        if y_lims[0] == None:
            if not np.isfinite(y_min) or y_min <= 0:
                y_llim = y_ulim * 1E-3
            else:
                y_llim = 0.5 * y_min
        else:
            if y_lims[0] <= 0:
                print('-- Warning: intercomparison y axis lower limit <= 0 although the scale is logarithmic. The limit has automatically been replaced')
                y_llim = 0.
            else:
                y_llim =  y_lims[0]
    
    # Get the y axis labels
    y_label = 'Norm. RC Signals [$m^{-1} sr^{-1}$]'
    
    return(y_llim, y_ulim, y_label)

def round_it(x, sig):
    
    if not np.isfinite(x):
        x = -999.
        sig = 3
        
    if x != 0.:
        x_out = np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)
    else:
        x_out = 0.
        
    return x_out