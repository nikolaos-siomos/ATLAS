#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:56:30 2026

@author: nikos
"""


import numpy as np

def convert_m_to_km(x):
   
    return 1E-3 * x

def get_x_label(vertical_scale):
   
    # Convert meters to kilometers and select ranges or heights for the x axis depending on the use_dis value 
    if vertical_scale == 'range':
        x_label = "Range from the lidar [km]"
        
    elif vertical_scale == 'height_agl':
        x_label = "Height [km asl]"
        
    elif vertical_scale == 'height_asl':
        x_label = "Height [km agl]"
        
    else:
        raise Exception("Vertical scale '{vertical_scale}' not supported")

    return x_label 

def multiply_y_values(sig, sig_err, coef):
    
    # Multiply the   
    y_vals  = coef * sig.copy()
    
    y_errs = coef * sig_err.copy()
    
    return(y_vals, y_errs)

def get_rayleigh_fit_y_limits(y1_vals, y2_vals, y_lims, wavelength, use_lin_scale):
    """
    Determine y-axis limits.

    Parameters
    ----------
    y1_vals : array-like
        Values used to estimate the upper y limit.
    y2_vals : array-like
        Values used to estimate the lower y limit.
    y_lims : list
        Either [] for automatic limits, or [lower, upper] for manual limits.
    wavelength : float
        Signal wavelength.
    use_lin_scale : bool
        If False, sanity checks are applied for logarithmic scale.

    Returns
    -------
    list
        [lower_limit, upper_limit]
    """

    scale_f = wavelength / 355.0
    scat_ratio_f = 2.5

    # Automatic limits
    if len(y_lims) == 0:
        auto = True
        y_max = np.nanmax(y1_vals)
        y_min = np.nanmin(y2_vals)

        y_llim = y_min / 2.0
        y_ulim = scat_ratio_f * scale_f * y_max

    # Manual limits
    else:
        auto = False
        y_llim = y_lims[0]
        y_ulim = y_lims[-1]

    # Sanity checks for logarithmic scale
    if not use_lin_scale:

        if y_llim <= 0:
            if not auto:
                print(
                    "-- Warning: y-axis lower limit <= 0 although the scale is "
                    "logarithmic. The limit has automatically been replaced."
                )
            y_llim = 1e-6

        if y_ulim <= 0:
            if not auto:
                print(
                    "-- Warning: y-axis upper limit <= 0 although the scale is "
                    "logarithmic. The limit has automatically been replaced."
                )
            y_ulim = 1.0

    return [y_llim, y_ulim]

def get_quicklook_x_axis(t_lims, t_tick, time):

    # Treat None or [] as "use full time range"
    if t_lims is None or len(t_lims) == 0:
        t_lims = [None, None]

    # If only one limit is provided, complete the pair
    if len(t_lims) == 1:
        t_lims = [t_lims[0], None]

    # Get the lowest time based on t_lims
    if t_lims[0] is None:
        ltime = time[0]
    else:
        lmins = int(str(t_lims[0])[:2]) * 60 + int(str(t_lims[0])[2:])
        ltime = time.astype("datetime64[D]")[0] + np.timedelta64(lmins, "m")

    # Get the highest time based on t_lims
    if t_lims[-1] is None:
        utime = time[-1]
    else:
        umins = int(str(t_lims[-1])[:2]) * 60 + int(str(t_lims[-1])[2:])
        utime = time.astype("datetime64[D]")[-1] + np.timedelta64(umins, "m")

    # Get indices inside selected time range
    valid = np.where((time >= ltime) & (time <= utime))[0]

    if len(valid) == 0:
        raise ValueError(
            f"No time values found within t_lims={t_lims}. "
            "Check that t_lims uses HHMM format."
        )

    x_lbin = valid[0]
    x_ubin = valid[-1]

    # Keep full time axis for plotting
    t_vals = time

    # Calculate x_tick if not provided
    if time.size / 15. > 10.:
        x_tick = np.round(time.size / 15., decimals=-1)
    else:
        x_tick = np.round(time.size / 15., decimals=0)

    if x_tick == 0:
        x_tick = 1.

    # Calculate t_tick in minutes if not provided
    if t_tick is None:
        mins = (
            (t_vals[x_ubin] - t_vals[x_lbin]).astype("timedelta64[m]")
            / np.timedelta64(1, "m")
        )

        if mins < 5:
            t_tick = 1.
        elif mins < 20:
            t_tick = 2.
        elif mins < 40.:
            t_tick = 4.
        elif mins < 120.:
            t_tick = 10.
        elif mins < 240.:
            t_tick = 20.
        elif mins < 480.:
            t_tick = 30.
        else:
            t_tick = 60.

    return x_lbin, x_ubin, x_tick, t_vals, t_tick

def get_quicklook_y_axis(heights, ranges, y_lims, use_dis):

    # Use Height or range above the lidar for the y axis  
    if use_dis:
        y_vals = 1E-3 * ranges
        
        y_label = 'Range above the lidar [km]'
        
    else:
        y_vals = 1E-3 * heights       
        
        y_label = 'Height above the lidar [km]'

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

def get_quicklook_z_axis(z_vals, y_vals, z_lims, use_log, z_max_zone, z_min_zone):

    # Treat None or [] as "auto limits"
    if z_lims is None or len(z_lims) == 0:
        z_lims = [None, None]

    if len(z_lims) == 1:
        z_lims = [z_lims[0], None]

    # Get the max signal bin and value
    mask_max_zone = (y_vals >= z_max_zone[0]) & (y_vals <= z_max_zone[1])

    z_vals_sm = np.nanmean(z_vals[:, mask_max_zone], axis=0)

    z_max = round_it(np.nanmax(z_vals_sm), 1)

    # Avoid division by zero
    if z_max == 0 or np.isnan(z_max):
        z_max = 1.0

    # Normalize with the max
    z_vals = z_vals / z_max

    # Get the signal upper and lower limits
    z_llim = z_lims[0]
    z_ulim = z_lims[-1]

    # Get lower limit automatically
    if use_log and z_llim is None:

        mask_min_zone = (y_vals >= z_min_zone[0]) & (y_vals <= z_min_zone[1])

        z_vals_sm = np.nanmean(z_vals[:, mask_min_zone], axis=0)

        z_llim = round_it(np.nanmin(z_vals_sm), 2)

    elif not use_log and z_llim is None:
        z_llim = 0.0

    # Get upper limit automatically
    if z_ulim is None:
        z_ulim = round_it(np.nanmax(z_vals), 1)

    return z_llim, z_ulim, z_vals

