#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 14:29:47 2025

@author: nikos
"""

import xarray as xr
import numpy as np
from datetime import datetime

def get_mid_time(times):
    
    xr_type = type(xr.DataArray())
    np_type = type(np.array([]))
    
    if type(times) == xr_type:
        times = times.values
    elif type(times) == np_type:
        pass
    else:
        raise Exception(f"-- Error: The times parameter must be a DataArray or a Numpy array. Type detected: {type(times)}")
    
    mid_time_np64 = times[0] + (times[-1] - times[0]) / 2.
    
    mid_time = np.datetime64(mid_time_np64, 'us').item()
    
    return(mid_time)

def date_time_wyoming(bnames):
    
    # dates = [name[:13].split('_')[0] for name in bnames]
    # times = [name[:13].split('_')[1] for name in bnames]
    
    datetime_str = [name[:13] for name in bnames]

    dts = np.array([datetime.strptime(dt_str,'%Y%m%d_%H%M') 
                    for dt_str in datetime_str])
    
    return(dts)

def date_time_ecmwf(bnames):
    
    datetime_str = [name[:12] for name in bnames]
    
    dts = np.array([datetime.strptime(dt_str,'%Y%m%d%H%M') 
                    for dt_str in datetime_str])
    
    # dates = np.array([dt.strftime('%Y%m%d') for dt in dts])
    # times = np.array([dt.strftime('%H%M') for dt in dts])
    
    return(dts)
    
def find_nearest_file(mtime, bnames, filetype = 'ecmwf'):
    
    if filetype == 'wyoming':
        dts = date_time_wyoming(bnames)
    elif filetype == 'ecmwf':
        dts = date_time_ecmwf(bnames)
    else:
        raise Exception("--Error: Filetype {filetype} not supported")
        
    delta_t = np.array([(dt - mtime).total_seconds() /3600. for dt in dts])

    ind_rs = np.argmin(np.abs(delta_t))
    
    selected_file = bnames[ind_rs]
    
    if not any(np.abs(delta_t) < 24):
        print(f"-- Warning: The nearest radiosonde in time {bnames[ind_rs]} was launched with a time difference of {np.round(delta_t[ind_rs],decimals=1)} hours with respect to the middle time of the measurement! Please provide a radiosond file with less than 18 hours temporal difference")
        rsonde_flag = "unusable"

    else:
        print(f'-- Selected radiosonde file: {bnames[ind_rs]}')
        rsonde_flag = "usable"
        
    return(selected_file, rsonde_flag)

def round_it(x, sig):
    
    if not np.isfinite(x):
        x = -999.
        sig = 3
        
    if x != 0.:
        x_out = np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)
    else:
        x_out = 0.
        
    return x_out
