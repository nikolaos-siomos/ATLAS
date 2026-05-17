#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 13:34:47 2022

@author: nick
"""

import sys
import numpy as np
from scipy import stats
import xarray as xr

def region(sig, x_vals, region, axis, squeeze = False):

    avg_height = (region[1] + region[0]) / 2.
    
    hwin = (region[1] - region[0]) / 2.
    
    # Get the reference height bin    
    avg_bin = get_avg_bin(x_vals = x_vals, avg_height = avg_height)

    # Get the reference window in bins    
    hwin_bin = get_hwin_bin(x_vals = x_vals, hwin = hwin)
        
    sig_sel = choose_from_axis(sig, axis, 
                               avg_bin - hwin_bin, 
                               avg_bin + hwin_bin + 1)

    if squeeze:
        avg = np.nanmean(sig_sel, axis = axis)
        std = np.nanstd(sig_sel, axis = axis)
        sem = std/np.sqrt(sig_sel.size)
    else:
        avg = np.nanmean(sig_sel, axis = axis, keepdims = True)
        std = np.nanstd(sig_sel, axis = axis, keepdims = True)
        sem = std/np.sqrt(sig_sel.size)
        
    return(avg, std, sem)

def scan(sig, atb, x_vals):
    
    min_win = 0.5
    
    max_win = 4.
    
    step = 0.1
    
    llim = step * np.ceil(x_vals[0] / step)
    
    ulim = step * np.floor(x_vals[-1] / step) - max_win
    
    edg = np.arange(llim, ulim + step, step)
    
    win = np.arange(min_win, max_win + step, step)
    
    avg_sig = np.nan * np.zeros((win.size,edg.size))
    avg_atb = np.nan * np.zeros((win.size,edg.size))
    std = np.nan * np.zeros((win.size,edg.size))
    sem = np.nan * np.zeros((win.size,edg.size))
    
    for i in range(win.size):

        for j in range(edg.size):
      
            mask_x = (x_vals >= edg[j]) & (x_vals < edg[j] + win[i])
            
            sig_sel = sig[mask_x]

            atb_sel = atb[mask_x]
            
            x_sel = x_vals[mask_x]
        
            # Calculate the mean value
            avg_sig[i,j] = np.nanmean(sig_sel)

            # Calculate the mean value
            avg_atb[i,j] = np.nanmean(atb_sel)
            
            # Calculate the standard deviation
            std[i,j] = np.nanstd(sig_sel)

            # Calculate the standard error
            sem[i,j] = np.nanstd(sig_sel) / np.sqrt(sig_sel.size)            
            
    avg_sig = xr.DataArray(avg_sig, dims = ['window', 'lower_limit'],
                           coords = [win, edg])

    avg_atb = xr.DataArray(avg_atb, dims = ['window', 'lower_limit'],
                           coords = [win, edg])
    
    std = xr.DataArray(std, dims = ['window', 'lower_limit'],
                       coords = [win, edg])
    
    sem = xr.DataArray(np.abs(sem/avg_atb), dims = ['window', 'lower_limit'],
                       coords = [win, edg])
    
    
    return(avg_sig, std, sem)

def get_avg_bin(x_vals, avg_height):

    if avg_height < x_vals[0]:
        raise Exception('-- Error: The height/distance provided for averaging is too low '+\
                        f'({avg_height}km) while the signal starts at {x_vals[0]}km')
        
    elif avg_height > x_vals[-1]:
        raise Exception('-- Error: Calibration height/distance provided for averaging is too high ' +
                        f'({avg_height}km) while the signal ends at {x_vals[-1]}km')
    else:
        avg_bin = np.where(x_vals >= avg_height)[0][0] 
        
    return(avg_bin)

def get_hwin_bin(x_vals, hwin):

    if hwin < (x_vals[1] - x_vals[0]):
        raise Exception('-- Error: The half calibration window provided '+\
                        f'({hwin}m) is smaller than the signal vertical step')
        
    else:
        hwin_bin = int(hwin / (x_vals[1] - x_vals[0]))
        
    return(hwin_bin)

def choose_from_axis(a, axis, start, stop):

    if axis <= a.ndim - 1:
    
        s = [slice(None) for i in range(a.ndim)]
        
        s[axis] = slice(start, stop)

    else:
        raise Exception('-- Error: The provided axis index is larger than '+
                        'the number of the axises of the array')
    
    return a[tuple(s)]