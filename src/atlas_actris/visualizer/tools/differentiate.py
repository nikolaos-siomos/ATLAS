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

def region(sig, x_vals, region, axis):

    diff_height = (region[1] + region[0]) / 2.
    
    hwin = (region[1] - region[0]) / 2.
    
    # Get the reference height bin    
    diff_bin = get_diff_bin(x_vals = x_vals, diff_height = diff_height)

    # Get the reference window in bins    
    hwin_bin = get_hwin_bin(x_vals = x_vals, hwin = hwin)
    
    sig_sel = choose_from_axis(sig, axis = 0, 
                               start = diff_bin - hwin_bin, 
                               stop = diff_bin + hwin_bin + 1)  

    x_sel = choose_from_axis(x_vals, axis = 0, 
                               start = diff_bin - hwin_bin, 
                               stop = diff_bin + hwin_bin + 1) 
    
    res = stats.linregress(x = x_sel, y = sig_sel)
        
    # Calculate the mean derivative
    mder = res[0]
    
    # Calculate the derivative error
    sder = res[4]

    # Calculate the derivative pval (<0.05 means the slope is significant)
    pval = res[3]
    
    return(mder, sder, pval)

def scan(sig, x_vals):
    
    min_win = 0.5
    
    max_win = 4.
    
    step = 0.1
    
    llim = step * np.ceil(x_vals[0] / step)
    
    ulim = step * np.floor(x_vals[-1] / step) - max_win
    
    edg = np.arange(llim, ulim + step, step)
    
    win = np.arange(min_win, max_win + step, step)
    
    slope = np.nan * np.zeros((win.size,edg.size))
    sderr = np.nan * np.zeros((win.size,edg.size))
    p_val = np.nan * np.zeros((win.size,edg.size))
    
    for i in range(win.size):
    
        for j in range(edg.size):
      
            mask_x = (x_vals >= edg[j]) & (x_vals < edg[j] + win[i])
            
            sig_sel = sig[mask_x]
            
            x_sel = x_vals[mask_x]
            
            res = stats.linregress(x = x_sel, y = sig_sel)

            res = stats.linregress(x = x_sel, y = sig_sel)
        
            # Calculate the mean derivative
            slope[i,j] = res[0]
            
            # Calculate the derivative error
            sderr[i,j] = res[4]
        
            # Calculate the derivative pval (<0.05 means the slope is significant)
            p_val[i,j] = res[3]
            
    slope = xr.DataArray(slope, dims = ['window', 'lower_limit'],
                         coords = [win, edg])
    
    sderr = xr.DataArray(sderr, dims = ['window', 'lower_limit'],
                         coords = [win, edg])
    
    p_val = xr.DataArray(p_val, dims = ['window', 'lower_limit'],
                         coords = [win, edg])
    
    return(slope, sderr, p_val)

def get_diff_bin(x_vals, diff_height):

    if diff_height < x_vals[0]:
        raise Exception(f'-- Error: Normalization height/distance  is too low ({diff_height}km) ' +
                        f'while the signal starts at {x_vals[0]}km')
        
    elif diff_height > x_vals[-1]:
        raise Exception(f'-- Error: Reference height/distance is too high ({diff_height}km) ' +
                        f'while the signal ends at {x_vals[-1]}km')
    else:
        diff_bin = np.where(x_vals >= diff_height)[0][0] 
        
    return(diff_bin)

def get_hwin_bin(x_vals, hwin):

    if hwin < (x_vals[1] - x_vals[0]):
        raise Exception(f'-- Error: The half reference window provided ({hwin}m) is ' +
                        'smaller than the signal vertical step')
        
    else:
        hwin_bin = int(hwin / (x_vals[1] - x_vals[0]))
        
    return(hwin_bin)

def choose_from_axis(a, axis, start, stop):

    if axis <= a.ndim:
    
        s = [slice(None) for i in range(a.ndim)]
        
        s[axis] = slice(start, stop)
        
        s = tuple(s)

    else:
        raise Exception('-- Error: The provided axis index is larger than the number ' +
                        'of the axises of the array')

    return a[s]