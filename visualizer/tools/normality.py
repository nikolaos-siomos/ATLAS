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
from scipy.stats import shapiro

def scan(sig, x_vals):
    
    min_win = 0.5
    
    max_win = 4.
    
    step = 0.1
    
    llim = step * np.ceil(x_vals[0] / step)
    
    ulim = step * np.floor(x_vals[-1] / step) - max_win
    
    edg = np.arange(llim, ulim + step, step)
    
    win = np.arange(min_win, max_win + step, step)
    
    p_val = np.nan * np.zeros((win.size,edg.size))
    
    for i in range(win.size):

        for j in range(edg.size):
      
            mask_x = (x_vals >= edg[j]) & (x_vals < edg[j] + win[i])
            
            sig_sel = sig[mask_x]

            # Calculate the mean value
            p_val[i,j] = shapiro(sig_sel)[1]

    p_val = xr.DataArray(p_val, dims = ['window', 'lower_limit'],
                         coords = [win, edg])
    
    return(p_val)