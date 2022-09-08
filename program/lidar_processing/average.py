#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 13:34:47 2022

@author: nick
"""

import sys
import numpy as np


def region(sig, x_vals, calibr, hwin, axis):

    # Get the reference height bin    
    calibr_bin = get_calibr_bin(x_vals = x_vals, calibr = calibr)

    # Get the reference window in bins    
    hwin_bin = get_hwin_bin(x_vals = x_vals, hwin = hwin)
        
    sig_sel = choose_from_axis(sig, axis, 
                               calibr_bin - hwin_bin, 
                               calibr_bin + hwin_bin + 1)

    avg = np.nanmean(sig_sel, axis = axis, keepdims = True)
        
    return(avg)

def get_calibr_bin(x_vals, calibr):

    if calibr < x_vals[0]:
        raise(f'-- Error: Calibration height/distance is too low ({calibr}km) ' +
              f'while the signal starts at {x_vals[0]}km')
        
    elif calibr > x_vals[-1]:
        raise(f'-- Error: Calibration height/distance is too high ({calibr}km) ' +
              f'while the signal ends at {x_vals[-1]}km')
    else:
        calibr_bin = np.where(x_vals >= calibr)[0][0] 
        
    return(calibr_bin)

def get_hwin_bin(x_vals, hwin):

    if hwin < (x_vals[1] - x_vals[0]):
        raise(f'-- Error: The half calibration window provided ({hwin}m) is ' +
              'smaller than the signal vertical step')
        
    else:
        hwin_bin = int(hwin * 1E-3 / (x_vals[1] - x_vals[0]))
        
    return(hwin_bin)

def choose_from_axis(a, axis, start, stop):

    if axis <= a.ndim:
    
        s = [slice(None) for i in range(a.ndim)]
        
        s[axis] = slice(start, stop)

    else:
        raise('-- Error: The provided axis index is larger than the number ' +
              'of the axises of the array')
    
    return a[s]