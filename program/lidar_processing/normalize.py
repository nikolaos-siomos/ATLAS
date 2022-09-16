#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 13:34:47 2022

@author: nick
"""

import sys
import numpy as np


def to_a_point(sig, sig_b, x_vals, norm, hwin, axis, axis_b = 0):

    # Get the reference height bin    
    norm_bin = get_norm_bin(x_vals = x_vals, norm = norm)

    # Get the reference window in bins    
    hwin_bin = get_hwin_bin(x_vals = x_vals, hwin = hwin)
        
    sig_sel = choose_from_axis(sig, axis, 
                               norm_bin - hwin_bin, 
                               norm_bin + hwin_bin + 1)

    sig_b_sel = choose_from_axis(sig_b, axis_b, 
                               norm_bin - hwin_bin, 
                               norm_bin + hwin_bin + 1)    

    sig_m = np.nanmean(sig_sel, axis = axis, keepdims=True)
    
    sig_b_m = np.nanmean(sig_b_sel, axis = axis_b, keepdims=True)
    
    norm_coef = sig_b_m / sig_m 

    # sig_n = sig * norm_coef

    return(norm_coef)

def get_norm_bin(x_vals, norm):

    if norm < x_vals[0]:
        raise(f'-- Error: Normalization height/distance  is too low ({norm}km) ' +
              f'while the signal starts at {x_vals[0]}km')
        
    elif norm > x_vals[-1]:
        raise(f'-- Error: Reference height/distance is too high ({norm}km) ' +
              f'while the signal ends at {x_vals[-1]}km')
    else:
        norm_bin = np.where(x_vals >= norm)[0][0] 
        
    return(norm_bin)

def get_hwin_bin(x_vals, hwin):

    if hwin < (x_vals[1] - x_vals[0]):
        raise(f'-- Error: The half reference window provided ({hwin}m) is ' +
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

def add_axis(a, axis):
    
    if axis <= a.ndm + 1:
        s = [slice(None) for i in range(a.ndim + 1)]
        
        s[axis] = np.newaxis
    else:
        raise('-- Error: The provided axis index is larger than the number ' +
              'of the axises of the array plus one')
    
    return a[s]