#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:41:20 2022

@author: nick
"""

import numpy as np

def sliding_average_2D(z_vals, y_vals, y_sm_lims, y_sm_hwin, expo):
    
    s_bin = np.where(y_vals > y_sm_lims[0])[0][0]
    e_bin = np.where(y_vals < y_sm_lims[-1])[0][-1]
    s_ihwin = int(1E-3 * y_sm_hwin[0]  / (y_vals[1] - y_vals[0]))
    e_ihwin = int(1E-3 * y_sm_hwin[-1] / (y_vals[1] - y_vals[0]))
    
    ihwins = np.zeros(y_vals.size)

    z_vals = z_vals.copy()
    z_vals_sm = z_vals.copy()
    
    if expo:
        ihwins[s_bin:e_bin] = \
            np.logspace(np.log2(s_ihwin), np.log2(e_ihwin),
                        e_bin - s_bin, base = 2)
    else:
        ihwins[s_bin:e_bin] = \
            np.linspace(s_ihwin, e_ihwin, e_bin - s_bin)
   
    ihwins[e_bin:] = np.nan
    
    for i in range(ihwins.size):
        
        if ihwins[i] > i:
            ihwin = i

        elif ihwins[i] > ihwins.size - i:
            ihwin = ihwins.size - i
        
        else:
            ihwin = ihwins[i]
            
        if ihwin > 0 and not np.isnan(ihwin):
            
            z_vals_sm[i,:] = np.nanmean(z_vals[i-int(ihwin):i+int(ihwin)+1,:], 
                                        axis = 0)
        
        if np.isnan(ihwin):
            z_vals_sm[i,:] = np.nan
    
    return(z_vals_sm)

def sliding_average_1D(y_vals, x_vals, x_sm_lims, x_sm_hwin, expo):
    
    s_bin = np.where(x_vals > x_sm_lims[0])[0][0]
    e_bin = np.where(x_vals < x_sm_lims[-1])[0][-1]
    s_ihwin = int(1E-3 * x_sm_hwin[0]  / (x_vals[1] - x_vals[0]))
    e_ihwin = int(1E-3 * x_sm_hwin[-1] / (x_vals[1] - x_vals[0]))
    
    ihwins = np.zeros(x_vals.size)
    
    y_vals = y_vals.copy()
    y_vals_sm = y_vals.copy()
    y_vals_sem = np.nan * np.zeros(y_vals.shape)
    
    if expo:
        ihwins[s_bin:e_bin] = \
            np.logspace(np.log2(s_ihwin), np.log2(e_ihwin),
                        e_bin - s_bin, base = 2)
    else:
        ihwins[s_bin:e_bin] = \
            np.linspace(s_ihwin, e_ihwin, e_bin - s_bin)
   
    ihwins[e_bin:] = np.nan

    for i in range(ihwins.size):
        
        if ihwins[i] > i:
            ihwin = i

        elif ihwins[i] > ihwins.size - i:
            ihwin = ihwins.size - i
        
        else:
            ihwin = ihwins[i]

        if ihwin > 0 and not np.isnan(ihwin):
            y_vals_sm[i] = np.nanmean(y_vals[i-int(ihwin):i+int(ihwin)+1])
            
            y_vals_sem[i] = np.nanstd(y_vals[i-int(ihwin):i+int(ihwin)+1]) /\
                np.sqrt(2. * ihwin + 1.)

        
        if np.isnan(ihwin):
            y_vals_sm[i] = np.nan
    
    return(y_vals_sm, y_vals_sem)