#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:41:20 2022

@author: nick
"""

import numpy as np

def sliding_average_2D(z_vals, y_vals, y_sm_lims, y_sm_win, expo = False):
    
    s_bin = np.where(y_vals >= y_sm_lims[0])[0][0]
    e_bin = np.where(y_vals <= y_sm_lims[-1])[0][-1]
    s_ihwin = int(1E-3 * y_sm_win[0]  / (2. * (y_vals[1] - y_vals[0])))
    e_ihwin = int(1E-3 * y_sm_win[-1] / (2. * (y_vals[1] - y_vals[0])))

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
    
    for i in range(s_bin, e_bin + 1):
        
        if ihwins[i] > i:
            ihwin = i

        elif ihwins[i] > ihwins.size - i:
            ihwin = ihwins.size - i
        
        else:
            ihwin = ihwins[i]
            
        if ihwin > 0 and not np.isnan(ihwin):
            
            z_vals_sm[:,i] = np.nanmean(z_vals[:,i-int(ihwin):i+int(ihwin)+1], 
                                        axis = 1)
        
        if np.isnan(ihwin):
            z_vals_sm[:,i] = np.nan
    
    return(z_vals_sm)

def sliding_average_1D(y_vals, x_vals, x_sm_lims, x_sm_win, 
                       expo = False, err_type = 'sem'):

    s_bin = np.where(x_vals >= x_sm_lims[0])[0][0]
    e_bin = np.where(x_vals <= x_sm_lims[-1])[0][-1]
    s_ihwin = int(1E-3 * x_sm_win[0]  / (2.*(x_vals[1] - x_vals[0])))
    e_ihwin = int(1E-3 * x_sm_win[-1] / (2.*(x_vals[1] - x_vals[0])))
    
    ihwins = np.zeros(x_vals.size)
    
    y_vals = y_vals.copy()
    y_vals_sm = y_vals.copy()
    y_vals_err = np.nan * np.zeros(y_vals.shape)
    
    if expo:
        ihwins[s_bin:e_bin] = \
            np.logspace(np.log2(s_ihwin), np.log2(e_ihwin),
                        e_bin - s_bin, base = 2)
    else:
        ihwins[s_bin:e_bin] = \
            np.linspace(s_ihwin, e_ihwin, e_bin - s_bin)
   
    ihwins[e_bin:] = np.nan

    for i in range(s_bin, e_bin + 1):
        
        if ihwins[i] > i:
            ihwin = i

        elif ihwins[i] > ihwins.size - i:
            ihwin = ihwins.size - i
        
        else:
            ihwin = ihwins[i]

        if ihwin > 0 and not np.isnan(ihwin):
            y_vals_sm[i] = np.nanmean(y_vals[i-int(ihwin):i+int(ihwin)+1])
            
            if err_type == 'sem':
                y_vals_err[i] = np.nanstd(y_vals[i-int(ihwin):i+int(ihwin)+1]) /\
                    np.sqrt(y_vals[i-int(ihwin):i+int(ihwin)+1].size)
            elif err_type == 'std':
                y_vals_err[i] = np.nanstd(y_vals[i-int(ihwin):i+int(ihwin)+1])
            else:
                raise Exception('Error type provided is not understoud. Please select one of: sem, std')

    return(y_vals_sm, y_vals_err)

def sliding_average_1D_fast(y_vals, x_vals, x_sm_lims, x_sm_win, expo = None, err_type = 'sem'):
    
    win = int(1E-3 * x_sm_win  / (x_vals[1] - x_vals[0]))
    
    if win % 2 == 0: win = win + 1
    
    s_bin = np.where(x_vals >= x_sm_lims[0])[0][0]
    e_bin = np.where(x_vals <= x_sm_lims[1])[0][-1]
    
    buf = int(np.floor(win / 2.))

    s_buf = s_bin + buf
    e_buf = e_bin - buf
    
    y_avg = y_vals.copy()
    del y_vals
    
    y_sqr = np.power(y_avg, 2)
    y_err = np.nan * np.zeros(y_avg.shape)

    y_avg[s_buf:e_buf] = np.convolve(y_avg[s_bin:e_bin], np.ones(win), 'valid') / win
    y_sqr[s_buf:e_buf] = np.convolve(y_sqr[s_bin:e_bin], np.ones(win), 'valid') / win
    
    
    if err_type == 'sem':
        y_err[s_buf:e_buf] = np.sqrt(y_sqr[s_buf:e_buf] - np.power(y_avg[s_buf:e_buf], 2)) / np.sqrt(win)
    elif err_type == 'std':
        y_err[s_buf:e_buf] = np.sqrt(y_sqr[s_buf:e_buf] - np.power(y_avg[s_buf:e_buf], 2))
    else:
        raise Exception('Error type provided is not understoud. Please select one of: sem, std')

    del y_sqr
    
    return(y_avg, y_err)

def sliding_average_2D_fast(z_vals, y_vals, y_sm_lims, y_sm_win, expo = None):
    
    win = int(1E-3 * y_sm_win  / (y_vals[1] - y_vals[0]))
    
    if win % 2 == 0: win = win + 1
    
    s_bin = np.where(y_vals >= y_sm_lims[0])[0][0]
    e_bin = np.where(y_vals <= y_sm_lims[1])[0][-1]
    
    buf = int(np.floor(win / 2.))

    s_buf = s_bin + buf
    e_buf = e_bin - buf
    
    z_avg = z_vals.copy()
    del z_vals
    
    z_sqr = np.power(z_avg, 2)
    # z_err = np.nan * np.zeros(z_avg.shape)

    for j in range(z_avg.shape[0]):
        z_avg[j,s_buf:e_buf] = np.convolve(z_avg[j,s_bin:e_bin], np.ones(win), 'valid') / win
        # z_sqr[s_buf:e_buf,j] = np.convolve(z_sqr[j,s_bin:e_bin], np.ones(win), 'valid') / win
        
        
        # if err_type == 'sem':
        #     z_err[s_buf:e_buf,j] = np.sqrt(z_sqr[j,s_buf:e_buf] - np.power(z_avg[j,s_buf:e_buf], 2)) / np.sqrt(win)
        # elif err_type == 'std':
        #     z_err[s_buf:e_buf,j] = np.sqrt(z_sqr[j,s_buf:e_buf] - np.power(z_avg[j,s_buf:e_buf], 2))
        # else:
        #     raise Exception('Error type provided is not understoud. Please select one of: sem, std')

    del z_sqr
    
    return(z_avg)