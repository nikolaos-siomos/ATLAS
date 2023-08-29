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
from scipy.stats import linregress, shapiro

def stats(y1, y2, x):
    
    min_win = 1.
    
    max_win = 4.
    
    step = 0.1
    
    llim = 2.
    
    ulim = 16.
    
    edg = np.arange(llim, ulim + step, step)
    
    win = np.arange(min_win, max_win + step, step)
    
    nder = np.nan * np.zeros((win.size,edg.size))
    mder = np.nan * np.zeros((win.size,edg.size))
    cder = np.nan * np.zeros((win.size,edg.size))
    rerr = np.nan * np.zeros((win.size,edg.size))
    rsem = np.nan * np.zeros((win.size,edg.size))
    mshp = np.nan * np.zeros((win.size,edg.size))
    coef = np.nan * np.zeros((win.size,edg.size))
    mneg = np.nan * np.zeros((win.size,edg.size))
    
    for i in range(win.size):
    
        for j in range(edg.size):
      
            mask_x = (x >= edg[j]) & (x < edg[j] + win[i])
            
            # mask_x_lh = (x >= edg[j]) & (x < edg[j] + win[i] / 2.)

            # mask_x_uh = (x >= edg[j] + win[i] / 2.) & (x < edg[j] + win[i])
            
            y1_sel = y1[mask_x]

            y2_sel = y2[mask_x]
            
            # The normalization coefficient 
            coef[i,j] = np.mean(y2_sel) / np.mean(y1_sel)
            
            # Calculate the residual from the molecular atmosphere
            res_sel = coef[i,j] * y1_sel - y2_sel
            
            x_sel = x[mask_x]
            
            hbin = int(x_sel.size / 2.)
            
            fit = linregress(x = x_sel, y = res_sel)

            fit_lh = linregress(x = x_sel[:hbin], y = res_sel[:hbin])
            
            fit_uh = linregress(x = x_sel[hbin:], y = res_sel[hbin:])     

            # Calculate the relative signal standar error derivative
            nder[i,j] = fit[0] / fit[4]
            
            # Calculate the derivative mask (pval < 0.05 means the slope is significant)
            mder[i,j] = fit[3] > 0.05

            # Check if the derivatives inside the 2 halves of the window are aslo not significant
            cder[i,j] = (fit_lh[3] > 0.05) & (fit_uh[3] > 0.05)
            
            # Calculate the standard error
            rsem[i,j] = np.nanstd(res_sel) / np.sqrt(res_sel.size) / np.nanmean(y2_sel)
            
            # Calculate the mean derivative
            rerr[i,j] = np.nanstd(y1_sel) / np.sqrt(res_sel.size) / np.nanmean(y1_sel)

            # Calculate the p value of the Saphiro-Wilkinson test (<0.05 means not normal)
            mshp[i,j] = shapiro(res_sel)[1] > 0.05

        
    for i in range(coef.shape[0]):
        for j in range(coef.shape[1]):    
            mneg[i,j] = ((coef[i,:j].copy() - coef[i,j].copy()) / coef[i,j].copy() <= rerr[i,:j]).all()
            
    nder = xr.DataArray(nder, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    mder = xr.DataArray(mder, dims = ['window', 'lower_limit'],
                        coords = [win, edg])

    cder = xr.DataArray(cder, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    rerr = xr.DataArray(rerr, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    rsem = xr.DataArray(rsem, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    mshp = xr.DataArray(mshp, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    coef = xr.DataArray(coef, dims = ['window', 'lower_limit'],
                        coords = [win, edg])

    mneg = xr.DataArray(mneg, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
        
    return(nder, mder, cder, rerr, rsem, mshp, coef, mneg)

def scan(mder, cder, rsem, mshp, mneg, dflt_region, auto_fit):
    
    mmol = (mder == True) & (cder == True) & (mshp == True) & (rsem <= 0.02) & (mneg == True)
    
    lower_limit = mmol.lower_limit.values
    
    window = mmol.window.values
        
    if auto_fit == True and mmol.any():    
    
        ulim = np.nan * mmol.copy()
        
        for i in range(ulim.window.size):
            ulim[i,:] = lower_limit + window[i]

        ulim = ulim.where(mmol == True)
    
        idx = ulim.argmax(dim = ('window', 'lower_limit'))
        
        ismol = True
        
        norm_region = [lower_limit[idx['lower_limit']], 
                       lower_limit[idx['lower_limit']] + window[idx['window']]] 
        
    else:
        
        if not mmol.any():

            print("-- Warning: The Rayleigh fit normalization region could not be automatically retrieved")
        
        idx_lower_limit = \
            np.where(lower_limit >= dflt_region[0])[0][0]

        idx_window = \
            np.where(window >= dflt_region[1] - dflt_region[0])[0][0]
        
        idx = dict(window = idx_window, lower_limit = idx_lower_limit)
        
        ismol = False

        norm_region = dflt_region
    
    return(norm_region, mmol, idx, ismol)   
    
    # edge_lower_limit = np.nan * np.zeros(mmol.window.size)



    # if mmol.any():
    #     for i in range(mmol.window.size):
    #         if mmol[i,:].any():
    #             edge_lower_limit[i] = \
    #                 lower_limit[mmol[i,:].values == True][-1]
    #     edge_upper_limit = edge_lower_limit + window
    #     print(edge_upper_limit)
       
    #     idx_lower_limit = \
    #         np.where(lower_limit == edge_lower_limit[np.nanargmax(edge_upper_limit)])[0][0]
        
    #     print(idx_lower_limit)
    #     print(np.where(window == window[np.nanargmax(edge_upper_limit)])[0])
    #     idx_window = \
    #         np.where(window == window[np.nanargmax(edge_upper_limit)])[0][0]
        
    #     idx_min = dict(idx_window, idx_lower_limit)
        
    #     ismol = True
        
    #     norm_region = [lower_limit[idx_lower_limit], 
    #                    lower_limit[idx_lower_limit] + window[idx_window]] 
        
    # else:
        
    #     idx_lower_limit = \
    #         np.where(lower_limit >= dflt_region[0])[0][0][0]

    #     idx_window = \
    #         np.where(window >= dflt_region[1] - dflt_region[0])[0][0][0]
        
    #     idx_min = dict(idx_window, idx_lower_limit)

    #     ismol = False

    #     norm_region = dflt_region

