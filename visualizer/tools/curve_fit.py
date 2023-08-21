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
    
    min_win = 0.5
    
    max_win = 4.
    
    step = 0.1
    
    llim = step * np.ceil(x[0] / step)
    
    ulim = step * np.floor(x[-1] / step) - max_win
    
    edg = np.arange(llim, ulim + step, step)
    
    win = np.arange(min_win, max_win + step, step)
    
    nder = np.nan * np.zeros((win.size,edg.size))
    mder = np.nan * np.zeros((win.size,edg.size))
    rerr = np.nan * np.zeros((win.size,edg.size))
    rsem = np.nan * np.zeros((win.size,edg.size))
    mshp = np.nan * np.zeros((win.size,edg.size))
    coef = np.nan * np.zeros((win.size,edg.size))
    mneg = np.nan * np.zeros((win.size,edg.size))
    
    for i in range(win.size):
    
        for j in range(edg.size):
      
            mask_x = (x >= edg[j]) & (x < edg[j] + win[i])
            
            y1_sel = y1[mask_x]

            y2_sel = y2[mask_x]
            
            # The normalization coefficient 
            coef[i,j] = np.mean(y2_sel) / np.mean(y1_sel)
            
            # Calculate the residual from the molecular atmosphere
            res_sel = coef[i,j] * y1_sel - y2_sel
            
            x_sel = x[mask_x]
            
            fit = linregress(x = x_sel, y = res_sel)

            # Calculate the relative signal standar error derivative
            nder[i,j] = fit[0] / fit[4]
            
            # Calculate the derivative mask (pval < 0.05 means the slope is significant)
            mder[i,j] = fit[3] > 0.05

            # Calculate the standard error
            rsem[i,j] = np.nanstd(res_sel) / np.sqrt(res_sel.size) / np.nanmean(y2_sel)
            
            # Calculate the mean derivative
            rerr[i,j] = np.nanstd(y1_sel) / np.sqrt(res_sel.size) / np.nanmean(y1_sel)

            # Calculate the p value of the Saphiro-Wilkinson test (<0.05 means not normal)
            mshp[i,j] = shapiro(res_sel)[1] > 0.05

        
    for i in range(coef.shape[0]):
        for j in range(coef.shape[1]):    
            mneg[i,j] = ((coef[i,:].copy() - coef[i,j].copy()) / coef[i,j].copy() <= rerr[i,:]).all()
            
    nder = xr.DataArray(nder, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    mder = xr.DataArray(mder, dims = ['window', 'lower_limit'],
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
        
    return(nder, mder, rerr, rsem, mshp, coef, mneg)

def scan(mder, rsem, mshp, mneg):
    
    mmol = (mder == True) & (mshp == True) & (rsem <= 0.02) & (mneg == True)
    
    if mmol.any():
    
        idx_min = rsem.where(mmol).argmin(dim = ['window', 'lower_limit'])
        
        ismol = True
    
    else:
        idx_min = rsem.argmin(dim = ['window', 'lower_limit'])

        ismol = False

    win = rsem[idx_min].window.values

    lower_limit = float(rsem[idx_min].lower_limit.values)
    
    norm_reg = [lower_limit, lower_limit + win]
    
    return(norm_reg, mmol, idx_min, ismol)       
    
    