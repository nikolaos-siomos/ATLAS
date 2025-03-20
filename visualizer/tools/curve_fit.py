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

def stats(y1, y2, x, min_win, max_win, step, llim, ulim, 
          cross_check_type = 'back', cross_check_crit = 'min', cross_check_all_points = True, 
          rsem_lim = 0.02, der_fac = 1., shp_lim = 0.05, crc_lim = np.nan, crc_fac = 1., 
          cancel_sem = False, cancel_der = False, cancel_sec = False,
          cancel_shp = False, cancel_crc = False):
    
    edg = np.arange(llim - step / 2., ulim + step / 2., step)
    
    win = np.arange(min_win, max_win + step, step)
    
    nder = np.nan * np.zeros((win.size,edg.size))
    rerr = np.nan * np.zeros((win.size,edg.size))
    rsem = np.nan * np.zeros((win.size,edg.size))
    mder = np.nan * np.zeros((win.size,edg.size), dtype = bool)
    msec = np.nan * np.zeros((win.size,edg.size), dtype = bool)
    mshp = np.nan * np.zeros((win.size,edg.size), dtype = bool)
    coef = np.nan * np.zeros((win.size,edg.size))
    
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
            
            # # Calculate the derivative mask (pval < 0.05 means the slope is significant)
            # mder[i,j] = fit[3] > der_lim

            # # Check if the derivatives inside the 2 halves of the window are aslo not significant
            # msec[i,j] = (fit_lh[3] > der_lim) & (fit_uh[3] > der_lim)
            

            # Calculate the derivative mask (pval < 0.05 means the slope is significant)
            mder[i,j] = np.abs(fit[0] / fit[4]) <= der_fac

            # Check if the derivatives inside the 2 halves of the window are aslo not significant
            msec[i,j] = np.abs(fit_lh[0] - fit_uh[0]) / (fit_lh[4] + fit_uh[4]) <= der_fac
            
            # Calculate the standard error
            rsem[i,j] = np.nanstd(res_sel) / np.sqrt(res_sel.size) / np.nanmean(y2_sel)
            
            # Calculate the mean derivative
            rerr[i,j] = np.nanstd(y1_sel) / np.sqrt(res_sel.size) / np.nanmean(y1_sel)

            # Calculate the p value of the Saphiro-Wilkinson test (<0.05 means not normal)
            mshp[i,j] = shapiro(res_sel)[1] > shp_lim

    # Calculate the standard error mask
    msem = (rsem <= rsem_lim)
    
    mtmp = tmp_mask(msem = msem, mder = mder, msec = msec, mshp = mshp,
                    cross_check_all_points = cross_check_all_points,
                    cancel_sem = cancel_sem, 
                    cancel_der = cancel_der, 
                    cancel_sec = cancel_sec,
                    cancel_shp = cancel_shp)
    
    mcrc = crc_check(coef = coef, rerr = rerr, mtmp = mtmp, 
                     crc_fac = crc_fac, crc_lim = crc_lim, 
                     cross_check_type = cross_check_type,
                     cross_check_crit = cross_check_crit)
    
    mfit = fit_mask(msem = msem, mder = mder, msec = msec, mshp = mshp, mcrc = mcrc,
                    cancel_sem = cancel_sem, 
                    cancel_der = cancel_der, 
                    cancel_sec = cancel_sec,
                    cancel_shp = cancel_shp, 
                    cancel_crc = cancel_crc)
        
    rsem = xr.DataArray(rsem, dims = ['window', 'lower_limit'],
                        coords = [win, edg])    
        
    nder = xr.DataArray(nder, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    msem = xr.DataArray(msem, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
        
    mder = xr.DataArray(mder, dims = ['window', 'lower_limit'],
                        coords = [win, edg])

    msec = xr.DataArray(msec, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    mshp = xr.DataArray(mshp, dims = ['window', 'lower_limit'],
                        coords = [win, edg])

    mcrc = xr.DataArray(mcrc, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    mfit = xr.DataArray(mfit, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    coef = xr.DataArray(coef, dims = ['window', 'lower_limit'],
                        coords = [win, edg])
    
    return(rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef)

def crc_check(coef, rerr, mtmp, crc_fac, crc_lim, cross_check_type, cross_check_crit):
    
    mcrc = np.zeros(coef.shape, dtype = bool)

    if not np.isfinite(crc_lim): 
        for i in range(coef.shape[0]):
            for j in range(coef.shape[1]): 
                if mtmp[i,j] == True:
                    if cross_check_type == 'back':
                        slicer = slice(None,j)
                    elif cross_check_type == 'forth':
                        slicer = slice(j,None)
                    elif cross_check_type == 'both':
                        slicer = slice(None,None)
                    else:
                        raise Exception(f"-- Error: The provided cross_check_type {cross_check_type} is wrong. Please select one of: back, forth, both")

                    if cross_check_crit == 'min':
                        crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] <= crc_fac * rerr[i,slicer]
                    elif cross_check_type == 'max':
                        crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] >= crc_fac * rerr[i,slicer]
                    elif cross_check_type == 'both':
                        crit = np.abs((coef[i,slicer] - coef[i,j]) / coef[i,j]) <= crc_fac * rerr[i,slicer]
                    else:
                        raise Exception(f"-- Error: The provided cross_check_crit {cross_check_crit} is wrong. Please select one of: min, max, both")

                    mcrc[i,j] = (crit[mtmp[i,slicer]]).all()

    else:
        for i in range(coef.shape[0]):
            for j in range(coef.shape[1]): 
                if mtmp[i,j] == True:
                    if cross_check_type == 'back':
                        slicer = slice(None,j)
                    elif cross_check_type == 'forth':
                        slicer = slice(j,None)
                    elif cross_check_type == 'both':
                        slicer = slice(None,None)
                    else:
                        raise Exception(f"-- Error: The provided cross_check_type {cross_check_type} is wrong. Please select one of: back, forth, both")

                    if cross_check_crit == 'min':
                        crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] <= crc_lim
                    elif cross_check_type == 'max':
                        crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] >= crc_lim
                    elif cross_check_type == 'both':
                        crit = np.abs((coef[i,slicer] - coef[i,j]) / coef[i,j]) <= np.abs(crc_lim)
                    else:
                        raise Exception(f"-- Error: The provided cross_check_crit {cross_check_crit} is wrong. Please select one of: min, max, both")

                    mcrc[i,j] = (crit[mtmp[i,slicer]]).all()
                    
    return(mcrc)
         
def tmp_mask(mder, msec, mshp, msem, cross_check_all_points,
             cancel_sem, cancel_der, cancel_sec,
             cancel_shp):

    if not cancel_sem:
        crit_sem = (msem == True)
    else:
        crit_sem = np.ones(msem.shape, dtype = bool)
        
    if not cancel_der:
        crit_der = (mder == True)
    else:
        crit_der = np.ones(mder.shape, dtype = bool)
        
    if not cancel_sec:
        crit_sec = (msec == True)
    else:
        crit_sec = np.ones(msec.shape, dtype = bool)
        
    if not cancel_shp:
        crit_shp = (mshp == True)
    else:
        crit_shp = np.ones(mshp.shape, dtype = bool)

    if cross_check_all_points:
        mtmp = np.ones(mder.shape, dtype = bool)
    else:
        mtmp = \
            (crit_sem) & (crit_der) & (crit_sec) & (crit_shp)
    
    return(mtmp)

def fit_mask(msem, mder, msec, mshp, mcrc,
             cancel_sem, cancel_der, cancel_sec,
             cancel_shp, cancel_crc):

    if not cancel_sem:
        crit_sem = (msem == True)
    else:
        crit_sem = np.ones(msem.shape, dtype = bool)
        
    if not cancel_der:
        crit_der = (mder == True)
    else:
        crit_der = np.ones(mder.shape, dtype = bool)
        
    if not cancel_sec:
        crit_sec = (msec == True)
    else:
        crit_sec = np.ones(msec.shape, dtype = bool)
        
    if not cancel_shp:
        crit_shp = (mshp == True)
    else:
        crit_shp = np.ones(mshp.shape, dtype = bool)
        
    if not cancel_crc:
        crit_crc = (mcrc == True)
    else:
        crit_crc = np.ones(mcrc.shape, dtype = bool)
        
    mfit = (crit_sem) & (crit_der) & (crit_sec) & (crit_shp) & (crit_crc)
    
    return(mfit)
       
def scan(mfit, dflt_region, auto_fit, prefered_range = 'far'):
    
    lower_limit = mfit.lower_limit.values
    
    window = mfit.window.values
        
    if auto_fit == True and mfit.any():    
    
        ulim = np.nan * mfit.copy()
    
        if prefered_range == 'far':
            for i in range(ulim.window.size):
                ulim[i,:] = lower_limit + window[i]
            ulim = ulim.where(mfit == True)
            idx = ulim.argmax(dim = ('window', 'lower_limit'))
            
        elif prefered_range == 'near':
            for i in range(ulim.window.size):
                ulim[i,:] = lower_limit - window[i]
            ulim = ulim.where(mfit == True)
            idx = ulim.argmin(dim = ('window', 'lower_limit'))
            
        else:
            raise Exception(f"-- Error: The provided prefered_range {prefered_range} is wrong. Select one of: near, far")
        
        fit = True
        
        norm_region = [lower_limit[idx['lower_limit']], 
                       lower_limit[idx['lower_limit']] + window[idx['window']]] 
        
    else:
        
        if not mfit.any():

            print("-- Warning: The normalization region could not be automatically retrieved")
        
        idx_lower_limit = \
            np.argmin(np.abs(lower_limit - dflt_region[0]))

        idx_window = \
            np.argmin(np.abs(window - dflt_region[1] + dflt_region[0]))
        
        idx = dict(window = idx_window, lower_limit = idx_lower_limit)
        
        fit = False

        norm_region = dflt_region
    
    return(norm_region, idx, fit)   


