#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 13:34:47 2022

@author: nick
"""

import numpy as np
import xarray as xr
from scipy.stats import linregress, shapiro
# from statsmodels.stats.diagnostic import acorr_ljungbox

default_norm_region = [5., 6.]

def statistics(y1, y2, x, keyw_args = {}, cross_check_type = 'apply_backwards', 
               cross_check_crit = 'above_negative_threshold', 
               cross_check_all_points = True, cancel_stats = []):
    
    """
    keyw_args: A dictionary of arguments that can be passed to configure the mask. These are namely:
        'fit_window' 
        'fit_window_step'
        'fit_mask_region'
        'rsem_threshold'
        'first_derivative_threshold'
        'second_derivative_threshold'
        'shapiro_wilk_threshold'
        'cross_criterion_threshold'
    
    cancel_stats: A list with the names of the statistical tests not to be considered. These are:
        'sem': Relative standard error of the mean
        'der': First derivative
        'sec': Second derivative (first derivative differences)
        'shp': Shapiro Wilk test
        'ccr': Cross criterion
        For example: cancel_stats = ['sem', 'ccr'] mean that the relative standard error of the mean and the cross criterion tests will not be considered in the mask
    
    cross_check_type: Configure where the cross criterion will be applied. Select one of the following sting values:
        'apply_backwards': Scan only regions that are below the normalization region (down to the beginning of the fit_mask_region)
        'apply_forwards': Scan only regions that are above the normalization region (uo to the end of the fit_mask_region)
        'apply_everywhere': Scan the whole profile (within the fit_mask region)
    
    cross_check_crit: Configure how the cross criterion will be applied. Select one of the following sting values:
        'above_negative_threshold': The minimum difference between the testing signal y1 and the reference signal y2 within the region where the cross criterion is applied must be smaller than the negative uncertainty of y1 multiplied by the cross_criterion_threshold 
        'below_positive_threshold': The maximum difference between the testing signal y1 and the reference signal y2 within the region where the cross criterion is applied must be larger than the positive uncertainty of y1 multiplied by the cross_criterion_threshold
        'between_both_thresholds': Both of the aforementioned thresholds will be applied. The absolute difference between the testing signal y1 and the reference signal y2 within the region where the cross criterion is applied must be smaller than the (positive) uncertainty of y1 multiplied by the cross_criterion_threshold
    """
    
    needed_keys = ['fit_mask_window', 
                   'fit_mask_window_step',
                   'fit_mask_region',
                   'rsem_threshold',
                   'first_derivative_threshold',
                   'second_derivative_threshold',
                   'shapiro_wilk_threshold',
                   'cross_criterion_threshold']
    
    default_values = {'fit_mask_window': [1., 8.], 
                      'fit_mask_window_step': 0.2,
                      'fit_mask_region': [2., 30.],
                      'rsem_threshold': 0.02,
                      'first_derivative_threshold': 2.,
                      'second_derivative_threshold': 2.,
                      'shapiro_wilk_threshold': 0.05,
                      'cross_criterion_threshold': 1.}
    
    for key in needed_keys:
        if key not in keyw_args.keys():
            keyw_args[key] = default_values[key]
            
    min_win = keyw_args['fit_mask_window'][0]
    max_win = keyw_args['fit_mask_window'][1]
    win_step = 1E-3 * keyw_args['fit_mask_window_step']
    mask_llim = keyw_args['fit_mask_region'][0]
    mask_ulim = keyw_args['fit_mask_region'][1]
    sem_threshold = keyw_args['rsem_threshold']
    der_threshold = keyw_args['first_derivative_threshold']
    sec_threshold = keyw_args['second_derivative_threshold'] 
    shp_threshold = keyw_args['shapiro_wilk_threshold']
    ccr_threshold = keyw_args['cross_criterion_threshold']
            
    resol = x[1] - x[0]
    
    mid = np.arange(mask_llim, mask_ulim + win_step / 2., win_step)
    
    win = np.arange(min_win, max_win + win_step / 2., win_step)
    
    vsem = np.nan * np.zeros((win.size,mid.size))
    vder = np.nan * np.zeros((win.size,mid.size))
    vsec = np.nan * np.zeros((win.size,mid.size))
    vshp = np.nan * np.zeros((win.size,mid.size))
    verr = np.nan * np.zeros((win.size,mid.size))
    vncf = np.nan * np.zeros((win.size,mid.size))
    vccr_min = np.nan * np.zeros((win.size,mid.size))
    vccr_max = np.nan * np.zeros((win.size,mid.size))
    

    for i in range(win.size):
    
        for j in range(mid.size):
      
            mask_x = (x >= mid[j] - win[i] / 2) & (x < mid[j] + win[i] / 2)
            
            mask_y = (y1[mask_x] == y1[mask_x]) & (y2[mask_x] == y2[mask_x])
            
            if np.sum(mask_y) >= 0.95 * win[i] / resol:
                
                y1_region = y1[mask_x][mask_y]
                y2_region = y2[mask_x][mask_y]
                
                x_region = x[mask_x][mask_y]

                # The normalization coefficient 
                vncf[i,j] = np.mean(y2_region) / np.mean(y1_region)
                
                # Calculate the residual from y2 (reference)
                residual_region = vncf[i,j] * y1_region - y2_region
                
                hbin = int(x_region.size / 2.)
                
                fit = linregress(x = x_region, y = residual_region)
    
                fit_lh = linregress(x = x_region[:hbin], y = residual_region[:hbin])
                
                fit_uh = linregress(x = x_region[hbin:], y = residual_region[hbin:])     
    
                # Calculate the standard error of the residual in the fit region 
                vsem[i,j] = np.nanstd(residual_region) / (np.sqrt(residual_region.size) * np.nanmean(y2_region))
                
                # Calculate the absolute relative standar error of the residual derivative
                vder[i,j] = np.abs(fit[0] / fit[4])
    
                # Calculate the absolute relative standard error of the first derivative difference of the residual between the first and second halves of the fit region
                vsec[i,j] = np.abs(fit_lh[0] - fit_uh[0]) / (fit_lh[4] + fit_uh[4])
                
                # Calculate the Shapiro-Wilk p value in the fitregion
                vshp[i,j] = shapiro(residual_region)[1]

                # Calculate the standard error of the signal in the fit region normalized 
                verr[i,j] = np.nanstd(y1_region) / (np.sqrt(y1_region.size) * np.nanmean(y1_region))
    
    for i in range(vncf.shape[0]):
    
        for j in range(1,vncf.shape[1]-1):
            
            if vncf[i,j] == vncf[i,j]:
      
                mask_x = (x >= mid[j] - win[i] / 2) & (x < mid[j] + win[i] / 2)
                    
                if cross_check_type == 'apply_backwards':
                    vccr_max[i,j] = np.max((1. - vncf[i,:j] / vncf[i,j]) / verr[i,:j])
                    vccr_min[i,j] = np.min((1. - vncf[i,:j] / vncf[i,j]) / verr[i,:j])
                elif cross_check_type == 'apply_forwards':
                    vccr_max[i,j] = np.max((1. - vncf[i,j:] / vncf[i,j]) / verr[i,j:])
                    vccr_min[i,j] = np.min((1. - vncf[i,j:] / vncf[i,j]) / verr[i,j:])
                elif cross_check_type == 'apply_everywhere':
                    vccr_max[i,j] = np.max((1. - vncf[i,:] / vncf[i,j]) / verr[i,:])
                    vccr_min[i,j] = np.min((1. - vncf[i,:] / vncf[i,j]) / verr[i,:])
                else:
                    raise Exception(f"-- Error: The provided cross_check_type {cross_check_type} is wrong. Please select one of: back, forth, both")
        
    # Calculate the standard error mask
    msem = (vsem <= sem_threshold)
    msem[msem != msem] = False
    if 'sem' in cancel_stats: msem = np.ones_like(msem, dtype = bool)
    
    # Calculate the derivative mask (pval < 0.05 means the slope is significant)
    mder = vder <= der_threshold
    mder[mder != mder] = False
    if 'der' in cancel_stats: mder = np.ones_like(mder, dtype = bool)

    # Check if the derivatives inside the 2 halves of the window are aslo not significant
    msec = vsec <= sec_threshold
    msec[msec != msec] = False
    if 'sec' in cancel_stats: msec = np.ones_like(msec, dtype = bool)

    # Check if the p value of the Saphiro-Wilkinson test is higher than 0.95 (<0.05 means not normal)
    mshp = vshp > shp_threshold
    mshp[mshp != mshp] = False
    if 'shp' in cancel_stats: mshp = np.ones_like(mshp, dtype = bool)
    
    if cross_check_crit == 'above_negative_threshold':
        mccr = vccr_min >= -ccr_threshold
    elif cross_check_type == 'below_positive_threshold':
        mccr = vccr_max <= ccr_threshold
    elif cross_check_type == 'between_both_thresholds':
        mccr = (vccr_min >= -ccr_threshold) & (vccr_max <= ccr_threshold)
    else:
        raise Exception(f"-- Error: The provided cross_check_crit {cross_check_crit} is wrong. Please select one of: above_negative_threshold, below_positive_threshold, between_both_thresholds")
    mccr[mccr != mccr] = False
    if 'ccr' in cancel_stats: mccr = np.ones_like(mccr, dtype = bool)

    from matplotlib import pyplot as plt
    plt.pcolormesh(vccr_min)
    plt.plot()
    
    mfit = (msem) & (mder) & (msec) & (mshp) & (mccr)
  
    # mtot = (msem) & ()
    # mtmp = tmp_mask(msem = msem, mder = mder, msec = msec, mshp = mshp,
    #                 cross_check_all_points = cross_check_all_points,
    #                 cancel_sem = cancel_sem, 
    #                 cancel_der = cancel_der, 
    #                 cancel_sec = cancel_sec,
    #                 cancel_shp = cancel_shp)
    
    # mccr = ccr_check(coef = coef, rerr = rerr, mtmp = mtmp, 
    #                  ccr_threshold = ccr_threshold, ccr_lim = ccr_lim, 
    #                  cross_check_type = cross_check_type,
    #                  cross_check_crit = cross_check_crit)
    
    # mfit = fit_mask(msem = msem, mder = mder, msec = msec, mshp = mshp, mccr = mccr,
    #                 cancel_sem = cancel_sem, 
    #                 cancel_der = cancel_der, 
    #                 cancel_sec = cancel_sec,
    #                 cancel_shp = cancel_shp, 
    #                 cancel_ccr = cancel_ccr)

    vncf = xr.DataArray(vncf, dims = ['window', 'lower_limit'],
                            coords = [win, mid])
        
    vsem = xr.DataArray(vsem, dims = ['window', 'lower_limit'],
                        coords = [win, mid])    
        
    vder = xr.DataArray(vder, dims = ['window', 'lower_limit'],
                        coords = [win, mid])
    
    vsec = xr.DataArray(vsec, dims = ['window', 'lower_limit'],
                        coords = [win, mid])
    
    vshp = xr.DataArray(vshp, dims = ['window', 'lower_limit'],
                        coords = [win, mid])

    vccr_min = xr.DataArray(vccr_min, dims = ['window', 'lower_limit'],
                            coords = [win, mid])

    vccr_max = xr.DataArray(vccr_max, dims = ['window', 'lower_limit'],
                            coords = [win, mid])
    
    msem = xr.DataArray(msem, dims = ['window', 'lower_limit'],
                        coords = [win, mid])
        
    mder = xr.DataArray(mder, dims = ['window', 'lower_limit'],
                        coords = [win, mid])

    msec = xr.DataArray(msec, dims = ['window', 'lower_limit'],
                        coords = [win, mid])
    
    mshp = xr.DataArray(mshp, dims = ['window', 'lower_limit'],
                        coords = [win, mid])

    mccr = xr.DataArray(mccr, dims = ['window', 'lower_limit'],
                        coords = [win, mid])
    
    mfit = xr.DataArray(mfit, dims = ['window', 'lower_limit'],
                        coords = [win, mid])
    

    
    masks = {'relative_sem' : msem,
             'first_derivative' : mder,
             'second_derivative' : msec,
             'shapiro_wilk' : mshp,
             'cross_criterion' : mccr,
             'total' : mfit}
    
    stats = {'relative_sem' : vsem,
             'first_derivative' : vder,
             'second_derivative' : vsec,
             'shapiro_wilk' : vshp,
             'cross_criterion_min' : vccr_min,
             'cross_criterion_max' : vccr_max,
             'normalization_factor' : vncf}
    
    return(stats, masks)

# def ccr_check(coef, rerr, mtmp, ccr_threshold, ccr_lim, cross_check_type, cross_check_crit):
    
#     mccr = np.zeros(coef.shape, dtype = bool)

#     if not np.isfinite(ccr_lim): 
#         for i in range(coef.shape[0]):
#             for j in range(coef.shape[1]): 
#                 if mtmp[i,j] == True:
#                     if cross_check_type == 'back':
#                         slicer = slice(None,j)
#                     elif cross_check_type == 'forth':
#                         slicer = slice(j,None)
#                     elif cross_check_type == 'both':
#                         slicer = slice(None,None)
#                     else:
#                         raise Exception(f"-- Error: The provided cross_check_type {cross_check_type} is wrong. Please select one of: back, forth, both")

#                     if cross_check_crit == 'min':
#                         crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] <= ccr_threshold * rerr[i,slicer]
#                     elif cross_check_type == 'max':
#                         crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] >= ccr_threshold * rerr[i,slicer]
#                     elif cross_check_type == 'both':
#                         crit = np.abs((coef[i,slicer] - coef[i,j]) / coef[i,j]) <= ccr_threshold * rerr[i,slicer]
#                     else:
#                         raise Exception(f"-- Error: The provided cross_check_crit {cross_check_crit} is wrong. Please select one of: min, max, both")

#                     mccr[i,j] = (crit[mtmp[i,slicer]]).all()

#     else:
#         for i in range(coef.shape[0]):
#             for j in range(coef.shape[1]): 
#                 if mtmp[i,j] == True:
#                     if cross_check_type == 'back':
#                         slicer = slice(None,j)
#                     elif cross_check_type == 'forth':
#                         slicer = slice(j,None)
#                     elif cross_check_type == 'both':
#                         slicer = slice(None,None)
#                     else:
#                         raise Exception(f"-- Error: The provided cross_check_type {cross_check_type} is wrong. Please select one of: back, forth, both")

#                     if cross_check_crit == 'min':
#                         crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] <= ccr_lim
#                     elif cross_check_type == 'max':
#                         crit = (coef[i,slicer] - coef[i,j]) / coef[i,j] >= ccr_lim
#                     elif cross_check_type == 'both':
#                         crit = np.abs((coef[i,slicer] - coef[i,j]) / coef[i,j]) <= np.abs(ccr_lim)
#                     else:
#                         raise Exception(f"-- Error: The provided cross_check_crit {cross_check_crit} is wrong. Please select one of: min, max, both")

#                     mccr[i,j] = (crit[mtmp[i,slicer]]).all()
                    
#     return(mccr)
         
# def tmp_mask(mder, msec, mshp, msem, cross_check_all_points,
#              cancel_sem, cancel_der, cancel_sec,
#              cancel_shp):

#     if not cancel_sem:
#         crit_sem = (msem == True)
#     else:
#         crit_sem = np.ones(msem.shape, dtype = bool)
        
#     if not cancel_der:
#         crit_der = (mder == True)
#     else:
#         crit_der = np.ones(mder.shape, dtype = bool)
        
#     if not cancel_sec:
#         crit_sec = (msec == True)
#     else:
#         crit_sec = np.ones(msec.shape, dtype = bool)
        
#     if not cancel_shp:
#         crit_shp = (mshp == True)
#     else:
#         crit_shp = np.ones(mshp.shape, dtype = bool)

#     if cross_check_all_points:
#         mtmp = np.ones(mder.shape, dtype = bool)
#     else:
#         mtmp = \
#             (crit_sem) & (crit_der) & (crit_sec) & (crit_shp)
    
#     return(mtmp)

# def fit_mask(msem, mder, msec, mshp, mccr,
#              cancel_sem, cancel_der, cancel_sec,
#              cancel_shp, cancel_ccr):

#     if not cancel_sem:
#         crit_sem = (msem == True)
#     else:
#         crit_sem = np.ones(msem.shape, dtype = bool)
        
#     if not cancel_der:
#         crit_der = (mder == True)
#     else:
#         crit_der = np.ones(mder.shape, dtype = bool)
        
#     if not cancel_sec:
#         crit_sec = (msec == True)
#     else:
#         crit_sec = np.ones(msec.shape, dtype = bool)
        
#     if not cancel_shp:
#         crit_shp = (mshp == True)
#     else:
#         crit_shp = np.ones(mshp.shape, dtype = bool)
        
#     if not cancel_ccr:
#         crit_ccr = (mccr == True)
#     else:
#         crit_ccr = np.ones(mccr.shape, dtype = bool)
        
#     mfit = (crit_sem) & (crit_der) & (crit_sec) & (crit_shp) & (crit_ccr)
    
#     return(mfit)
       
def scan(masks, prefered_range = 'far'):
    
    mfit = masks['total']
    lower_limit = mfit.lower_limit.copy()
    
    window = mfit.window.copy()
    
    # mask = ((masks['relative_sem'].values) & (masks['first_derivative'].values) & (masks['second_derivative'].values) &\
    #  (masks['shapiro_wilk'].values) & (masks['cross_criterion'].values))
    # from matplotlib import pyplot as plt
    print(masks['shapiro_wilk'])
    raise Exception
    # plt.pcolormesh(mask)
                    
    if mfit.any():
        # ulim = np.nan * mfit.copy()
        win_mid = (window / 2. + lower_limit).reshape('window', 'lower_limit')
        win_mid = win_mid.where(mfit == True)
        if prefered_range == 'far':
            
            # for i in range(ulim.window.size):
            #     ulim[i,:] = lower_limit + window[i]
            # ulim = ulim.where(mfit == True)
            idx = win_mid.argmax(dim = ('window', 'lower_limit'))
            
        elif prefered_range == 'near':
            # win_mid = (window / 2. + lower_limit).reshape('window', 'lower_limit')

            # for i in range(ulim.window.size):
            #     ulim[i,:] = lower_limit - window[i]
            # ulim = ulim.where(mfit == True)
            idx = win_mid.argmin(dim = ('window', 'lower_limit'))
            
        else:
            raise Exception(f"-- Error: The provided prefered_range {prefered_range} is wrong. Select one of: near, far")
             
        print(idx)
        raise Exception
        auto_norm_region = [lower_limit.loc[idx.lower_limit] - window[idx.window] / 2., 
                            lower_limit[idx.lower_limit] + window[idx.window] / 2.] 
    else:
        auto_norm_region = None
        idx = None
    
    return(auto_norm_region, idx)   

def norm_region_index(norm_region, idx, masks):
    
    if idx == None:
        mfit = masks['total']
        lower_limit = mfit.lower_limit.copy().values
        window = mfit.window.copy().values
        
        norm_region_mid = (norm_region[0] + norm_region[1]) / 2.
        norm_region_width = norm_region[1] - norm_region[0]
        
        lower_limit = np.where(lower_limit >= norm_region_mid)[0][0]
        window = np.where(window >= norm_region_width)[0][0]
    
    return(idx)

def metrics_norm_region(idx, stats, masks):
    
    stats_norm_region = {}
    masks_norm_region = {}
    
    if idx != None:
        for key in stats.keys():
            stats_norm_region[key] = stats[key][idx].values
        
        for key in masks.keys():
            masks_norm_region[key] = masks[key][idx].values

    return(stats_norm_region, masks_norm_region)
    
def select_norm_region(auto_norm_region, user_norm_region):
    
    if user_norm_region != None:
        norm_region = user_norm_region
        norm_region_flag = 'external'
    
    elif user_norm_region == None and auto_norm_region != None:
        norm_region = auto_norm_region
        norm_region_flag = 'auto'
    
    else:
        norm_region = default_norm_region
        norm_region_flag = 'default'
        
    return(norm_region, norm_region_flag)


