#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 13:34:47 2022

@author: nick
"""

import numpy as np
import xarray as xr
from scipy.stats import linregress, shapiro
from visualizer.smoothing import sliding_average_1D_fast as smooth_1D
from utils.error_classes import CustomWarning

default_norm_region = [4., 10.]

needed_keys = [
    'fit_mask_window', 
    'fit_mask_window_step',
    'fit_mask_region',
    'rsem_threshold',
    'first_derivative_threshold',
    'second_derivative_threshold',
    'shapiro_wilk_threshold',
    'cross_criterion_threshold',
    'durbin_watson_threshold',
    'residual_extinction_threshold'
    ]

default_values = {
    'fit_mask_window': [1., 8.], 
    'fit_mask_window_step': 0.2,
    'fit_mask_region': [4., 30.],
    'rsem_threshold': 0.02,
    'first_derivative_threshold': 2.,
    'second_derivative_threshold': 2.,
    'shapiro_wilk_threshold': 0.05,
    'cross_criterion_threshold': 1,
    'durbin_watson_threshold': 0.5,
    'residual_extinction_threshold': 10.
    }

allowed_ccr_directions = [
    "apply_backwards",
    "apply_forwards",
    "apply_everywhere"
    ]

allowed_ccr_types = [
    "above_negative_threshold",
    "below_positive_threshold",
    "between_both_thresholds"
    ]

def norm_region_index(norm_region, idx, masks):
    
    if idx == None:
        mfit = masks['total']
        middle_point = mfit.middle_point.copy().values
        window = mfit.window.copy().values
        
        norm_region_mid = (norm_region[0] + norm_region[1]) / 2.
        norm_region_width = norm_region[1] - norm_region[0]
        
        middle_point = np.where(middle_point >= norm_region_mid)[0][0]
        window = np.where(window >= norm_region_width)[0][0]
    
    return(idx)

def metrics_norm_region(idx, stats, masks):
    
    stats_norm_region = {}
    masks_norm_region = {}
    
    if idx != None:
        for key in stats.keys():
            stats_norm_region[key] = stats[key][idx[0],idx[1]].values
        
        for key in masks.keys():
            masks_norm_region[key] = masks[key][idx[0],idx[1]].values

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

def idx_from_region(norm_region, middle_point, window):

    middle_point_user = (norm_region[0] + norm_region[1]) / 2.
    
    window_user = norm_region[1] - norm_region[0]
    
    middle_point_user = np.argmin(np.abs(middle_point.values - middle_point_user))
    
    window_idx_user = np.argmin(np.abs(window.values - window_user))
    
    auto_idx = [window_idx_user, middle_point_user]
    
    return(auto_idx)

def durbin_watson(residuals):
    
    residuals = np.asarray(residuals)
    
    diff = np.diff(residuals)
    
    return np.sum(diff ** 2) / np.sum(residuals ** 2)
    
    
def cross_crit_x_mask(x, start, end, cursor, 
                      ccr_direction, allowed_ccr_directions):
    
    if ccr_direction == 'apply_backwards':
        mask_x = (x <= cursor) & (x >= start)
    elif ccr_direction == 'apply_forwards':
        mask_x = (x >= cursor) & (x <= end)
    elif ccr_direction == 'apply_everywhere':
        mask_x = (x >= start) & (x <= end)
    else:
        raise Exception(f"-- Error: The provided ccr_direction {ccr_direction} is wrong. Please select one of: {allowed_ccr_directions}")

    return(mask_x)

def calculate_ccr_values(c_n, S, S_m, sigma):
    
    vccr_neg = np.nanmax((S_m - c_n * S) / (c_n * sigma))
    vccr_pos = np.nanmax((c_n * S - S_m) / (c_n * sigma))
    
    return(vccr_neg, vccr_pos)

def calculate_ccr_arrays(x, mid, win, c_n, S, S_m, sigma,
                         ccr_direction, allowed_ccr_directions):
    
    vccr_neg = np.nan * np.zeros((win.size,mid.size))
    vccr_pos = np.nan * np.zeros((win.size,mid.size))
    
    for j in range(1,mid.size-1):
        mask_x = cross_crit_x_mask(
            x = x, 
            cursor = mid[j],
            start = mid[0],
            end = mid[-1],
            ccr_direction = ccr_direction, 
            allowed_ccr_directions = allowed_ccr_directions
            )
        
        for i in range(win.size):
            if c_n[i,j] == c_n[i,j]:
                vccr_neg[i,j], vccr_pos[i,j] = calculate_ccr_values(
                    c_n = c_n[i,j], 
                    S = S[mask_x], 
                    S_m = S_m[mask_x], 
                    sigma = sigma[mask_x]
                    )
                
    return(vccr_neg, vccr_pos)

def calculate_non_ccr_arrays(x, win, mid, y1, y2):
    
    vsig = np.nan * np.zeros((win.size,mid.size))
    vmol = np.nan * np.zeros((win.size,mid.size))
    verr = np.nan * np.zeros((win.size,mid.size))
    vsem = np.nan * np.zeros((win.size,mid.size))
    vder = np.nan * np.zeros((win.size,mid.size))
    vsec = np.nan * np.zeros((win.size,mid.size))
    vshp = np.nan * np.zeros((win.size,mid.size))
    vstd = np.nan * np.zeros((win.size,mid.size))
    vdbw = np.nan * np.zeros((win.size,mid.size))
    vncf = np.nan * np.zeros((win.size,mid.size))
    
    data_pack = region_calculations(
                x = x,
                y1 = y1,
                y2 = y2,
                win = win,
                mid = mid
                )
    
    # Load the mean RC signal
    vsig = data_pack["y1_rc_avg"]
    
    # Load the mean molecular attenuated backscattersignal
    vmol = data_pack["y2_rc_avg"]
    
    # Load the standard deviation of the signal in the fit region normalized 
    vstd = data_pack["y1_rc_std"]

    # Load the standard error of the mean of the signal in the fit region normalized 
    verr = data_pack["y1_rc_sem"]
    
    # The normalization coefficient 
    vncf = vmol / vsig
   
    # Load the standard error of the residual in the fit region 
    vsem = data_pack["residual_rc_sem"] / data_pack["y2_rc_avg"]
    
    # Load the absolute relative standar error of the residual derivative
    vder = np.abs(data_pack["slope_err"] / data_pack["slope"])

    # Load the absolute relative standard error of the first derivative difference of the residual between the first and second halves of the fit region
    vsec = np.abs(data_pack["slope_lhalf_err"] - data_pack["slope_rhalf_err"]) / (data_pack["slope_lhalf"] + data_pack["slope_rhalf"])
    
    # Load the Shapiro-Wilk p value in the fitregion
    vshp = data_pack["shapiro_wilk"]
    
    # Load the Durbin-Watson test values for autocorrelation
    vdbw = data_pack["durbin_watson"]

    # Load 
    vext = data_pack["residual_extinction"]

    # Load 
    vexe = data_pack["residual_extinction_error"]
        
    stats = {
        'relative_sem' : vsem,
        'first_derivative' : vder,
        'second_derivative' : vsec,
        'shapiro_wilk' : vshp,
        'durbin_watson' : vdbw,
        'normalization_factor' : vncf,
        'mean_signal' : vsig,
        'mean_molecular_signal' : vmol,
        'signal_std' : vstd,
        'signal_sem' : verr,
        'residual_extinction' : vext,
        'residual_extinction_error' : vexe
        }
    
    for key in stats.keys():
        stats[key] = xr.DataArray(stats[key], dims = ['window', 'middle_point'],
                                  coords = [win, mid])
    
    return(stats)
                
def get_array_coords(x, keyw_args):
  
    min_win = keyw_args['fit_mask_window'][0]
    max_win = keyw_args['fit_mask_window'][1]
    win_step = keyw_args['fit_mask_window_step']
    mask_llim = keyw_args['fit_mask_region'][0]
    mask_ulim = keyw_args['fit_mask_region'][1] 
    
    mid = np.arange(mask_llim, mask_ulim + win_step / 2., win_step)
    
    win = np.arange(min_win, max_win + win_step / 2., win_step)
    
    return(win, mid)
    
def region_calculations(x, y1, y2, win, mid):
    
    data_keys = [
        "y1_rc_avg",
        "y2_rc_avg" ,
        "y1_rc_std",
        "y2_rc_std",
        "y1_rc_sem",
        "y2_rc_sem",
        "norm_factor_rc",
        "residual_rc_sem",
        # "y1_raw_avg",
        # "y2_raw_avg",
        # "y1_raw_std",
        # "y2_raw_std",
        # "y1_raw_sem",
        # "y2_raw_sem",
        # "norm_factor_raw",
        # "residual_raw_sem",
        "slope",
        "slope_err",
        "slope_lhalf",
        "slope_lhalf_err",
        "slope_rhalf",
        "slope_rhalf_err",
        "shapiro_wilk",
        "durbin_watson",
        "residual_extinction",
        "residual_extinction_error"
        ]
    
    
    data_pack = np.nan * np.zeros((len(data_keys), win.size, mid.size))

    for i in range(win.size):
    
        for j in range(mid.size):
      
            mask_x = (x >= mid[j] - win[i] / 2) & (x < mid[j] + win[i] / 2)
            
            mask_y = (y1[mask_x] == y1[mask_x]) & (y2[mask_x] == y2[mask_x])
                        
            if np.sum(mask_y) >= 0.95 * win[i] / (x[1] - x[0]):
                
                # RC section
                x_region = x[mask_x][mask_y]
                y1_rc_region = y1[mask_x][mask_y]
                y2_rc_region = y2[mask_x][mask_y]
                
                y1_rc_avg = np.nanmean(y1_rc_region)
                y2_rc_avg = np.nanmean(y2_rc_region)
            
                y1_rc_std = np.nanstd(y1_rc_region)
                y2_rc_std = np.nanstd(y2_rc_region)
            
                y1_rc_sem = y1_rc_std / np.sqrt(y2_rc_region.size)
                y2_rc_sem = y2_rc_std / np.sqrt(y2_rc_region.size)
            
                norm_factor_rc = y2_rc_avg / y1_rc_avg
                
                residual_rc_region = norm_factor_rc * y1_rc_region - y2_rc_region
                
                residual_rc_sem = np.nanstd(residual_rc_region) / np.sqrt(residual_rc_region.size)
                
                # # Raw section
                # x2_region = np.power(x_region, 2)
                
                # y1_raw_region = y1_rc_region / x2_region
                # y2_raw_region = y2_rc_region / x2_region
                
                # y1_raw_avg = np.nanmean(y1_raw_region)
                # y2_raw_avg = np.nanmean(y2_raw_region)
                
                # y1_raw_std = np.nanstd(y1_raw_region)
                # y2_raw_std = np.nanstd(y2_raw_region)
            
                # y1_raw_sem = y1_raw_std / np.sqrt(y2_raw_region.size)
                # y2_raw_sem = y2_raw_std / np.sqrt(y2_raw_region.size)
                
                # norm_factor_raw = y2_raw_avg / y1_raw_avg
                
                # residual_raw_region = norm_factor_raw * y1_raw_region - y2_raw_region
                
                # residual_raw_sem = np.nanstd(residual_raw_region) / np.sqrt(residual_raw_region.size)

                fit, fit_lh, fit_rh = \
                    derivative_caluclations(
                        x_region = x_region, 
                        residual_region = residual_rc_region
                        )
                    
                shp = shapiro(residual_rc_region)[1]
                
                # if i == 1:
                    # plt.hist(residual_rc_region)
                    # plt.title(f"RC win: {win[i]} km, mid : {np.round(mid[j],1)} km, SPH: {np.round(shp,3)}")
                    # plt.xlim([-5E-6,5E-6])
                    # plt.show()

                
                # Durbin-Watson test for autocorrelation
                dbw = durbin_watson(residual_rc_region)
                    
                # Aerosol extinction error calculation
                ext, ext_err = residual_extinction_caluclations(
                    x_region = x_region, 
                    y1_region = y1_rc_region,
                    y2_region = y2_rc_region,
                    norm_factor = norm_factor_rc
                    )

                data_pack[0,i,j] = y1_rc_avg
                data_pack[1,i,j] = y2_rc_avg 
                data_pack[2,i,j] = y1_rc_std
                data_pack[3,i,j] = y2_rc_std
                data_pack[4,i,j] = y1_rc_sem
                data_pack[5,i,j] = y2_rc_sem
                data_pack[6,i,j] = norm_factor_rc
                data_pack[7,i,j] = residual_rc_sem
                # data_pack[8,i,j] = y1_raw_avg
                # data_pack[9,i,j] = y2_raw_avg
                # data_pack[10,i,j] = y1_raw_std
                # data_pack[11,i,j] = y2_raw_std
                # data_pack[12,i,j] = y1_raw_sem
                # data_pack[13,i,j] = y2_raw_sem
                # data_pack[14,i,j] = norm_factor_raw
                # data_pack[15,i,j] = residual_raw_sem
                data_pack[8,i,j] = fit[4]
                data_pack[9,i,j] = fit[0]
                data_pack[10,i,j] = fit_lh[4]
                data_pack[11,i,j] = fit_lh[0]
                data_pack[12,i,j] = fit_rh[4]
                data_pack[13,i,j] = fit_rh[0]
                data_pack[14,i,j] = shp
                data_pack[15,i,j] = dbw
                data_pack[16,i,j] = ext
                data_pack[17,i,j] = ext_err
                
    data_dict = {}
    
    for i in range(len(data_keys)):
        data_dict[data_keys[i]] = data_pack[i,:,:]
        
    return(data_dict)
    
def derivative_caluclations(x_region, residual_region):
    
    hbin = int(x_region.size / 2.)
    
    fit = linregress(x = x_region, y = residual_region)
   
    fit_lh = linregress(x = x_region[:hbin], y = residual_region[:hbin])
    
    fit_uh = linregress(x = x_region[hbin:], y = residual_region[hbin:]) 

    return(fit, fit_lh, fit_uh)    

    
def residual_extinction_caluclations(x_region, y1_region, y2_region, norm_factor):
    
    # mask_non_zero = (y1_region != 0.) & (y2_region != 0.)
    
    # y1 = y1_region[mask_non_zero]
    # y2 = y2_region[mask_non_zero]
    
    y_ratio = norm_factor * y1_region / y2_region #resverse the ratio to avoid having y1 in the denominator that can be 0 due to errors -- > sign adapted for res_ext from - to +
        
    fit = linregress(x = x_region, y = y_ratio)
    
    # y_ratio_m = norm_factor * np.mean(y1_region) / np.mean(y2_region)
    
    if np.mean(y1_region) > 0.:
        # res_ext = 1e3 * (fit[0] / y_ratio_m) / 2.
        res_ext = 1e3 * fit[0] / 2. 
        res_ext_err = 1e3 * fit[4] / 2.

    else:
        res_ext = np.nan
        res_ext_err = np.nan
        
    
    return(res_ext, res_ext_err)    

def add_ccr_arrays(stats, x, y1, y1_err, y2, ccr_direction, allowed_ccr_directions):
    
    vsig = stats["mean_signal"]
    vmol = stats["mean_molecular_signal"]
    vstd = stats["signal_std"]

    vncf = stats["normalization_factor"]
    
    mid = vsig.middle_point.values
    win = vsig.window.values
    
    if y1_err is None or y1 is None:
        
        if y1 != None or y1_err != None:
            CustomWarning("The y1_err and y1_avg must be either provided together ir not at all. Providing one of them is equivalent to providing none of them")
        
        vccr_neg, vccr_pos = calculate_ccr_arrays(
            x = mid, 
            mid = mid, 
            win = win, 
            c_n = vncf.values, 
            S = vsig[0,:].values, 
            S_m = vmol[0,:].values, 
            sigma = vstd[0,:].values,
            ccr_direction = ccr_direction, 
            allowed_ccr_directions = allowed_ccr_directions
            )

    else:     
        vccr_neg, vccr_pos = calculate_ccr_arrays(
            x = x, 
            mid = mid, 
            win = win, 
            c_n = vncf.values, 
            S = y1, 
            S_m = y2, 
            sigma = y1_err,
            ccr_direction = ccr_direction, 
            allowed_ccr_directions = allowed_ccr_directions
            )

    stats["cross_criterion_neg"] = vccr_neg
    stats["cross_criterion_pos"] = vccr_pos
    
    for key in ["cross_criterion_neg", "cross_criterion_pos"]:
        stats[key] = xr.DataArray(stats[key], dims = ['window', 'middle_point'],
                                  coords = [win, mid])
        
                # if i == 0 and mid[j] > 4. and mid[j] < 10:
                #     plt.plot(x[mask_x],(c * y1[mask_x] - y2[mask_x])/(c * y1_err[mask_x]))
                #     plt.plot(x[mask_x],(y2[mask_x]))
                #     plt.title(f"win: {win[i]}, mid: {mid[j]}")
                #     plt.show()

    return(stats)

def create_masks(stats, keyw_args, cancel_stats, ccr_type,):
    
    # Unpack keyword arguments
    sem_threshold = keyw_args['rsem_threshold']
    der_threshold = keyw_args['first_derivative_threshold']
    sec_threshold = keyw_args['second_derivative_threshold'] 
    shp_threshold = keyw_args['shapiro_wilk_threshold']
    ccr_threshold = keyw_args['cross_criterion_threshold']
    dbw_threshold = keyw_args['durbin_watson_threshold']
    ext_threshold = keyw_args['residual_extinction_threshold']

    win = stats['mean_signal'].window.values
    mid = stats['mean_signal'].middle_point.values
    
    vsem = stats['relative_sem'].values
    vder = stats['first_derivative'].values
    vsec = stats['second_derivative'].values
    vshp = stats['shapiro_wilk'].values
    vccr_neg = stats['cross_criterion_neg'].values
    vccr_pos = stats['cross_criterion_pos'].values
    vdbw = stats['durbin_watson'].values
    vsig = stats['mean_signal'].values
    vext = stats['residual_extinction'].values
    vexe = stats['residual_extinction_error'].values

    # Calculate the standard error mask
    if 'sem' in cancel_stats: 
        msem = np.ones_like(vsem, dtype = bool)
    else:
        msem = (vsem <= sem_threshold)
        msem[msem != msem] = False
    
    # Calculate the derivative mask (pval < 0.05 means the slope is significant)
    if 'der' in cancel_stats: 
        mder = np.ones_like(vder, dtype = bool)
    else:
        mder = vder <= der_threshold
        mder[mder != mder] = False
    
    # Check if the derivatives inside the 2 halves of the window are aslo not significant
    if 'sec' in cancel_stats: 
        msec = np.ones_like(vsec, dtype = bool)
    else:
        msec = vsec <= sec_threshold
        msec[msec != msec] = False
    
    # Check if the p value of the Saphiro-Wilkinson test is higher than 0.95 (<0.05 means not normal)
    if 'shp' in cancel_stats: 
        mshp = np.ones_like(vshp, dtype = bool)
    else:
        mshp = vshp > shp_threshold
        mshp[mshp != mshp] = False
    
    # Cross check test calculations
    if 'ccr' in cancel_stats: 
        mccr = np.ones_like(vccr_pos, dtype = bool)
    else:
        if ccr_type == 'above_negative_threshold':
            # mccr = vccr_min >= -ccr_threshold
            mccr = vccr_neg <= ccr_threshold
        elif ccr_type == 'below_positive_threshold':
            # mccr = vccr_max <= ccr_threshold
            mccr = vccr_pos <= ccr_threshold
        elif ccr_type == 'between_both_thresholds':
            mccr = (vccr_neg <= ccr_threshold) & (vccr_pos <= ccr_threshold)
        else:
            raise Exception(f"-- Error: The provided ccr_type {ccr_type} is wrong. Please select one of: {allowed_ccr_types}")    
        mccr[mccr != mccr] = False
    
    # Check for autocorrelation
    if 'dbw' in cancel_stats: 
        mdbw = np.ones_like(vdbw, dtype = bool)
    else:
        mdbw = (vdbw >= dbw_threshold[0]) & (vdbw <= dbw_threshold[1])
        mdbw[mdbw != mdbw] = False
    
    # Check for negative signal values
    if 'pos' in cancel_stats: 
        mpos = np.ones_like(vsig, dtype = bool)
    else:
        mpos = vsig > 0.
        mpos[mpos != mpos] = False
    
    if 'ext' in cancel_stats: 
        mext = np.ones_like(vext, dtype = bool)
    else:
        # mext = (vext <= 0.) & (vext >= -ext_threshold)
        # mext = (vext - vexe <= 0.) & (vext + vexe >= -ext_threshold)
        mext = np.abs(vext) <= ext_threshold
        mpos[mext != mext] = False
    
    # mpsn = (vpsn > 0.9) & (vpsn > 1.1)
    # mpsn[mpsn != mpsn] = False
    # if 'psn' in cancel_stats: mpsn = np.ones_like(mpsn, dtype = bool)
    
    mfit = (msem) & (mder) & (msec) & (mshp) & (mccr) & (mpos) & (mdbw) & (mext)

    masks = {
        'relative_sem' : msem,
        'first_derivative' : mder,
        'second_derivative' : msec,
        'shapiro_wilk' : mshp,
        'cross_criterion' : mccr,
        'durbin_watson' : mdbw,
        'is_positive' : mpos,
        'residual_extinction' : mext,
        'total' : mfit
        }
   
    for key in masks.keys():
        masks[key] = xr.DataArray(masks[key], dims = ['window', 'middle_point'],
                                  coords = [win, mid])
        
    return(masks)

def statistics(y1, y2, x, y1_err = None, y1_avg = None, keyw_args = {}, 
               ccr_direction = 'apply_backwards', 
               ccr_type = 'above_negative_threshold', 
               cancel_stats = []):
    
    """
    keyw_args: A dictionary of arguments that can be passed to configure the mask. These are namely:
        'fit_mask_window' 
        'fit_mask_window_step'
        'fit_mask_region'
        'rsem_threshold'
        'first_derivative_threshold'
        'second_derivative_threshold'
        'shapiro_wilk_threshold'
        'cross_criterion_threshold'
        'durbin_watson_threshold'
        'is_positive'
        'residual_extinction'
        'check_poisson'
    
    cancel_stats: A list with the names of the statistical tests not to be considered. These are:
        'sem': Relative standard error of the mean
        'der': First derivative
        'sec': Second derivative (first derivative differences)
        'shp': Shapiro Wilk test
        'ccr': Cross criterion
        'dbw': Durbin Watson
        'pos': Positive signal
        'ext': Residual extinction
        'psn': Poisson distribution
        For example: cancel_stats = ['sem', 'ccr'] mean that the relative standard error of the mean and the cross criterion tests will not be considered in the mask
    
    ccr_direction: Configure where the cross criterion will be applied. Select one of the following sting values:
        'apply_backwards': Scan only regions that are below the normalization region (down to the beginning of the fit_mask_region)
        'apply_forwards': Scan only regions that are above the normalization region (uo to the end of the fit_mask_region)
        'apply_everywhere': Scan the whole profile (within the fit_mask region)
    
    cross_ccr_type: Configure how the cross criterion will be applied. Select one of the following sting values:
        'above_negative_threshold': The minimum difference between the testing signal y1 and the reference signal y2 within the region where the cross criterion is applied must be smaller than the negative uncertainty of y1 multiplied by the cross_criterion_threshold 
        'below_positive_threshold': The maximum difference between the testing signal y1 and the reference signal y2 within the region where the cross criterion is applied must be larger than the positive uncertainty of y1 multiplied by the cross_criterion_threshold
        'between_both_thresholds': Both of the aforementioned thresholds will be applied. The absolute difference between the testing signal y1 and the reference signal y2 within the region where the cross criterion is applied must be smaller than the (positive) uncertainty of y1 multiplied by the cross_criterion_threshold
    """
    
    # Fill in keyword arguments with default values if missing
    for key in needed_keys:
        if key not in keyw_args.keys():
            keyw_args[key] = default_values[key]
            
    # Get mask coords            
    win, mid = get_array_coords(x, keyw_args)

    # Get smoothed signal and uncertainty (150 m window) to use for the ccr test
    y1_avg, y1_err = \
        smooth_1D(
            y_vals = y1, 
            x_vals = x,
            x_sm_lims = keyw_args['fit_mask_region'],
            x_sm_win = 150.,
            expo = False,
            err_type = 'std'
            )
  
    # Calculate values used for masking
    stats = calculate_non_ccr_arrays(
        x = x, 
        win = win, 
        mid = mid, 
        y1 = y1, 
        y2 = y2
        )

    # Calculate valuesused for cross_crit masking
    stats = add_ccr_arrays(
        stats = stats, 
        x = x, 
        y1 = y1_avg, 
        y1_err = y1_err, 
        y2 = y2, 
        ccr_direction = ccr_direction, 
        allowed_ccr_directions = allowed_ccr_directions
        )


    masks = create_masks(
        stats = stats, 
        keyw_args = keyw_args, 
        cancel_stats = cancel_stats, 
        ccr_type = ccr_type
        )
    

    # from matplotlib import pyplot as plt
    # (vccr_neg).plot()
    # plt.show()
    # (vccr_pos).plot()
    # plt.show()
    # mccr.plot()
    # plt.show()
    # vsig.plot()
    # plt.show()
    # verr.plot()
    # plt.show()
    # vmol.plot()
    # plt.show()
    # vncf.plot()
    # plt.show()
    # (stats['residual_extinction']).plot()
    # plt.show()
    
    
    return(stats, masks)

       
def scan(masks, prefered_range = 'far', default_norm_region = None, auto_fit = True):
        
    mfit = masks['total']

    window = mfit.window

    middle_point = mfit.middle_point    
    
    window_matrix = window * xr.ones_like(mfit)
    
    window_matrix = window_matrix.where(mfit,np.nan)

    middle_point_matrix = xr.ones_like(mfit) * middle_point
    
    middle_point_matrix = middle_point_matrix.where(mfit,np.nan)

    if mfit.any() and auto_fit == True:

        norm_region_flag = "auto"
        
        if prefered_range == 'far':
            
            middle_point_max = middle_point_matrix.max().values
            middle_point_max_idx = np.where(middle_point == middle_point_max)[0][0]
            
            window_max = window_matrix.isel({'middle_point' : middle_point_max_idx}).max().values
            window_max_idx = np.where(window == window_max)[0][0]
            
            auto_norm_region = [middle_point_max - window_max / 2., 
                                middle_point_max + window_max / 2.]
            
            auto_idx = [window_max_idx, middle_point_max_idx]


        elif prefered_range == 'near':
            
            middle_point_min = middle_point_matrix.min().values
            middle_point_min_idx = np.where(middle_point == middle_point_min)[0][0]
            
            window_max = window_matrix.isel({'middle_point' : middle_point_min_idx}).max().values
            window_max_idx = np.where(window == window_max)[0][0]
            
            auto_norm_region = [middle_point_min - window_max / 2., 
                                middle_point_min + window_max / 2.]
            
            auto_idx = [window_max_idx, middle_point_min_idx]
            
        else:
            raise Exception(f"-- Error: The provided prefered_range {prefered_range} is wrong. Select one of: near, far")
    
    else:        
        norm_region_flag = "external"

        auto_norm_region = default_norm_region
        
        auto_idx = idx_from_region(default_norm_region, middle_point, window)
            
    return(auto_norm_region, norm_region_flag, auto_idx)   
    