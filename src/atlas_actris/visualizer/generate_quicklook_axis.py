#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 13:10:34 2026

@author: nikos
"""

import numpy as np
from utils.toolbox import round_it

def get_quicklook_x_limits(x_vals, x_lims):

    # Get the altitude/distance lower limit and bin
    x_lbin = np.where(x_vals >= x_lims[0])[0][0]
    
    if x_lbin > 0:
        x_lbin = x_lbin - 1
    
    x_llim = x_lims[0]

    # Get the altitude/distance upper limit and bin
    x_ubin = np.where(x_vals <= x_lims[-1])[0][-1] 
    
    if x_ubin < x_vals.size:
        x_ubin = x_ubin + 1

    x_ulim = x_lims[-1]

    return [x_llim, x_ulim]

def get_quicklook_y_limits(y_vals, x_vals, y_lims, use_log):
    """
    Determine quicklook y-axis limits.

    Parameters
    ----------
    y_vals : array-like
        Signal values, preferably already normalized.
        Shape usually (time, bins).
    x_vals : array-like
        Range/height values with shape usually (bins,).
    y_lims : list
        [] for automatic limits, [lower] for fixed lower and automatic upper,
        or [lower, upper] for fixed limits.
    use_log : bool
        Whether the y-axis is logarithmic.

    Returns
    -------
    y_llim : float
        Lower y-axis limit.
    y_ulim : float
        Upper y-axis limit.
    """

    # Treat None or [] as automatic limits
    if y_lims is None or len(y_lims) == 0:
        y_lims = [None, None]

    elif len(y_lims) == 1:
        y_lims = [y_lims[0], None]

    elif len(y_lims) != 2:
        raise ValueError("y_lims must be [], [lower], or [lower, upper].")

    y_llim, y_ulim = y_lims

    # Automatic lower limit
    if y_llim is None:

        if use_log:
            y_llim = 1e-2

        else:
            y_llim = 0.0

    # Automatic upper limit
    if y_ulim is None:
        y_ulim = 1.

    # Final sanity checks for log scale
    if use_log:
        if y_llim <= 0 or np.isnan(y_llim):
            y_llim = 1e-2

        if y_ulim <= 0 or np.isnan(y_ulim):
            y_ulim = 1.0

    return [y_llim, y_ulim]

def normalize_quicklook_y(y_vals, x_vals, y_max_zone):
    """
    Normalize quicklook signal using the maximum mean signal
    inside a selected x-zone.

    Parameters
    ----------
    y_vals : array-like
        Signal values with shape usually (time, bins).
    x_vals : array-like
        Range/height values with shape usually (bins,).
    y_max_zone : list or tuple
        Zone [x_min, x_max] used to estimate the normalization maximum.

    Returns
    -------
    y_vals_nrm : array-like
        Normalized signal.
    y_max : float
        Normalization factor.
    """

    mask_max_zone = (x_vals >= y_max_zone[0]) & (x_vals <= y_max_zone[1])

    y_vals_sm = np.nanmean(y_vals[:, mask_max_zone], axis=0)

    y_max = round_it(np.nanmax(y_vals_sm), 1)

    # Avoid division by zero or invalid normalization
    if y_max == 0 or np.isnan(y_max):
        y_max = 1.0

    y_vals_nrm = y_vals / y_max

    return y_vals_nrm
