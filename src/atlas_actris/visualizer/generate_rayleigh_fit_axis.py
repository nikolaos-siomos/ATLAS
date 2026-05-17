#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:56:30 2026

@author: nikos
"""


import numpy as np

def get_rayleigh_fit_y_limits(y1_vals, y2_vals, y_lims, wavelength, use_lin_scale):
    """
    Determine y-axis limits.

    Parameters
    ----------
    y1_vals : array-like
        Values used to estimate the upper y limit.
    y2_vals : array-like
        Values used to estimate the lower y limit.
    y_lims : list
        Either [] for automatic limits, or [lower, upper] for manual limits.
    wavelength : float
        Signal wavelength.
    use_lin_scale : bool
        If False, sanity checks are applied for logarithmic scale.

    Returns
    -------
    list
        [lower_limit, upper_limit]
    """

    scale_f = wavelength / 355.0
    scat_ratio_f = 2.5

    # Automatic limits
    if len(y_lims) == 0:
        auto = True
        y_max = np.nanmax(y1_vals)
        y_min = np.nanmin(y2_vals)

        y_llim = y_min / 2.0
        y_ulim = scat_ratio_f * scale_f * y_max

    # Manual limits
    else:
        auto = False
        y_llim = y_lims[0]
        y_ulim = y_lims[-1]

    # Sanity checks for logarithmic scale
    if not use_lin_scale:

        if y_llim <= 0:
            if not auto:
                print(
                    "-- Warning: y-axis lower limit <= 0 although the scale is "
                    "logarithmic. The limit has automatically been replaced."
                )
            y_llim = 1e-6

        if y_ulim <= 0:
            if not auto:
                print(
                    "-- Warning: y-axis upper limit <= 0 although the scale is "
                    "logarithmic. The limit has automatically been replaced."
                )
            y_ulim = 1.0

    return [y_llim, y_ulim]