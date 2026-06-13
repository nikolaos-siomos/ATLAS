#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:56:30 2026

@author: nikos
"""


import numpy as np


def _finite_values(vals):
    vals = np.asarray(vals, dtype=float)
    return vals[np.isfinite(vals)]


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

    # Manual limits
    if len(y_lims) != 0:
        auto = False
        y_llim = y_lims[0]
        y_ulim = y_lims[-1]

        if not np.isfinite(y_llim) or not np.isfinite(y_ulim):
            print(
                "-- Warning: provided y_lims contain NaN/Inf. "
                "Falling back to automatic y-axis limits."
            )
            y_lims = []

    # Automatic limits
    if len(y_lims) == 0:
        auto = True

        y1_finite = _finite_values(y1_vals)
        y2_finite = _finite_values(y2_vals)

        if y1_finite.size == 0 or y2_finite.size == 0:
            print(
                "-- Warning: cannot determine Rayleigh-fit y-axis limits because "
                "Y1 or Y2 has no finite values. Using fallback limits."
            )

            if use_lin_scale:
                return [0.0, 1.0]
            else:
                return [1e-6, 1.0]

        y_max = np.nanmax(y1_finite)
        y_min = np.nanmin(y2_finite)

        y_llim = y_min / 2.0
        y_ulim = scat_ratio_f * scale_f * y_max

    # Final finite check
    if not np.isfinite(y_llim) or not np.isfinite(y_ulim):
        print(
            "-- Warning: calculated Rayleigh-fit y-axis limits are NaN/Inf. "
            "Using fallback limits."
        )

        if use_lin_scale:
            return [0.0, 1.0]
        else:
            return [1e-6, 1.0]

    # Avoid equal or reversed limits
    if y_ulim <= y_llim:
        print(
            "-- Warning: calculated Rayleigh-fit y-axis limits are invalid "
            f"({y_llim}, {y_ulim}). Using expanded fallback limits."
        )

        if use_lin_scale:
            center = y_llim
            spread = max(abs(center) * 0.1, 1.0)
            return [center - spread, center + spread]
        else:
            return [1e-6, max(1.0, y_ulim)]

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
