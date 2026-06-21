#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 13:34:47 2022

@author: nick
"""

import sys
import numpy as np


def _as_1d_float_array(x_vals):
    """Return x_vals as a 1D float numpy array without changing length/order."""
    x_vals = np.asarray(x_vals, dtype=float)

    if x_vals.ndim != 1:
        raise ValueError(
            "x_vals must be a 1D coordinate array. "
            f"Got shape {x_vals.shape}."
        )

    return x_vals


def _finite_coordinate_limits(x_vals):
    """Return first and last finite coordinates for range checks."""
    x_vals = _as_1d_float_array(x_vals)
    finite = np.isfinite(x_vals)

    if not np.any(finite):
        raise ValueError("x_vals contains no finite coordinates.")

    x_finite = x_vals[finite]

    return x_finite[0], x_finite[-1]


def _coord_step(x_vals):
    """
    Estimate the coordinate step from finite adjacent coordinate pairs.

    This avoids failures when the vertical scale starts or ends with NaNs, while
    preserving the original bin indexing used by downstream plotting/export.
    """
    x_vals = _as_1d_float_array(x_vals)

    dx_vals = np.diff(x_vals)
    dx_vals = dx_vals[np.isfinite(dx_vals) & (dx_vals > 0.0)]

    if dx_vals.size == 0:
        finite = np.isfinite(x_vals)
        raise ValueError(
            "Cannot calculate coordinate step because x_vals has no positive "
            "finite adjacent differences. "
            f"finite_x={np.count_nonzero(finite)}/{x_vals.size}"
        )

    dx = np.nanmedian(dx_vals)

    if not np.isfinite(dx) or dx <= 0.0:
        raise ValueError(f"Invalid coordinate step calculated from x_vals: {dx}")

    return dx


def to_a_point(sig, sig_b, x_vals, region, axis, axis_b = 0):

    norm_height = (region[1] + region[0]) / 2.
    
    hwin = (region[1] - region[0]) / 2.
    
    # Get the reference height bin    
    norm_bin = get_norm_bin(x_vals = x_vals, norm_height = norm_height)

    # Get the reference window in bins    
    hwin_bin = get_hwin_bin(x_vals = x_vals, hwin = hwin)

    sig_sel = choose_from_axis(sig, axis, 
                               norm_bin - hwin_bin, 
                               norm_bin + hwin_bin + 1)

    sig_b_sel = choose_from_axis(sig_b, axis_b, 
                                 norm_bin - hwin_bin, 
                                 norm_bin + hwin_bin + 1)    

    sig_m = np.nanmean(sig_sel, axis = axis, keepdims=True)
    
    sig_b_m = np.nanmean(sig_b_sel, axis = axis_b, keepdims=True)
    
    norm_coef = sig_b_m / sig_m 

    # sig_n = sig * norm_coef

    return(norm_coef, norm_bin)


def get_norm_bin(x_vals, norm_height):

    x_vals = _as_1d_float_array(x_vals)
    x_first, x_last = _finite_coordinate_limits(x_vals)

    if norm_height < x_first:
        raise Exception(f'-- Error: Normalization height/distance  is too low ({norm_height}km) ' +
                        f'while the signal starts at {x_first}km')
        
    elif norm_height > x_last:
        raise Exception(f'-- Error: Reference height/distance is too high ({norm_height}km) ' +
                        f'while the signal ends at {x_last}km')
    else:
        mask = np.isfinite(x_vals) & (x_vals >= norm_height)

        if not np.any(mask):
            raise Exception(
                f'-- Error: Could not find a finite normalization bin at or above {norm_height}km'
            )

        norm_bin = np.where(mask)[0][0]
        
    return(norm_bin)


def get_hwin_bin(x_vals, hwin):

    dx = _coord_step(x_vals)

    if hwin < dx:
        raise Exception(f'-- Error: The half reference window provided ({hwin}km) is ' +
                        'smaller than the signal vertical step')
        
    else:
        hwin_bin = int(hwin / dx)
        
    return(hwin_bin)


def choose_from_axis(a, axis, start, stop):

    if axis <= a.ndim:
    
        s = [slice(None) for i in range(a.ndim)]
        
        s[axis] = slice(start, stop)
        
        s = tuple(s)

    else:
        raise Exception('-- Error: The provided axis index is larger than the number ' +
                        'of the axises of the array')

    return a[s]


def add_axis(a, axis):
    
    if axis <= a.ndm + 1:
        s = [slice(None) for i in range(a.ndim + 1)]
        
        s[axis] = np.newaxis
    else:
        raise('-- Error: The provided axis index is larger than the number ' +
              'of the axises of the array plus one')
    
    return a[s]
