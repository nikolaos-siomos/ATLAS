#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:37:48 2023

@author: nick
"""

import numpy as np
from . import normalize
from ..tools.smoothing import sliding_average_1D_fast as smooth_1D_fast
from ..tools.smoothing import sliding_average_2D_fast as smooth_2D_fast            

def process(x, y, iters, smooth, x_sm_lims, x_sm_win, expo, region):
            
    # Averaging sector signals
    y_m = np.nanmean(y[:iters,:], axis = 0) 

    # Smoothing
    if smooth:
        if not isinstance(x_sm_win,list):
            from ..tools.smoothing import sliding_average_1D_fast as smooth_1D
            from ..tools.smoothing import sliding_average_2D_fast as smooth_2D
        else:
            from ..tools.smoothing import sliding_average_1D as smooth_1D
            from ..tools.smoothing import sliding_average_2D as smooth_2D
    
        # Smoothing averaged sectors
        y_m_sm, _ = \
            smooth_1D(y_vals = y_m, 
                      x_vals = x,
                      x_sm_lims = x_sm_lims,
                      x_sm_win = x_sm_win,
                      expo = expo)

        
        y_m_sm, _ = \
            smooth_1D_fast(y_vals = y_m_sm, 
                           x_vals = x,
                           x_sm_lims = [x_sm_lims[1] - 0.3,20. + 0.3],
                           x_sm_win = 500.,
                           expo = expo)

        # Smoothing unaveraged sectors
        y_sm = \
            smooth_2D(z_vals = y, 
                      y_vals = x,
                      y_sm_lims = x_sm_lims,
                      y_sm_win = x_sm_win,
                      expo = expo)

        y_sm = \
            smooth_2D_fast(z_vals = y_sm, 
                           y_vals = x,
                           y_sm_lims = [x_sm_lims[1] - 0.3,20. + 0.3],
                           y_sm_win = 500.,
                           expo = expo)
        
    else:
        y_m_sm = y_m
        
        y_sm = y

    # Store the extra iteration if it exists 
    if y_sm.shape[0] == iters + 1:
        y_extra = y[iters,:]
        y_extra_sm = y_sm[iters,:]
        extra = True
    else:
        y_extra = []
        y_extra_sm = []
        extra = False
        
    # Minimum value per bin of all iterations
    y_l_sm = np.nanmin(y_sm[:iters,:], axis = 0)

    # Maximum value per bin of all iterations
    y_u_sm = np.nanmax(y_sm[:iters,:], axis = 0)

    # Normalization for smoothed signals
    coef, norm_bin = normalize.to_a_point(sig = y_m, 
                                          sig_b = np.ones(x.shape), 
                                          x_vals = x,
                                          region = region,
                                          axis = 0)
    
    return(coef, y_m, y_sm, y_m_sm, y_l_sm, y_u_sm, y_extra, y_extra_sm, extra)