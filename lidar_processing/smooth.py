"""
@authors: N. Siomos & P. Paschou

Smoothing of signals (optional derivation) by applying polynomial fit  
"""
import sys
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval, polyder
import xarray as xr
import matplotlib.pyplot as plt

def get_win(x, lims, wins, step):
    '''

    Parameters
    ----------
    x : (M,) array
        singal range values.
    lims : array
        range limits of each signal node.
    wins : array
        window length limits of each signal node.
    step : float
        resolution.

    Returns
    -------
    ihwin : (M,) array
        The half window bins for each range .

    '''
    win = np.zeros(x.shape)
    
    # set a variable window value per range bin (0 below first range limit)
    mask_win = (x > lims[0]) & (x <= lims[-1])
    win[mask_win] = np.interp(x[mask_win], lims, wins)
        
    # Limit half window values to avoid setting
    # a window broader than the signal limits
    mask_down = (x - win/2. < step)
    mask_up = (x[-1] - x - win/2. < step)
    win[mask_down] = 2.*(x[mask_down] - step)
    win[mask_up] = 2.*(x[-1] - x[mask_up])

    # Convert float variable window to even half bin window
    ihwin = np.floor(win/(2.*step)).astype(int)

    return(ihwin)

def smoothing(signal, ihwin, order, deriv_order=0, itype='sig'):
    
    iters = signal.iters.values
    time = signal.time.values
    if itype == 'sig':
        channel_id = signal.channel.values
        channel_nm = 'channel'
    if itype == 'prod':
        channel_id = signal.product.values
        channel_nm = 'product'
        
    x = signal.range.values
    y = signal.values
    
    y_sm = np.copy(y)
    
    ind = np.arange(0, x.size, 1)
    slices = [slice(i - ihwin[i], i + ihwin[i]+1) for i in ind]
    
    # window points less than polynomial order
    mask_ord = (2*ihwin + 1 >= order)
    
    
    for n in range(iters.size):
        # n = 0 # Apply the smoothing only in the measured signal.
        for k in range(time.size):    
            for i in ind[mask_ord == True]:
                arg = y[n, k, :, slices[i]]
                # Mask to exclude signal slices with spare NaNs
                mask_ch = (arg == arg).sum(axis=1)
                # ch_ind = np.where(mask_ch != 0)[0] # exclude channels slices full of NaNs
                ch_ind = np.where(mask_ch == arg.shape[1])[0]
                
                p = polyfit(x[slices[i]], arg[ch_ind, :].T, order)
                
                if deriv_order >0: # calculate or not the derivative
                    p = polyder(p, deriv_order)
                
                y_sm[n, k, ch_ind, i] = polyval(x[i], p) 
    
    signal = xr.DataArray(y_sm, 
                          coords=[iters, time, channel_id, x],
                          dims=['iters', 'time', channel_nm, 'range'])
    
    return(signal)

def smoothing_1d(sig, sig_range, ihwin, order, deriv_order):
    
    sig_sm = np.copy(sig)
    
    ind = np.arange(0, sig_range.size, 1)
    slices = [slice(i - ihwin[i], i + ihwin[i]+1) for i in ind]
    
    # window points less than polynomial order and signal values not NaN
    mask_ord = ((2*ihwin + 1) >= order) & (sig == sig)    
    
    for i in ind[mask_ord == True]:
        
        # mask the sig val for spare NaNs
        mask_nan = sig[slices[i]] == sig[slices[i]]
        x = sig_range[slices[i]][mask_nan==True]
        y = sig[slices[i]][mask_nan==True]

        p = polyfit(x, y.T, order)
        #p = polyfit(sig_range[slices[i]], sig[slices[i]].T, order)

        if deriv_order > 0: # calculate or not the derivative
            p = polyder(p, deriv_order)
        
        sig_sm[i] = polyval(sig_range[i], p) 
    
    return(sig_sm)