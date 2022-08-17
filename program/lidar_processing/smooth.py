"""
@authors: N. Siomos & P. Paschou

Smoothing of signals (optional derivation) by applying polynomial fit  
"""
import sys
import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval, polyder
import xarray as xr
import matplotlib.pyplot as plt

def get_win(ranges, flt_node_ranges, flt_node_windows, resol):
    
    '''
    General:
        Calculates the half smoothing windows 
        
    Input:
        ranges :
            A 2D xarray with dimensions (time, channel) that includes 
            the singal range values per channel and bin
            
        flt_node_ranges : numpy array
            A numpy array with the range limits of each smoothing node
            
        flt_node_windows : 
            A numpy array with the window length of each smoothing node
            
        resol:
            A pandas series with the range resolution in m per channel. 
            The index should correspond to the channel dimension of sig
                

    Output
        half_flt_window : (M,) array
            An xarray with the half window in bins for each range and channel

    '''

    channels = ranges.channel.values
    
    # Use the uppermost signal range in the last smoothing node if it is higher 
    if flt_node_ranges and sig.range.values[-1] > lims[-1]:
        lims[-1] = sig.range.values[-1]

    half_flt_window = np.nan * ranges.copy()
    
        
    for ch in channels:

        ch_d = dict(channel = ch)
    
    
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
    half_flt_window = np.floor(win/(2.*step)).astype(int)

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