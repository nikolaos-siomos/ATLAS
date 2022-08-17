#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 17:06:17 2022

@author: nick
"""

import numpy as np

import xarray as xr

def detect_overflows(sig, dead_time, daq_range, iscr = False):

    """
    General:
        Detects and removes lidar profiles with photon values above the
        maximum allowed countrate or analog values above the data 
        acquisition range
        
    Input:
        sig: 
            A 2D or 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, ...). 
            
        dead_time: 
            A pandas series with the dead time per channel in nanoseconds. 
            The index should correspond to the channel dimension of sig            

        DAQ_range: 
            A pandas series with the data acquisition range of the analog
            channels.
            The index should correspond to the channel dimension of sig   
        
        iscr:
            A boolean scalar. If set to False, saturated or cliped values are
            just detected (not removed). By stting it to True the problematic
            timeframes will be removed
            
    Returns:
        
        sig_out: 
            An xarray in the same shape as sig with the dead time corrected
            signals
            
    """
    
    channels = sig.channel.values
    
    time = sig.time.values

    mask = xr.DataArray(np.zeros([time.size, channels.size], dtype = bool),
                        dims = ['time','channel'],
                        coords = [time, channels])
    
    sig_out = np.nan * sig.copy()
    
    for ch in channels:

        ch_d = dict(channel = ch)
        
        if ch[2] == 'p': #3rd digit of channel name is the acquisition mode (a or p)
        
            max_countrate = 1000./dead_time[ch]

            if (sig.loc[ch_d].values >= max_countrate).any():
                print(f"-- Warning: Channel {ch} - Photon signal countrate values above the maximum allowed value were detected! ")

                mask.loc[ch_d] = \
                    (sig.loc[ch_d] >= max_countrate).any(dim = 'bins').values
                                            

        if ch[2] == 'a': #3rd digit of channel name is the acquisition mode (a or p)
        
            if (sig.loc[ch_d].values >= 0.95*daq_range[ch]).any():
                print(f"-- Warning: Channel {ch} - Analog signal mV values too close or above the data acqusition range were detected! ")

                mask.loc[ch_d] = \
                    (sig.loc[ch_d] >= 0.95 * daq_range[ch])\
                        .any(dim = 'bins').values

    if iscr: # Remove the problematic profiles
    
        mask_t = mask.any(dim = 'channel').values
        
        print(f"-- Warning: {np.sum(mask_t)} profiles with overflows have been removed ")

        time_d = dict(time = time[~mask_t])
        
        sig_out = sig.loc[time_d]    

    return(sig_out)