#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 17:06:17 2022

@author: nick
"""

import numpy as np

import xarray as xr

def detect_overflows(sig, dead_time, daq_range):

    """
    General:
        Detects lidar profiles with photon values close to the
        maximum allowed countrate or analog values close to the data 
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
            
    """
    
    channels = sig.channel.values
    
    for ch in channels:

        ch_d = dict(channel = ch)
        
        if ch[6] == 'p': #7th digit of channel name is the acquisition mode (a or p)
        
            max_countrate = 0.95 * 1000./dead_time[ch]

            if (sig.loc[ch_d].values >= max_countrate).any():
                print(f"-- Warning: Channel {ch} - Photon signal countrate values above the 95% of the maximum allowed value were detected! ")


        if ch[6] == 'a': #7th digit of channel name is the acquisition mode (a or p)
 
            max_mV = 0.95*daq_range[ch]

            if (sig.loc[ch_d].values >= max_mV).any():
                print(f"-- Warning: Channel {ch} - Analog signal mV values above the 95% of the data acqusition range were detected! ")   

    return()