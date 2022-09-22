#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:27:01 2022

@author: nick
"""

import numpy as np

def check_channels(sel_channels, all_channels):
    
    if not isinstance(sel_channels,type(None)):
        sel_channels = np.array(sel_channels)
        all_channels = np.array(all_channels)

        missing_ch = [ch not in all_channels for ch in sel_channels]
        
        if any(missing_ch):
            raise Exception("-- Error: The following provided channels do not exist: "+\
                            f"{sel_channels[missing_ch]} \n Please select one of:"+\
                            f"{all_channels}")
                
        channels = sel_channels
    
    else:
        channels = all_channels

    return(channels)