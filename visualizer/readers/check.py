#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:27:01 2022

@author: nick
"""

import numpy as np

def check_channels(sel_channels, all_channels, exclude_field_type,
                   exclude_scattering_type, exclude_detection_mode,
                   exclude_channel_subtype):
    
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
    
    # Make a channel maskand exclude channels that don't fit the criteria
    mask = np.array([ch[0] not in exclude_field_type and
                     ch[1] not  in exclude_scattering_type and
                     ch[2] not  in exclude_detection_mode and
                     ch[3] not  in exclude_channel_subtype 
                     for ch in channels])
    
    if all(~mask):
        raise Exception('-- Error: The provided channel filtering arguments are too strict and exclude all channels. Please revide the following arguments: exclude_field_type, exclude_scattering_type, exclude_detection_mode, exclude_channel_subtype, channels')
    
    channels = channels[mask]

    return(channels)

def check_channels_no_exclude(sel_channels, all_channels):
    
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