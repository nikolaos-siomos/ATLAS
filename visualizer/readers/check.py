#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:27:01 2022

@author: nick
"""

import numpy as np

def check_channels(sel_channels, all_channels, exclude_telescope_type,
                   exclude_channel_type, exclude_acquisition_mode,
                   exclude_channel_subtype):
    
    if not isinstance(sel_channels,type(None)):
        if not isinstance(sel_channels,list):
            sel_channels = np.array([sel_channels])
        else:
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
    
    if not isinstance(exclude_telescope_type,list):
        exclude_telescope_type = [exclude_telescope_type]

    if not isinstance(exclude_channel_type,list):
        exclude_channel_type = [exclude_channel_type]
        
    if not isinstance(exclude_acquisition_mode,list):
        exclude_acquisition_mode = [exclude_acquisition_mode]

    if not isinstance(exclude_channel_subtype,list):
        exclude_channel_subtype = [exclude_channel_subtype]
            
    # Make a channel maskand exclude channels that don't fit the criteria
    mask = np.array([ch[4] not in exclude_telescope_type and
                     ch[5] not  in exclude_channel_type and
                     ch[6] not  in exclude_acquisition_mode and
                     ch[7] not  in exclude_channel_subtype 
                     for ch in channels])

    if all(~mask):
        raise Exception('-- Error: The provided channel filtering arguments are too strict and exclude all channels. Please revide the following arguments: exclude_telescope_type, exclude_channel_type, exclude_acquisition_mode, exclude_channel_subtype, channels')
    
    channels = channels[mask]

    return(channels)

def check_channels_no_exclude(sel_channels, all_channels):
    
    if not isinstance(sel_channels,type(None)):
        if not isinstance(sel_channels,list):
            sel_channels = np.array([sel_channels])
        else:
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