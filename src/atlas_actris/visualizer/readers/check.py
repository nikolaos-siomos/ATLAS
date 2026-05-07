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

def find_rt_channels(ch_r, ch_t, channels):
    
    if ch_r == None or ch_t == None:
        channels_r = []
        channels_t = []
        ch_r_all = np.array([ch for ch in channels if ch[7] == 'r'])
        ch_t_all = np.array([ch for ch in channels if ch[7] == 't'])
        if len(ch_r_all) == 0:
            print("-- Warning: No relfected channels were detected. Please make sure that the channel_subtype in set correctly in the configuration file")
        if len(ch_t_all) == 0:
            print("-- Warning: No transmitted channels were detected. Please make sure that the channel_subtype in set correctly in the configuration file")
        
        for ch_r_i in ch_r_all:
            for ch_t_i in ch_t_all:
                if ch_r_i[4]  == ch_t_i[4] and ch_r_i[6]  == ch_t_i[6] and \
                    ch_r_i[:4]  == ch_r_i[:4]:
                        channels_r.extend([ch_r_i])
                        channels_t.extend([ch_t_i])
    else:
        channels_r = ch_r
        channels_t = ch_t
        
    return(channels_r, channels_t)