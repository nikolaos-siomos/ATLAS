#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:22:31 2022

@author: nick
"""

import numpy as np
import sys

def trim_channels(cfg, sig, shots, channel_info, meas_type):
    
    """Channels are selected based on the channel_id variable in the configuration
    file. All IDs in the config_file must be correspondent to the IDs in
    the raw files. If an unkown ID is included in the config_file then an error
    is raised. For polarization calibration measurements, channels that are 
    neither co- nor cross- polar will be automatically removed"""
    

    channel_subtype_cfg = cfg.channels.channel_subtype.values

    channel_ind_cfg = cfg.channels.index.values
            
    if meas_type == 'pcb':
        
        channel_ind_cfg = channel_ind_cfg[(channel_subtype_cfg == 'r') | 
                                          (channel_subtype_cfg == 't')]
        
    channel_ind_raw = channel_info.index.values
    
    if any(channel_ind_cfg_i not in channel_ind_raw 
           for channel_ind_cfg_i in channel_ind_cfg):
        print("-- Warning: Some of the recorder channel IDs defined in the configuration file do not match with the respective IDs in the raw files. These channels will be skipped! ")
    
    channel_ind_com = [channel_ind_cfg_i for channel_ind_cfg_i in channel_ind_cfg 
                       if channel_ind_cfg_i in channel_ind_raw]

    
    sig = sig.loc[dict(channel = channel_ind_com)]
    
    shots = shots.loc[dict(channel = channel_ind_com)]
    
    channel_info = channel_info.loc[channel_ind_com,:]
    
    cfg.channels = cfg.channels.loc[channel_ind_com,:]
        
    return(sig, shots, channel_info, cfg)

def merge_config(cfg, meas_info, channel_info):

    for key in meas_info.index:
        if key not in cfg.meas.index:
            cfg.meas.loc[key] = meas_info.loc[key]    

    for key in channel_info.columns:
        if key not in cfg.channels.columns:
            cfg.channels.loc[:,key] = channel_info.loc[:,key].values    
    
    return(cfg)
    
def fill_defaults(cfg):
    
    if 'channel_bandwidth' not in cfg.channels:
        cfg.channels.loc[:,'channel_bandwidth'] = 1.  
    
    if 'zenith_angle' not in cfg.meas:
        cfg.meas['zenith_angle'] = 0.  
        
    if 'azimuth_angle' not in cfg.meas:
        cfg.meas['azimuth_angle'] = 0.  
        
    cfg.channels.loc[:,'laser_repetition_rate'] = cfg.meas['laser_A_repetition_rate']
    if (cfg.channels.loc[:,'laser'] == 2).any():
        cfg.channels.loc[:,'laser_repetition_rate'][cfg.channels.loc[:,'laser'] == 2] = cfg.meas['laser_B_repetition_rate']
    if (cfg.channels.loc[:,'laser'] == 3).any():
        cfg.channels.loc[:,'laser_repetition_rate'][cfg.channels.loc[:,'laser'] == 3] = cfg.meas['laser_C_repetition_rate']
        
    acquisition_type = np.empty(cfg.channels.acquisition_mode.size,dtype=object)
    acquisition_type[cfg.channels.acquisition_mode.values.astype(int) == 0] = 'a'
    acquisition_type[cfg.channels.acquisition_mode.values.astype(int) == 1] = 'p' 
  
    cfg.channels.loc[:,'acquisition_type'] = acquisition_type  

    cfg.channels.loc[:,'background_low'] = float(cfg.meas.loc['altitude']) + cfg.channels.loc[:,'background_low_bin'].astype(float) * cfg.channels.loc[:,'range_resolution'].astype(float)
    cfg.channels.loc[:,'background_high'] = float(cfg.meas.loc['altitude']) + cfg.channels.loc[:,'background_high_bin'].astype(float) * cfg.channels.loc[:,'range_resolution'].astype(float)
        
    return(cfg)

def unit_conv_bits_to_mV(channel_info, signal, shots):

    """Converts analog signals from bits to mV"""
    
    if len(signal) > 0:
        
        mask_an = channel_info.acquisition_mode.values == 0
        
        channel_id_an = channel_info.index.values[mask_an]
        
        data_acquisition_range = channel_info.data_acquisition_range
        
        analog_to_digital_resolution = channel_info.analog_to_digital_resolution
        
        for ch in channel_id_an:
            ch_d = dict(channel = ch)
            # analog conversion (to mV)
            signal.loc[ch_d] = signal.loc[ch_d]*data_acquisition_range.loc[ch]/(shots.loc[ch_d]*(np.power(2,analog_to_digital_resolution.loc[ch])-1.))
    
    return(signal) 

def screen_low_shots(time_info, channel_info, signal, shots):

    """Replaces values in all bins with nans if the numer of shots is
    lower than 10% of the maximum """
    
    
    if len(signal) > 0:
        time = signal.copy().time.values
                
        channel_id = channel_info.index.values
        
        time_ls = []
                
        for ch in channel_id:
            ch_d = dict(channel = ch)
            
            shots_ch = shots.loc[ch_d].values.astype(float)
            
            if not np.isnan(shots_ch).all():
                
                mask_t = shots_ch < 0.9 * np.nanmax(shots_ch)
                            
                time_ls.extend(time[mask_t])
    
                if len(time[mask_t]) > 0:                
                    print('-- Warning: The following files have less shots than 90% of the naximum number of shots encountered and will be screened out:')
                    for t in time[mask_t]:
                        ind_t = np.where(t == time)[0][0]
                        print(f'    Filename: {time_info.filename[ind_t]} | Channel: {ch}')
                    
                    sel = dict(time = time[mask_t], channel = ch)
                    
                    signal.loc[sel] = np.nan
            
                    shots.loc[sel] = np.nan
        
        # if le


        
            
    return(signal) 