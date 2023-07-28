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

def merge_config(cfg, system_info, channel_info):
    
    """ The optional fields are that are not provided in the configuration file
    get values from the headers of the raw files. An error is raised if the
    information is not available even from the raw file, which indicates either
    a corrupted input file or a bug in the parser"""

    for key in system_info.index:
        
        if key not in cfg.system.index and key not in system_info.index :
            raise Exception(f"-- Error: The optional field {key} was not provided in the configuration file and was also not found in the raw files. Please verify that the raw file format is correct and contact the developing team if necessary")
        
        elif key not in cfg.system.index and key in system_info.index:
            cfg.system.loc[key] = system_info.loc[key]    

    for key in channel_info.columns:
        
        if key not in cfg.channels.columns and key not in channel_info.columns:
            raise Exception(f"-- Error: The optional field {key} of the Channels Section was not provided in the configuration file and was also not found in the raw files. Please verify that the raw file format is correct and contact the developing team if necessary")
        
        elif key not in cfg.channels.columns and key in channel_info.columns:
            cfg.channels.loc[:,key] = channel_info.loc[:,key].values    
    
    return(cfg)
    
def fill_defaults(cfg):
    
    """The partially optional values are filled with their default values 
    directly at the cfg object. 
    In addition:
        -- the aquisition_mode is converted from integer to string
        -- the laser repetition rate is assigned to every channel
        -- the background_low and background_high bins which are need to ensure
           compatibility with the SCC are added to the cdf object
    """
    
    # Channel defaults
    empty =  cfg.channels.index.size * ['_']

    # Acquisition type
    mask_an = cfg.channels.acquisition_mode.values.astype(int) == 0
    mask_pc = cfg.channels.acquisition_mode.values.astype(int) == 1
    if 'acquisition_type' not in cfg.channels:
        cfg.channels.loc[:,'acquisition_type'] = empty
    cfg.channels.loc[:,'acquisition_type'][mask_an] = 'a'
    cfg.channels.loc[:,'acquisition_type'][mask_pc] = 'p'
    
    # Dead time check
    if 'dead_time' not in cfg.channels:
        cfg.channels.loc[:,'dead_time'] = empty
    mask_dt_an = (cfg.channels.loc[:,'dead_time'] == '_') & (cfg.channels.loc[:,'acquisition_type'] == 'a')
    mask_dt_pc = (cfg.channels.loc[:,'dead_time'] == '_') & (cfg.channels.loc[:,'acquisition_type'] == 'p')
    cfg.channels.loc[:,'dead_time'][mask_dt_an] = np.nan
    cfg.channels.loc[:,'dead_time'][mask_dt_pc] = 3.70
        
    # Trigger delay check
    if 'daq_trigger_offset' not in cfg.channels:
        cfg.channels.loc[:,'daq_trigger_offset'] = empty
    mask_trigger = (cfg.channels.loc[:,'daq_trigger_offset'] == '_')
    cfg.channels.loc[:,'daq_trigger_offset'][mask_trigger] = 0   
    
    # Dead time correction type check
    if 'dead_time_correction_type' not in cfg.channels:
        cfg.channels.loc[:,'dead_time_correction_type'] = empty   
    mask_dt_cor_an = (cfg.channels.loc[:,'dead_time_correction_type'] == '_') & (cfg.channels.loc[:,'acquisition_type'] == 'a')
    mask_dt_cor_pc = (cfg.channels.loc[:,'dead_time_correction_type'] == '_') & (cfg.channels.loc[:,'acquisition_type'] == 'p')
    
    cfg.channels.loc[:,'dead_time_correction_type'][mask_dt_cor_an] = np.nan
    cfg.channels.loc[:,'dead_time_correction_type'][mask_dt_cor_pc] = 0
        
    # Background Low Bin and Height
    if 'background_low_bin' not in cfg.channels.columns.values:
        cfg.channels.loc[:,'background_low_bin'] = empty
    mask_lbin_tr = (cfg.channels.loc[:,'background_low_bin'] == '_') & (cfg.channels.loc[:,'daq_trigger_offset'].astype(int) <= -400)
    mask_lbin_fr = (cfg.channels.loc[:,'background_low_bin'] == '_') & (cfg.channels.loc[:,'daq_trigger_offset'].astype(int) >  -400)
    cfg.channels.loc[:,'background_low_bin'][mask_lbin_tr] = 100
    cfg.channels.loc[:,'background_low_bin'][mask_lbin_fr] = cfg.channels.loc[:,'bins'][mask_lbin_fr] - 600
    cfg.channels.loc[:,'background_low'] = float(cfg.system.loc['altitude']) + np.cos(np.deg2rad(float(cfg.system['zenith_angle']))) * cfg.channels.loc[:,'background_low_bin'].astype(float) * cfg.channels.loc[:,'range_resolution'].astype(float)
  
    # Background High Bin
    if 'background_high_bin' not in cfg.channels.columns.values:
        cfg.channels.loc[:,'background_high_bin'] = empty
    mask_ubin_tr = (cfg.channels.loc[:,'background_high_bin'] == '_') & (cfg.channels.loc[:,'daq_trigger_offset'].astype(int) <= -400)
    mask_ubin_fr = (cfg.channels.loc[:,'background_high_bin'] == '_') & (cfg.channels.loc[:,'daq_trigger_offset'].astype(int) >  -400)
    cfg.channels.loc[:,'background_high_bin'][mask_ubin_tr] = -cfg.channels.loc[:,'daq_trigger_offset'][mask_ubin_tr].astype(int) - 100
    cfg.channels.loc[:,'background_high_bin'][mask_ubin_fr] = cfg.channels.loc[:,'bins'][mask_ubin_fr] - 100
    cfg.channels.loc[:,'background_high'] = float(cfg.system.loc['altitude']) + np.cos(np.deg2rad(float(cfg.system['zenith_angle']))) * cfg.channels.loc[:,'background_high_bin'].astype(float) * cfg.channels.loc[:,'range_resolution'].astype(float)
                        
    # Laser repetiotion rate - assign each laser value to the respective channel 
    mask_laser_A = (cfg.channels.loc[:,'laser'] == 1)
    mask_laser_B = (cfg.channels.loc[:,'laser'] == 2)
    mask_laser_C = (cfg.channels.loc[:,'laser'] == 3)
    if 'laser_repetition_rate' not in cfg.channels.columns.values:
        cfg.channels.loc[:,'laser_repetition_rate'] = empty
        cfg.channels.loc[:,'laser_repetition_rate'][mask_laser_A] = cfg.system['laser_A_repetition_rate']
        if 'laser_B_repetition_rate' in cfg.system.index:
            cfg.channels.loc[:,'laser_repetition_rate'][mask_laser_B] = cfg.system['laser_B_repetition_rate']
        if 'laser_C_repetition_rate' in cfg.system.index:
            cfg.channels.loc[:,'laser_repetition_rate'][mask_laser_C] = cfg.system['laser_C_repetition_rate']

    # Emitted wavelength
    if 'emitted_wavelength' not in cfg.channels:
        cfg.channels.loc[:,'emitted_wavelength'] = empty
    mask_uv = (cfg.channels.loc[:,'emitted_wavelength'] == '_') & \
        (cfg.channels.loc[:,'detected_wavelength'] >= 340.) & \
            (cfg.channels.loc[:,'detected_wavelength'] < 520.)
    mask_vs = (cfg.channels.loc[:,'emitted_wavelength'] == '_') & \
        (cfg.channels.loc[:,'detected_wavelength'] >= 520.) & \
            (cfg.channels.loc[:,'detected_wavelength'] < 1000.)    
    mask_ir = (cfg.channels.loc[:,'emitted_wavelength'] == '_') & \
        (cfg.channels.loc[:,'detected_wavelength'] >= 1000.)
    cfg.channels.loc[:,'emitted_wavelength'][mask_uv] = 354.717
    cfg.channels.loc[:,'emitted_wavelength'][mask_vs] = 532.075
    cfg.channels.loc[:,'emitted_wavelength'][mask_ir] = 1064.150

    # Channel bandwidth
    if 'channel_bandwidth' not in cfg.channels:
        cfg.channels.loc[:,'channel_bandwidth'] = empty
    mask_bdw = (cfg.channels.loc[:,'channel_bandwidth'] == '_')
    cfg.channels.loc[:,'channel_bandwidth'][mask_bdw] = 1.   

    # for ch in cfg.channels.index:
    #     print(cfg.channels.loc[ch,:])
    #     raise Exception()
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