#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:22:31 2022

@author: nick
"""

import numpy as np
import sys
import datetime as dt

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
    
    """ The optional fields that are not provided in the configuration file
    get values from the headers of the raw files. An error is raised if the
    information is not available even from the raw file, which indicates either
    a corrupted input file or a bug in the parser"""

    types_from_raw = get_from_raw('System')
    for key in types_from_raw.keys():
        
        if key not in cfg.system.index and key not in system_info.index :
            raise Exception(f"-- Error: The optional field {key} of the System Section was not provided in the configuration file and was also not found in the raw files. Please verify that the raw file format is correct and contact the developing team if necessary")
        
        elif key not in cfg.system.index and key in system_info.index:
            cfg.system.loc[key] = flexible_type(system_info.loc[key], 
                                                scalar_type =  types_from_raw[key])
            
        elif key in cfg.system.index:
            cfg.system.loc[key] = flexible_type(cfg.system.loc[key], 
                                                scalar_type =  types_from_raw[key])
            

    types_from_raw = get_from_raw('Channels')
    for key in types_from_raw.keys():
        
        if key not in cfg.channels.columns and key not in channel_info.columns:
            raise Exception(f"-- Error: The optional field {key} of the Channels Section was not provided in the configuration file and was also not found in the raw files. Please verify that the raw file format is correct and contact the developing team if necessary")
        
        elif key not in cfg.channels.columns and key in channel_info.columns:
            cfg.channels.loc[:,key] = channel_info.loc[:,key].values.astype(types_from_raw[key])   
            
        elif key in cfg.channels.columns:
            mask_empty = (cfg.channels.loc[:,key] == '_')
            cfg.channels.loc[mask_empty,key] = channel_info.loc[mask_empty,key].values   
            cfg.channels.loc[:,key] = cfg.channels.loc[:,key].astype(types_from_raw[key])
    
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
    mask_lbin_ng = (cfg.channels.loc[:,'background_low_bin'].astype(int) < 0)
    cfg.channels.loc[:,'background_low_bin'][mask_lbin_ng] = (cfg.channels.loc[:,'bins'][mask_lbin_ng] + cfg.channels.loc[:,'background_low_bin'][mask_lbin_ng].astype(int)).astype(object)
  
    # Background High Bin
    if 'background_high_bin' not in cfg.channels.columns.values:
        cfg.channels.loc[:,'background_high_bin'] = empty
    mask_ubin_tr = (cfg.channels.loc[:,'background_high_bin'] == '_') & (cfg.channels.loc[:,'daq_trigger_offset'].astype(int) <= -400)
    mask_ubin_fr = (cfg.channels.loc[:,'background_high_bin'] == '_') & (cfg.channels.loc[:,'daq_trigger_offset'].astype(int) >  -400)
    cfg.channels.loc[:,'background_high_bin'][mask_ubin_tr] = -cfg.channels.loc[:,'daq_trigger_offset'][mask_ubin_tr].astype(int) - 100
    cfg.channels.loc[:,'background_high_bin'][mask_ubin_fr] = cfg.channels.loc[:,'bins'][mask_ubin_fr] - 100
    mask_ubin_ng = (cfg.channels.loc[:,'background_high_bin'].astype(float) < 0)
    cfg.channels.loc[:,'background_high_bin'][mask_ubin_ng] = (cfg.channels.loc[:,'bins'][mask_ubin_ng] + cfg.channels.loc[:,'background_high_bin'][mask_ubin_ng].astype(int)).astype(object)
                      
    # Laser repetiotion rate - assign each laser value to the respective channel 
    mask_laser_A = (cfg.channels.loc[:,'laser'] == 1)
    mask_laser_B = (cfg.channels.loc[:,'laser'] == 2)
    mask_laser_C = (cfg.channels.loc[:,'laser'] == 3)
    if 'laser_repetition_rate' not in cfg.channels.columns.values:
        cfg.channels.loc[:,'laser_repetition_rate'] = empty
    cfg.channels.loc[:,'laser_repetition_rate'][mask_laser_A] = cfg.system['laser_A_repetition_rate']
    if not np.isnan(cfg.system['laser_B_repetition_rate']):
        cfg.channels.loc[:,'laser_repetition_rate'][mask_laser_B] = cfg.system['laser_B_repetition_rate']
    if not np.isnan('laser_C_repetition_rate' in cfg.system.index):
        cfg.channels.loc[:,'laser_repetition_rate'][mask_laser_C] = cfg.system['laser_C_repetition_rate']

    # Emitted wavelength
    if 'emitted_wavelength' not in cfg.channels:
        cfg.channels.loc[:,'emitted_wavelength'] = empty
    em_wv = cfg.channels.loc[:,'emitted_wavelength'].copy()
    dt_wv = cfg.channels.loc[:,'detected_wavelength'].copy()
    mask_uv = (em_wv == '_') & (dt_wv >= 340.) & (dt_wv < 520.)
    mask_vs = (em_wv == '_') & (dt_wv >= 520.) & (dt_wv < 1000.)    
    mask_ir = (em_wv == '_') & (dt_wv >= 1000.)
    if mask_uv.any():
        emitted_uv = dt_wv[mask_uv].iloc[np.argmin(np.abs(dt_wv[mask_uv] -  355.))]
        cfg.channels.loc[:,'emitted_wavelength'][mask_uv] = emitted_uv
    if mask_vs.any():
        emitted_vs = dt_wv[mask_vs].iloc[np.argmin(np.abs(dt_wv[mask_vs] -  532.))]
        cfg.channels.loc[:,'emitted_wavelength'][mask_vs] = emitted_vs
    if mask_ir.any():
        emitted_ir = dt_wv[mask_ir].iloc[np.argmin(np.abs(dt_wv[mask_ir] - 1064.))]
        cfg.channels.loc[:,'emitted_wavelength'][mask_ir] = emitted_ir

    # Convert to the correct field type for the partially optional fields
    default_types = get_default_types()
    for key in default_types.keys():
        cfg.channels.loc[:,key] = cfg.channels.loc[:,key].astype(default_types[key])

    # Baground Low and High required by the SCC
    cfg.channels.loc[:,'background_low'] = cfg.system.loc['altitude'] + np.cos(np.deg2rad(cfg.system['zenith_angle'])) * cfg.channels.loc[:,'background_low_bin'].values * cfg.channels.loc[:,'range_resolution'].values
    cfg.channels.loc[:,'background_high'] = cfg.system.loc['altitude'] + np.cos(np.deg2rad(cfg.system['zenith_angle'])) * cfg.channels.loc[:,'background_high_bin'].values * cfg.channels.loc[:,'range_resolution'].values

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
                
                mask_t = (shots_ch < 0.9 * np.nanmax(shots_ch)) & \
                    (np.nanmax(shots_ch) > 20)
                
                time_ls.extend(time[mask_t])
    
                if len(time[mask_t]) > 0:                
                    print('-- Warning: The following files have less shots than 90% of the naximum number of shots encountered and will be screened out:')
                    for t in time[mask_t]:
                        ind_t = np.where(t == time)[0][0]
                        print(f'    Filename: {time_info.filename[ind_t]} | Channel: {ch}')
                    
                    sel = dict(time = time[mask_t], channel = ch)
                    
                    signal.loc[sel] = np.nan
            
                    shots.loc[sel] = np.nan        
            
    return(signal) 

def get_from_raw(section):
    
    if section == 'System':
        fields_from_raw = dict(altitude = 'float', 
                               latitude = 'float',
                               longitude = 'float', 
                               zenith_angle = 'float',
                               azimuth_angle = 'float',
                               laser_A_repetition_rate = 'float',
                               laser_B_repetition_rate = 'float',
                               laser_C_repetition_rate = 'float')
                        
    elif section == 'Channels':
        fields_from_raw = dict(acquisition_mode = 'int', 
                               laser = 'int',
                               bins = 'int', 
                               range_resolution = 'float', 
                               data_acquisition_range = 'float',
                               analog_to_digital_resolution = 'float', 
                               detected_wavelength = 'float',
                               channel_bandwidth = 'float')
        
    else:
        raise Exception(f"Section {section} not understood. Please select among System or Channels")
    
    return(fields_from_raw)

def get_default_types():
        
    default_types = dict(dead_time = 'float', 
                         daq_trigger_offset = 'int', 
                         background_low_bin = 'int',
                         background_high_bin = 'int',
                         acquisition_mode = 'str',     
                         dead_time_correction_type = 'float', 
                         emitted_wavelength = 'float',
                         laser_repetition_rate = 'float')
    
    return(default_types)


def flexible_type(val, scalar_type):
    
    if scalar_type == 'str':
        val = str(val)
        
    elif scalar_type == 'float':
        val = float(val)
    
    elif scalar_type == 'int':
        val = int(val)
        
    else:
        raise Exception(f"The provided scalar type {scalar_type} is not understood. Please select among str, int, float") 
      
    return(val)

def slice_in_time(sig_raw, shots, time_info, slice_reg):
    
    s_slice = slice_reg[0]
    e_slice = slice_reg[1]
    
    if s_slice != None and e_slice != None:
    
        s_date = time_info.index[0].date()
        e_date = time_info.index[1].date()
        
        s_year = s_date.year
        s_month = s_date.month
        s_day = s_date.day
        e_day = e_date.day
    
        if e_day - s_day > 1:
            raise Exception("-- Error: Slicing measurements that cover more than two consequetive days is not supported. You have to split the measurement manually for now")
      
        s_hour = int(s_slice[:2])
        s_mins = int(s_slice[2:])
        e_hour = int(e_slice[:2])
        e_mins = int(e_slice[2:])    
    
        s_dt = dt.datetime(s_year, s_month, s_day, s_hour, s_mins)
        if e_hour + e_mins / 60.  > s_hour + s_mins / 60.:
            e_dt = dt.datetime(s_year, s_month, s_day, e_hour, e_mins)
        else:
            e_dt = dt.datetime(s_year, s_month, s_day + 1, e_hour, e_mins)
    
        time_info = time_info.loc[s_dt:e_dt].copy()
        sig_raw = sig_raw.loc[s_dt:e_dt,:,:].copy()
        shots = shots.loc[s_dt:e_dt,:].copy()
        
    return(sig_raw, shots, time_info)