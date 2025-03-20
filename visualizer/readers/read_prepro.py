#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:16:49 2025

@author: nikos
"""
import xarray as xr
import netCDF4 as nc
import numpy as np

def unpack(input_file):
    
    data = xr.open_dataset(input_file)

    profiles = dict()
    metadata = dict()
    
    # Extract Quicklook and Rayleigh parameters
    try: 
        sig = data.Range_Corrected_Signals
        sig = sig.copy().where(sig != nc.default_fillvals['f8'])

        profiles['sig'] = sig
        
        profiles['ranges'] = data.Range_levels
        profiles['heights'] = data.Height_levels
        
        # Extract date info
        metadata['start_date'] = data.RawData_Start_Date
        metadata['start_time'] = data.RawData_Start_Time_UT
        metadata['stop_time'] = data.RawData_Stop_Time_UT
        
    except: pass
        
    # Extract Quicklook parameters
    try: 
        metadata['Raw_Data_Start_Time'] = data['Raw_Data_Start_Time']
        metadata['Raw_Data_Stop_Time'] = data['Raw_Data_Stop_Time']
        metadata['meas_id'] = data.Measurement_ID

    except: pass

    # Extract Rayleigh parameters
    try: 
        atb = data.Attenuated_Backscatter
        
        profiles['atb'] = atb
    
    except: pass
    
    # Extract Quadrant Telecover parameters
    try: 
        sig_n = data.Range_Corrected_Signals_North_Sector
        sig_e = data.Range_Corrected_Signals_East_Sector
        sig_s = data.Range_Corrected_Signals_South_Sector
        sig_w = data.Range_Corrected_Signals_West_Sector
        
        sig_n = sig_n.copy().where(sig_n != nc.default_fillvals['f8'])
        sig_e = sig_e.copy().where(sig_e != nc.default_fillvals['f8'])
        sig_s = sig_s.copy().where(sig_s != nc.default_fillvals['f8'])
        sig_w = sig_w.copy().where(sig_w != nc.default_fillvals['f8'])
        
        ranges = data.Range_levels
        heights = data.Height_levels
      
        start_time_n_dt = data['Raw_Data_Start_Time_North_Sector']
        start_time_e_dt = data['Raw_Data_Start_Time_East_Sector']
        start_time_s_dt = data['Raw_Data_Start_Time_South_Sector']
        start_time_w_dt = data['Raw_Data_Start_Time_West_Sector']

        stop_time_n_dt = data['Raw_Data_Stop_Time_North_Sector']
        stop_time_e_dt = data['Raw_Data_Stop_Time_East_Sector']
        stop_time_s_dt = data['Raw_Data_Stop_Time_South_Sector']
        stop_time_w_dt = data['Raw_Data_Stop_Time_West_Sector']
                
        sampling_time_per_quadrant = np.min([
            np.min((stop_time_n_dt - start_time_n_dt).values),
            np.min((stop_time_e_dt - start_time_e_dt).values),
            np.min((stop_time_s_dt - start_time_s_dt).values),
            np.min((stop_time_w_dt - start_time_w_dt).values)
                ])
        
        profiles['sig_n'] = sig_n
        profiles['sig_e'] = sig_e
        profiles['sig_s'] = sig_s
        profiles['sig_w'] = sig_w
        
        profiles['ranges'] = ranges
        profiles['heights'] = heights

        metadata['start_date'] = data.RawData_Start_Date
        metadata['start_time'] = data.RawData_Start_Time_UT
        metadata['stop_time'] = data.RawData_Stop_Time_UT
        
        metadata['sampling_time_per_quadrant'] = sampling_time_per_quadrant        

    except: pass

    # Extract Ring Telecover parameters
    try: 
        sig_i = data.Range_Corrected_Signals_Inner_Sector
        sig_o = data.Range_Corrected_Signals_Outer_Sector
        
        sig_i = sig_i.copy().where(sig_i != nc.default_fillvals['f8'])
        sig_o = sig_o.copy().where(sig_o != nc.default_fillvals['f8'])
           
        start_time_i_dt = data['Raw_Data_Start_Time_Inner_Sector']
        start_time_o_dt = data['Raw_Data_Start_Time_Outer_Sector']

        stop_time_i_dt = data['Raw_Data_Stop_Time_Inner_Sector']
        stop_time_o_dt = data['Raw_Data_Stop_Time_Outer_Sector']
                
        sampling_time_per_ring = np.min([
            np.min((stop_time_i_dt - start_time_i_dt).values),
            np.min((stop_time_o_dt - start_time_o_dt).values)
            ])
        
        profiles['sig_i'] = sig_i
        profiles['sig_o'] = sig_o

        metadata['start_date'] = data.RawData_Start_Date
        metadata['start_time'] = data.RawData_Start_Time_UT
        metadata['stop_time'] = data.RawData_Stop_Time_UT
        
        metadata['sampling_time_per_ring'] = sampling_time_per_ring  
                     
    except: pass
    
    # Extract Pol. Cal. parameters
    try:
        sig_ray = data.Range_Corrected_Signals_Rayleigh
        sig_ray = sig_ray.copy().where(sig_ray != nc.default_fillvals['f8'])
        
        sig_m45 = data.Range_Corrected_Signals_minus_45
        sig_m45 = sig_m45.copy().where(sig_m45 != nc.default_fillvals['f8'])\
            .mean(dim = 'time_m45')
        
        sig_p45 = data.Range_Corrected_Signals_plus_45
        sig_p45 = sig_p45.copy().where(sig_p45 != nc.default_fillvals['f8'])\
            .mean(dim = 'time_p45')
                    
        profiles['sig_ray'] = sig_ray
        profiles['sig_m45'] = sig_m45
        profiles['sig_p45'] = sig_p45
        
        profiles['mldr'] = data.Molecular_Linear_Depolarization_Ratio
        
        profiles['ranges_ray'] = data.Range_levels_Rayleigh
        profiles['ranges_cal'] = data.Range_levels_Calibration

        profiles['heights_ray'] = data.Height_levels_Rayleigh
        profiles['heights_cal'] = data.Height_levels_Calibration

        metadata['start_date_ray'] = data.RawData_Start_Date_Rayleigh
        metadata['start_time_ray'] = data.RawData_Start_Time_UT_Rayleigh
        metadata['stop_time_ray'] = data.RawData_Stop_Time_UT_Rayleigh
        
        metadata['start_date_cal'] = data.RawData_Start_Date_Calibration
        metadata['start_time_cal'] = data.RawData_Start_Time_UT_Calibration
        metadata['stop_time_cal'] = data.RawData_Stop_Time_UT_Calibration
        
        metadata['start_date'] = data.RawData_Start_Date_Calibration
        metadata['start_time'] = data.RawData_Start_Time_UT_Calibration
        metadata['stop_time'] = data.RawData_Stop_Time_UT_Calibration
        
        metadata['dwl_ray'] = data.Detected_Wavelength_Rayleigh
        metadata['ewl_ray'] = data.Emitted_Wavelength_Rayleigh
        metadata['bdw_ray'] = data.Channel_Bandwidth_Rayleigh 
    
        metadata['dwl_cal'] = data.Detected_Wavelength_Calibration
        metadata['ewl_cal'] = data.Emitted_Wavelength_Calibration
        metadata['bdw_cal'] = data.Channel_Bandwidth_Calibration    
        
        metadata['laser_pointing_angle_ray'] = data.Laser_Pointing_Angle_Rayleigh
        metadata['laser_pointing_angle_cal'] = data.Laser_Pointing_Angle_Calibration
    
        metadata['meas_id'] = data.Measurement_ID_Calibration
        
        try: metadata['scc_channel_id'] = data.channel_ID_Calibration
        except: metadata['scc_channel_id'] = ''

    except: pass
        
    # Extract common Rayleigh - Telecover parameters
    try:
        # Extract Channel info
        metadata['dwl'] = data.Detected_Wavelength
        metadata['ewl'] = data.Emitted_Wavelength
        metadata['bdw'] = data.Channel_Bandwidth    
        metadata['dead_time'] = data.Dead_Time
        metadata['daq_trigger_offset'] = data.DAQ_Trigger_Offset
        metadata['background_low_bin'] = data.Background_Low_Bin
        metadata['background_high_bin'] = data.Background_High_Bin
        metadata['raw_data_range_resolution'] = data.Raw_Data_Range_Resolution
        
        metadata = atlas_to_scc_triggering(metadata)
        
        # Extract Other info
        metadata['laser_pointing_angle'] = data.Laser_Pointing_Angle
        metadata['meas_id'] = data.Measurement_ID
        
        try: metadata['scc_channel_id'] = data.channel_ID
        except: metadata['scc_channel_id'] = ''

    except: pass

    # Extract channels
    metadata['atlas_channel_id'] = data.channel.values
            
    # Extract SCC info
    metadata['station_id'] = data.Station_ID.lower()
        
    metadata['station_name'] = data.Station_Name

    try: metadata['lidar_id'] = data.Lidar_ID
    except: metadata['lidar_id'] = ''
    
    try: metadata['lidar_name'] = data.Lidar_Name
    except: metadata['lidar_name'] = ''
    
    try: metadata['version_id'] = data.Version_ID
    except: metadata['version_id'] = ''
    
    try: metadata['config_id'] = data.Configuration_ID
    except: metadata['config_id'] = ''
    
    try: metadata['config_name'] = data.Configuration_Name
    except: metadata['config_name'] = ''
    
    # Exctract molecular data
    try: metadata['mol_method'] = data.Molecular_atmosphere_method
    except: metadata['mol_method'] = ''
    
    try: metadata['rs_station_name'] = data.Sounding_Station_Name
    except: metadata['rs_station_name'] = ''
    
    try: metadata['rs_start_date'] = data.Sounding_Start_Date
    except: metadata['rs_start_date'] = ''
    
    try: metadata['rs_start_time'] = data.Sounding_Start_Time_UT
    except: metadata['rs_start_time'] = ''
  
    try: metadata['wmo_id'] = data.WMO_Station_Number
    except: metadata['wmo_id'] = ''

    try: metadata['wban_id'] = data.WBAN_Station_Number
    except: metadata['wban_id']  = ''

    return(profiles, metadata)

def atlas_to_scc_triggering(metadata):
    background_low_bin = metadata['background_low_bin']
    background_high_bin = metadata['background_high_bin']
    daq_trigger_offset = metadata['daq_trigger_offset']
    raw_data_range_resolution = metadata['raw_data_range_resolution']
    
    background_mode = np.nan * daq_trigger_offset.copy().astype(object)
    background_low = np.nan * daq_trigger_offset.copy().astype(object)
    background_high = np.nan * daq_trigger_offset.copy().astype(object)
    first_signal_rangebin = np.nan * daq_trigger_offset.copy().astype(object)
    trigger_delay = np.nan * daq_trigger_offset.copy().astype(object)

    mask_pretrg = daq_trigger_offset.values < -50
    mask_isdelay = (daq_trigger_offset.values >= -50) & (daq_trigger_offset.values < 0)

    if mask_pretrg.any():
        background_mode[mask_pretrg] = 'Pre-Trigger'
        background_low[mask_pretrg] = background_low_bin[mask_pretrg]
        background_high[mask_pretrg] = background_high_bin[mask_pretrg]
        first_signal_rangebin[mask_pretrg] = -daq_trigger_offset[mask_pretrg]
        trigger_delay[mask_pretrg] = -999.
    if (~mask_pretrg).any():
        background_mode[~mask_pretrg] = 'Far Field'
        background_mode[~mask_pretrg] = raw_data_range_resolution[~mask_pretrg] * background_low_bin[~mask_pretrg]
        background_high[~mask_pretrg] = raw_data_range_resolution[~mask_pretrg] * background_high_bin[~mask_pretrg]
        if mask_isdelay.any():
            first_signal_rangebin[mask_isdelay] = -daq_trigger_offset[mask_isdelay] 
            trigger_delay[mask_isdelay] = -999.
        if (~mask_isdelay).any():
            first_signal_rangebin[~mask_isdelay] = -999.
            trigger_delay[~mask_isdelay] = daq_trigger_offset[~mask_isdelay] * raw_data_range_resolution[~mask_isdelay] * 20. / 3.

    metadata['background_mode'] = background_mode
    metadata['background_low'] = background_low
    metadata['background_high'] = background_high
    metadata['first_signal_rangebin'] = first_signal_rangebin
    metadata['trigger_delay'] = trigger_delay

    return(metadata)
