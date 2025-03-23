#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 22:34:11 2022

@author: nick
"""

import numpy as np

telescope_map = {'n' : 'Near Range',
                 'f' : 'Far Range',
                 'x' : ''}

mode_map = {'a' : 'analog',
            'p' : 'photon'}

type_map = {'p' : 'Co-polar',
            'c' : 'Cross-polar',
            't' : 'Total',
            'v' : 'Vibrational Raman',
            'r' : 'Rotational Raman',
            'a' : 'Cabannes',
            'f' : 'Fluorescence'}

subtype_map = {'r' : 'Reflected',
               't' : 'Transmitted',
               'n' : 'N2',
               'o' : 'O2',
               'w' : 'H20',
               'c' : 'CH4',
               'l' : 'Low hat',
               'h' : 'High hat',
               'a' : 'Mie',
               'm' : 'Molecular',
               'b' : 'Broadband',
               's' : 'Spectral',
               'x' : ''}

def rayleigh(channel, metadata, norm_region):
    
    start_date = metadata['start_date']
    start_time = metadata['start_time'] 
    stop_time = metadata['stop_time'] 

    wave = metadata['dwl'].copy().loc[channel].values
    lidar_name = metadata['lidar_name'] 
    station_name = metadata['station_name'] 
    meas_id = metadata['meas_id'] 
    st_name = metadata['rs_station_name'] 
    rs_start_date = metadata['rs_start_date']
    rs_start_time = metadata['rs_start_time'] 
    wmo_id = metadata['wmo_id'] 
        
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    start = f'{start_time[:2]}:{start_time[2:4]}'
    
    duration = np.round(float(stop_time[:2]) - float(start_time[:2]) + float(stop_time[2:4]) / 60. - float(start_time[2:4]) / 60., decimals = 1)
    wave = int(np.round(metadata['dwl'].copy().loc[channel].values,decimals=0))
    
    rs_date = f'{rs_start_date[6:]}.{rs_start_date[4:6]}.{rs_start_date[:4]}'    
    rs_start = f'{rs_start_time[:2]}:{rs_start_time[2:4]}'
    
    ray_l = np.round(norm_region[0], decimals = 1)
    ray_u = np.round(norm_region[1], decimals = 1)
    
    expression = f'{telescope_map[channel[4]]} {subtype_map[channel[7]]} {type_map[channel[5]]}'
    
    header = f'station ID = {meas_id[8:11]}\n'+\
        f'system = {lidar_name} - {station_name}\n' +\
            f'signal = {wave}, {expression}, {mode_map[channel[6]]}, dark-subtracted\n' +\
                f'date of measurement, time, duration of measurement = {date}, {start}UTC, {duration} s\n' +\
                    f'location, WMO radiosonde station ID, date of radiosonde = {st_name}, {wmo_id}, {rs_date}, {rs_start}UT\n' +\
                        f'lower and upper Rayleigh height limits = {ray_l}, {ray_u}\n' +\
                            'range, attnRayleighBSC, RangeCorrectedSignal'
                        
    return(header)

def telecover(channel, metadata, iters, extra_sec):
    
    lidar_name = metadata['lidar_name'] 
    station_name = metadata['station_name'] 
    meas_id = metadata['meas_id'] 
    
    start_date = metadata['start_date']
    start_time = metadata['start_time'] 
    stop_time = metadata['start_time'] 
    
    duration = np.round(float(stop_time[:2]) - float(start_time[:2]) + float(stop_time[2:4]) / 60. - float(start_time[2:4]) / 60., decimals = 1)
    
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    start = f'{start_time[:2]}:{start_time[2:4]}'
    
    wave = int(np.round(metadata['dwl'].copy().loc[channel].values,decimals=0))
    
    expression = f'{telescope_map[channel[4]]} {type_map[channel[5]]} {subtype_map[channel[7]]}'
    
    secs = [key for key in extra_sec.keys()]
    
    extra = [f'{key}{iters+1}' for key in extra_sec.keys() if extra_sec[key]]
    
    iter_num = np.arange(1,iters+1,1)
    
    combo = []
    for num in iter_num:
        for sec in secs:
            combo.append(f'{sec}{num}')
    
    sec_text = ', '.join(combo+extra) 
    
    header = f'station ID = {meas_id[8:11]}\n'+\
        f'system = {lidar_name} - {station_name}\n' +\
            f'signal = {wave}, {expression}, {mode_map[channel[6]]}, dark-subtracted\n' +\
                f'date of measurement, time, duration of measurement = {date}, {start}UTC, {duration} s\n' +\
                    f'range, {sec_text}'
                        
    return(header)

def polarisation_calibration(channel_r, channel_t, metadata, K, G_R, G_T, H_R, H_T):
        
    start_date_ray = metadata['start_date_ray']
    start_time_ray = metadata['start_time_ray'] 
    stop_time_ray = metadata['stop_time_ray'] 

    duration_ray = np.round(float(stop_time_ray[:2]) - float(start_time_ray[:2]) + float(stop_time_ray[2:4]) / 60. - float(start_time_ray[2:4]) / 60., decimals = 1)
    
    start_date_cal = metadata['start_date_cal']
    start_time_cal = metadata['start_time_cal'] 
    stop_time_cal = metadata['stop_time_cal'] 

    duration_cal = np.round(float(stop_time_cal[:2]) - float(start_time_cal[:2]) + float(stop_time_cal[2:4]) / 60. - float(start_time_cal[2:4]) / 60., decimals = 1)
    
    lidar_name = metadata['lidar_name'] 
    station_name = metadata['station_name'] 
    meas_id = metadata['meas_id'] 
    
    wave = int(np.round(metadata['dwl_cal'].copy().loc[channel_r].values,decimals=0))
    
    expression_r = f'R: {telescope_map[channel_r[4]]} {type_map[channel_r[5]]} {mode_map[channel_r[6]]}'
    expression_t = f'T: {telescope_map[channel_t[4]]} {type_map[channel_t[5]]} {mode_map[channel_t[6]]}'
    
    print()
    header = f'station ID = {meas_id[8:11]}\n'+\
        f'system = {lidar_name} - {station_name}\n' +\
            f'signal = {wave}, {expression_r}, {expression_t}, dark-subtracted\n' +\
                f'date of calibration measurement, time, duration of measurement = {start_date_cal}, {start_time_cal}UTC, {duration_cal} s\n' +\
                    f'date of Rayleigh measurement, time, duration of measurement = {start_date_ray}, {start_time_ray}UTC, {duration_ray} s\n' +\
                        f'GR, GT, HR, HT, K = {G_R} {G_T} {H_R} {H_T} {K}\n' +\
                            'range, ITplus45, IRplus45, ITminus45, IRminus45, ITRayleigh, IRRayleigh'
                        
    return(header)