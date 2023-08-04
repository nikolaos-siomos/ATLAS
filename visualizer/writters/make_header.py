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

def rayleigh(start_date, start_time, start_time_sec, stop_time_sec,
             wave, lidar, loc, meas_id, channel, norm_region,
             st_name, rs_start_date, rs_start_time, wmo_id, wban_id):
    
    duration = np.sum(stop_time_sec - start_time_sec)
    
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    start = f'{start_time[:2]}:{start_time[2:4]}'
    
    wv = int(np.round(wave,decimals=0))
    
    rs_date = f'{rs_start_date[6:]}.{rs_start_date[4:6]}.{rs_start_date[:4]}'    
    rs_start = f'{rs_start_time[:2]}:{rs_start_time[2:4]}'
    
    ray_l = np.round(norm_region[0], decimals = 1)
    ray_u = np.round(norm_region[1], decimals = 1)
    
    expression = f'{telescope_map[channel[4]]} {subtype_map[channel[7]]} {type_map[channel[5]]}'
    
    header = f'station ID = {meas_id[8:11]}\n'+\
        f'system = {lidar} - {loc}\n' +\
            f'signal = {wv}, {expression}, {mode_map[channel[6]]}, dark-subtracted\n' +\
                f'date of measurement, time, duration of measurement = {date}, {start}UTC, {duration} s\n' +\
                    f'location, WMO radiosonde station ID, date of radiosonde = {st_name}, {wmo_id}, {rs_date}, {rs_start}UT\n' +\
                        f'lower and upper Rayleigh height limits = {ray_l}, {ray_u}\n' +\
                            'range, attnRayleighBSC, RangeCorrectedSignal'
                        
    return(header)

def telecover(start_date, start_time, sampling,
              wave, lidar, loc, meas_id, channel, iters, extra_sec):
    
    duration = iters * sampling
    
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    start = f'{start_time[:2]}:{start_time[2:4]}'
    
    wv = int(np.round(wave,decimals=0))
    
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
        f'system = {lidar} - {loc}\n' +\
            f'signal = {wv}, {expression}, {mode_map[channel[6]]}, dark-subtracted\n' +\
                f'date of measurement, time, duration of measurement = {date}, {start}UTC, {duration} s\n' +\
                    f'range, {sec_text}'
                        
    return(header)

def polarisation_calibration(cal_start_date, cal_start_time, cal_duration,
                             ray_start_date, ray_start_time, ray_duration,
                             G_R, G_T, H_R, H_T, K,
                             wave, lidar, loc, meas_id, channel_r, channel_t):
        
    cal_date = f'{cal_start_date[6:]}.{cal_start_date[4:6]}.{cal_start_date[:4]}'
    cal_start = f'{cal_start_time[:2]}:{cal_start_time[2:4]}'

    ray_date = f'{ray_start_date[6:]}.{ray_start_date[4:6]}.{ray_start_date[:4]}'
    ray_start = f'{ray_start_time[:2]}:{ray_start_time[2:4]}'
        
    wv = int(np.round(wave,decimals=0))
    
    expression_r = f'R: {telescope_map[channel_r[4]]} {type_map[channel_r[5]]} {mode_map[channel_r[6]]}'
    expression_t = f'T: {telescope_map[channel_r[4]]} {type_map[channel_r[5]]} {mode_map[channel_t[6]]}'
    
    print()
    header = f'station ID = {meas_id[8:11]}\n'+\
        f'system = {lidar} - {loc}\n' +\
            f'signal = {wv}, {expression_r}, {expression_t}, dark-subtracted\n' +\
                f'date of calibration measurement, time, duration of measurement = {cal_date}, {cal_start}UTC, {cal_duration} s\n' +\
                    f'date of Rayleigh measurement, time, duration of measurement = {ray_date}, {ray_start}UTC, {ray_duration} s\n' +\
                        f'GR, GT, HR, HT, K = {G_R} {G_T} {H_R} {H_T} {K}\n' +\
                            'range, ITplus45, IRplus45, ITminus45, IRminus45, ITRayleigh, IRRayleigh'
                        
    return(header)