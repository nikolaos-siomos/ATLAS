#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 22:34:11 2022

@author: nick
"""

import numpy as np

def quicklook(start_date, start_time, end_time, lidar, channel, 
              zan, loc, smooth, sm_lims, sm_win, sm_expo):
    
    zan = np.round(float(zan), decimals = 1)
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
    
    title = f'{lidar} {loc} {channel} - Quicklook'+sm_part+'\n'+\
            f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+f'{zan}'+r'$^{o}$ off-zenith'
                        
    return(title)

def rayleigh(start_date, start_time, end_time, lidar, channel, 
             zan, loc, ewl, dwl, bdw, smooth, sm_lims, sm_win, sm_expo,
             mol_method, st_name, rs_start_date, rs_start_time, wmo_id, wban_id):
    
    zan = np.round(float(zan), decimals = 1)
    
    ewl = np.round(float(ewl), decimals = 2)

    dwl = np.round(float(dwl), decimals = 2)
    
    bdw = np.round(float(bdw), decimals = 2)
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    mol_part = mol_text(mol_method, st_name, wmo_id, wban_id, rs_start_date, rs_start_time)
    
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'

    title = f'{lidar} {loc} {channel} - Rayleigh Fit' + sm_part + '\n'+\
                f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+f'{zan}'+r'$^{o}$ off-zenith' + mol_part + '\n'+\
                    f'Emitted WL: {ewl}nm, Received WL: {dwl}nm, Bandwidth: {bdw}nm'
    
    
    return(title)

def telecover(start_date, start_time, end_time, lidar, channel, 
              zan, loc, iters, sampling, smooth, sm_lims, sm_win, sm_expo):
    
    zan = np.round(float(zan), decimals = 1)

    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
        
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    title = f'{lidar} {loc} {channel} - Telecover Test'+sm_part+'\n'+\
            f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+f'{zan}'+r'$^{o}$ off-zenith' + f' - Iterations: {iters}, Sampling Time per Sector: {sampling}s'
                
    return(title)

def polarization_calibration(start_date_cal, start_time_cal, end_time_cal,
                             start_date_ray, start_time_ray, end_time_ray,  
                             lidar, channel_r, channel_t, zan, loc, 
                             ewl, dwl, bdw, smooth, sm_lims, sm_win, sm_expo,
                             mol_method, st_name, rs_start_date, rs_start_time,
                             wmo_id, wban_id):
    
    zan = np.round(float(zan), decimals = 1)
    
    ewl = np.round(float(ewl), decimals = 2)

    dwl = np.round(float(dwl), decimals = 2)
    
    bdw = np.round(float(bdw), decimals = 2)
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    date_cal = f'{start_date_cal[6:]}.{start_date_cal[4:6]}.{start_date_cal[:4]}'
    
    start_cal = f'{start_time_cal[:2]}:{start_time_cal[2:4]}:{start_time_cal[4:6]}'

    end_cal = f'{end_time_cal[:2]}:{end_time_cal[2:4]}:{end_time_cal[4:6]}'

    date_ray = f'{start_date_ray[6:]}.{start_date_ray[4:6]}.{start_date_ray[:4]}'
    
    start_ray = f'{start_time_ray[:2]}:{start_time_ray[2:4]}:{start_time_ray[4:6]}'

    end_ray = f'{end_time_ray[:2]}:{end_time_ray[2:4]}:{end_time_ray[4:6]}'
         
    rs_date = f'{rs_start_date[6:]}.{rs_start_date[4:6]}.{rs_start_date[:4]}'
    
    rs_start = f'{rs_start_time[:2]}:{rs_start_time[2:4]}'
    
    mol_part = mol_text(mol_method, st_name, wmo_id, wban_id, rs_date, rs_start)
        
    title = f'{lidar} {loc} {channel_r} to {channel_t} - Pol. Calibration'+sm_part+'\n'+\
            f'Calibration on {date_cal} from {start_cal} to {end_cal} UTC, '+\
            r'$\nearrow$' + f' {zan}' + r'$^{o}$ off-zenith'+ '\n'+\
            f'Rayleigh on {date_ray} from {start_ray} to {end_ray} UTC, '+\
            r'$\nearrow$' + f' {zan}' + r'$^{o}$ off-zenith' + mol_part + '\n'+\
            f'Emitted WL: {ewl}nm, Received WL: {dwl}nm, Bandwidth: {bdw}nm'  

    return(title)

def intercomparison(start_date, start_time, end_time, lidar_1, lidar_2, 
                    channel_1, channel_2, zan, lat, lon, elv):
    
    zan = np.round(float(zan), decimals = 1)
    
    lat = np.round(float(lat), decimals = 3)

    lon = np.round(float(lon), decimals = 3)

    elv = np.round(float(elv), decimals = 0)
    
    loc = f'lat: {lat}$^o$, lon: {lon}$^o$, alt: {elv} m'
        
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
        
    title = f'{lidar_1} {channel_1} vs {lidar_2} {channel_2} at {loc}\n'+\
                f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+\
                    f' {zan}'+r'$^{o}$ off-zenith' 
    
    return(title)

def sm_text(smooth, sm_lims, sm_win, sm_expo):
    
    sm_llim = np.round(float(sm_lims[0]), decimals = 3)

    sm_ulim = np.round(float(sm_lims[-1]), decimals = 3)

    if smooth == True:
        if isinstance(sm_win,list):
            sm_lwin = np.round(float(sm_win[0]), decimals = 0)
            sm_uwin = np.round(float(sm_win[-1]), decimals = 0)
            
            if sm_lwin > sm_uwin:
                change = 'Decrease'
            else:
                change = 'Increase'

            if sm_expo == True:
                sm_type = 'Exp.'                
            else:
                sm_type = 'Lin.'
                
            sm_part = f' - Smoothing: {sm_llim} to {sm_ulim} km, Win.: {sm_lwin}m to {sm_uwin}m, {change}: {sm_type}'

        else:
            sm_win = np.round(float(sm_win), decimals = 0)
            
            sm_part = f' - Smoothing: {sm_llim} to {sm_ulim} km, Win.: {sm_win}m'
   
    else:
        sm_part = ' - No Smoothing'
    
    return(sm_part)

def mol_text(mol_method, st_name, wmo_id, wban_id, rs_start_date, rs_start_time):
    
    rs_date = f'{rs_start_date[6:]}.{rs_start_date[4:6]}.{rs_start_date[:4]}'
    
    rs_start = f'{rs_start_time[:2]}:{rs_start_time[2:4]}'    
    
    if mol_method == 'Radiosonde': 
        mol_part = f' - {mol_method} {st_name} {rs_date} {rs_start}UT {wmo_id} {wban_id} '
    else: mol_part = f' - {mol_method}'
    
    return(mol_part)