#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 22:34:11 2022

@author: nick
"""

import numpy as np

def quicklook(start_date, start_time, end_time, lidar, channel, 
              station_id, lidar_id, version_id, config_id, config_name, scc_id, 

              zan, loc, smooth, sm_lims, sm_win, sm_expo):
    
    zan = np.round(float(zan), decimals = 1)
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)

    dateloc_part = dateloc_text(start_date, start_time, end_time, zan)
    
    channel_part = channel_text(lidar, loc, channel, scc_id)

    config_part =  config_text(config_id, config_name)

    title = channel_part + ' - ' + sm_part + '\n'+\
        config_part + ' - ' + dateloc_part
                        
    return(title)

def rayleigh(start_date, start_time, end_time, lidar, channel,
             station_id, lidar_id, version_id, config_id, config_name, scc_id, 
             zan, loc, ewl, dwl, bdw, smooth, sm_lims, sm_win, sm_expo,
             mol_method, st_name, rs_start_date, rs_start_time, wmo_id, wban_id):
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    mol_part = mol_text(mol_method, st_name, wmo_id, wban_id, rs_start_date, rs_start_time)
            
    if_part = if_text(ewl, dwl, bdw)
    
    dateloc_part = dateloc_text(start_date, start_time, end_time, zan)

    channel_part = channel_text(lidar, loc, channel, scc_id)
    
    config_part =  config_text(config_id, config_name)
    
    title = channel_part + ' - ' + dateloc_part + ' - ' + sm_part + '\n'+\
                config_part + ' - ' + mol_part + ' - ' + if_part 
    
    return(title)

def rayleigh_mask(start_date, start_time, end_time, lidar, channel,
                  station_id, lidar_id, version_id, config_id, config_name, scc_id, 
                  zan, loc, ewl, dwl, bdw, smooth, sm_lims, sm_win, sm_expo,
                  mol_method, st_name, rs_start_date, rs_start_time, wmo_id, wban_id):
    
    zan = np.round(float(zan), decimals = 1)
    
    ewl = np.round(float(ewl), decimals = 2)

    dwl = np.round(float(dwl), decimals = 2)
    
    bdw = np.round(float(bdw), decimals = 2)
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
    
    if_part = if_text(ewl, dwl, bdw)
        
    mol_part = mol_text(mol_method, st_name, wmo_id, wban_id, rs_start_date, rs_start_time)
    
    dateloc_part = dateloc_text(start_date, start_time, end_time, zan)

    channel_part = channel_text(lidar, loc, channel, scc_id)
    
    config_part =  config_text(config_id, config_name)

    title = channel_part + ' - ' + sm_part + '\n'+\
                config_part + ' - ' + dateloc_part + '\n'+\
                    mol_part+ ' - ' + if_part
    
    return(title)


def telecover(start_date, start_time, end_time, lidar, channel, 
              station_id, lidar_id, version_id, config_id, config_name, scc_id, 
              zan, loc, iters, sampling, smooth, sm_lims, sm_win, sm_expo):
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)

    dateloc_part = dateloc_text(start_date, start_time, end_time, zan)
    
    iter_part = iter_text(iters, sampling)
    
    channel_part = channel_text(lidar, loc, channel, scc_id)

    config_part =  config_text(config_id, config_name)
        
    title = channel_part + ' - ' + iter_part + ' - ' + sm_part+'\n'+\
        config_part + ' - ' + dateloc_part
                
    return(title)

def polarization_calibration(start_date_cal, start_time_cal, end_time_cal,
                             start_date_ray, start_time_ray, end_time_ray,  
                             lidar, channel_r, channel_t, 
                             station_id, lidar_id, version_id, config_id, 
                             config_name, scc_id_r, scc_id_t, 
                             zan_cal, zan_ray, loc, 
                             ewl, dwl, bdw, smooth, sm_lims, sm_win, sm_expo,
                             mol_method, st_name, rs_start_date, rs_start_time,
                             wmo_id, wban_id):
    
    if_part = if_text(ewl, dwl, bdw)

    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    mol_part = mol_text(mol_method, st_name, wmo_id, wban_id, rs_start_date, rs_start_time)
        
    channel_part = channel_text_ratio(lidar, loc, channel_r, channel_t, 
                                      scc_id_r, scc_id_t)
    
    config_part = config_text(config_id, config_name)

    dateloc_part_cal = dateloc_text(start_date_cal, start_time_cal, end_time_cal, zan_cal)
   
    dateloc_part_ray = dateloc_text(start_date_ray, start_time_ray, end_time_ray, zan_ray)
    
    title = channel_part + ' - ' + sm_part + '\n'+\
            'Calibration: ' + dateloc_part_cal +' - '+\
            'Rayleigh: ' + dateloc_part_ray + '\n' +\
             config_part + ' - ' + mol_part + ' - ' + if_part

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
                    f' {zan}'+r'$^{o}$ ZA' 
    
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
                sm_type = 'Exp'                
            else:
                sm_type = 'Lin'
                
            sm_part = f'Smoothing: {sm_llim} to {sm_ulim} km, Win: {sm_lwin}m to {sm_uwin}m, {change}: {sm_type}'

        else:
            sm_win = np.round(float(sm_win), decimals = 0)
            
            sm_part = f'Smoothing: {sm_llim} to {sm_ulim} km, Win: {sm_win}m'
   
    else:
        sm_part = 'No Smoothing'
    
    return(sm_part)

def channel_text(lidar, loc, channel, scc_id):
    
    channel_part = f'{lidar} {loc} {channel} ({scc_id})'.strip()

    return(channel_part)

def channel_text_ratio(lidar, loc, channel_r, channel_t, scc_id_r, scc_id_t):
    
    channel_part = f'{lidar} {loc} {channel_r} ({scc_id_r}) to {channel_t} ({scc_id_t})'.strip()

    return(channel_part)

def config_text(config_id, config_name):
    
    config_part = f'Config {config_id}: {config_name}'.strip()

    return(config_part)

def if_text(ewl, dwl, bdw):
    
    ewl = np.round(float(ewl), decimals = 2)

    dwl = np.round(float(dwl), decimals = 2)
    
    bdw = np.round(float(bdw), decimals = 2)
    
    if_part = f'Emitted WL: {ewl}nm, Received WL: {dwl}nm, Bandwidth: {bdw}nm'.strip()
    
    return(if_part)

def dateloc_text(start_date, start_time, end_time, zan):
    
    zan = np.round(float(zan), decimals = 1)

    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
    
    dateloc_part = f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+f'{zan}'+r'$^{o}$ ZA'.strip()
    
    return(dateloc_part)

def iter_text(iters, sampling):
    
    iter_part = f'Iterations: {iters}, Sampling Time per Sector: {sampling}s'.strip()
    
    return(iter_part)

def mol_text(mol_method, st_name, wmo_id, wban_id, rs_start_date, rs_start_time):
    
    rs_date = f'{rs_start_date[6:]}.{rs_start_date[4:6]}.{rs_start_date[:4]}'
    
    rs_start = f'{rs_start_time[:2]}:{rs_start_time[2:4]}'    
    
    if mol_method == 'Radiosonde': 
        mol_part = f'{mol_method} {st_name} {rs_date} {rs_start}UT {wmo_id} {wban_id}'.strip()
    else: mol_part = f'{mol_method}'
    
    return(mol_part)