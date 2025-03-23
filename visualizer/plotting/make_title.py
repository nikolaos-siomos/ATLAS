#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 22:34:11 2022

@author: nick
"""

import numpy as np

def quicklook(channel, metadata, args):
    
    start_date = metadata['start_date']
    start_time = metadata['start_time'] 
    stop_time = metadata['stop_time'] 
    lidar_name = metadata['lidar_name'] 
    config_id = metadata['config_id']
    config_name = metadata['config_name']
    scc_channel_id = metadata['scc_channel_id'].copy().loc[channel].values
    laser_pointing_angle = metadata['laser_pointing_angle']
    station_name = metadata['station_name']
    smooth = args['smooth']
    sm_lims = args['smoothing_range']
    sm_win = args['smoothing_window']
    sm_expo = args['smooth_exponential']

    laser_pointing_angle = np.round(float(laser_pointing_angle), decimals = 1)
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)

    dateloc_part = dateloc_text(start_date, start_time, stop_time, laser_pointing_angle)
    
    channel_part = channel_text(lidar_name, station_name, channel, scc_channel_id)

    config_part =  config_text(config_id, config_name)

    title = channel_part + ' - ' + sm_part + '\n'+\
        config_part + ' - ' + dateloc_part
                        
    return(title)

def rayleigh(channel, metadata, args, is_mask = False):
    
    start_date = metadata['start_date']
    start_time = metadata['start_time'] 
    stop_time = metadata['stop_time'] 
    lidar_name = metadata['lidar_name'] 
    config_id = metadata['config_id']
    config_name = metadata['config_name']
    scc_channel_id = metadata['scc_channel_id'].copy().loc[channel].values
    laser_pointing_angle = metadata['laser_pointing_angle']
    station_name = metadata['station_name']
    dwl = metadata['dwl'].copy().loc[channel].values
    ewl = metadata['ewl'].copy().loc[channel].values
    bdw = metadata['bdw'].copy().loc[channel].values
    smooth = args['smooth']
    sm_lims = args['smoothing_range']
    sm_win = args['smoothing_window']
    sm_expo = args['smooth_exponential']
    mol_method = metadata['mol_method']
    rs_station_name = metadata['rs_station_name']
    rs_start_date = metadata['rs_start_date']
    rs_start_time = metadata['rs_start_time'] 
    wmo_id = metadata['wmo_id']
    wban_id = metadata['wban_id']
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    mol_part = mol_text(mol_method, rs_station_name, wmo_id, wban_id, rs_start_date, rs_start_time)
            
    if_part = if_text(ewl, dwl, bdw)
    
    dateloc_part = dateloc_text(start_date, start_time, stop_time, laser_pointing_angle)

    channel_part = channel_text(lidar_name, station_name, channel, scc_channel_id)
    
    config_part =  config_text(config_id, config_name)
    
    if is_mask == False:
        title = channel_part + ' - ' + dateloc_part + ' - ' + sm_part + '\n'+\
                    config_part + ' - ' + mol_part + ' - ' + if_part 
    else:
        title = channel_part + ' - ' + sm_part + '\n'+\
                    config_part + ' - ' + dateloc_part + '\n'+\
                        mol_part+ ' - ' + if_part

    return(title)

def telecover(channel, metadata, args, iters, is_ring = False):

    start_date = metadata['start_date']
    start_time = metadata['start_time'] 
    stop_time = metadata['stop_time'] 
    lidar_name = metadata['lidar_name'] 
    config_id = metadata['config_id']
    config_name = metadata['config_name']
    scc_channel_id = metadata['scc_channel_id'].copy().loc[channel].values
    laser_pointing_angle = metadata['laser_pointing_angle']
    station_name = metadata['station_name']
    dwl = metadata['dwl'].copy().loc[channel].values
    ewl = metadata['ewl'].copy().loc[channel].values
    bdw = metadata['bdw'].copy().loc[channel].values
    smooth = args['smooth']
    sm_lims = args['smoothing_range']
    sm_win = args['smoothing_window']
    sm_expo = args['smooth_exponential']
    
    if is_ring == False:
        sampling_time_per_sector = metadata['sampling_time_per_quadrant']
    else:
        sampling_time_per_sector = metadata['sampling_time_per_ring']
        
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)

    dateloc_part = dateloc_text(start_date, start_time, stop_time, laser_pointing_angle)
    
    iter_part = iter_text(iters, sampling_time_per_sector)
    
    channel_part = channel_text(lidar_name, station_name, channel, scc_channel_id)

    config_part =  config_text(config_id, config_name)
    
    if_part = if_text(ewl, dwl, bdw)
        
    title = channel_part + ' - ' + iter_part + ' - ' + sm_part+'\n'+\
        config_part + ' - ' + dateloc_part + ' - ' + if_part
                
    return(title)

def polarization_calibration(channel_r, channel_t, metadata, args):
    
    start_date_cal = metadata['start_date_cal']
    start_time_cal = metadata['start_time_cal'] 
    stop_time_cal = metadata['stop_time_cal'] 
    start_date_ray = metadata['start_date_ray']
    start_time_ray = metadata['start_time_ray'] 
    stop_time_ray = metadata['stop_time_ray'] 
    lidar_name = metadata['lidar_name'] 
    config_id = metadata['config_id']
    config_name = metadata['config_name']
    scc_channel_id_r = metadata['scc_channel_id'].copy().loc[channel_r].values
    scc_channel_id_t = metadata['scc_channel_id'].copy().loc[channel_t].values
    laser_pointing_angle_cal = metadata['laser_pointing_angle_cal']
    laser_pointing_angle_ray = metadata['laser_pointing_angle_ray']
    station_name = metadata['station_name']
    dwl_r = metadata['dwl_cal'].copy().loc[channel_r].values
    ewl_r = metadata['ewl_cal'].copy().loc[channel_r].values
    bdw_r = metadata['bdw_cal'].copy().loc[channel_r].values
    dwl_t = metadata['dwl_cal'].copy().loc[channel_t].values
    ewl_t = metadata['ewl_cal'].copy().loc[channel_t].values
    bdw_t = metadata['bdw_cal'].copy().loc[channel_t].values
    smooth = args['smooth']
    sm_lims = args['smoothing_range']
    sm_win = args['smoothing_window']
    sm_expo = args['smooth_exponential']
    mol_method = metadata['mol_method']
    rs_station_name = metadata['rs_station_name']
    rs_start_date = metadata['rs_start_date']
    rs_start_time = metadata['rs_start_time'] 
    wmo_id = metadata['wmo_id']
    wban_id = metadata['wban_id']
    
    if_part_r = if_text(ewl_r, dwl_r, bdw_r, label = 'Channel R')

    if_part_t = if_text(ewl_t, dwl_t, bdw_t, label = 'Channel T')

    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    mol_part = mol_text(mol_method, rs_station_name, wmo_id, wban_id, rs_start_date, rs_start_time)
        
    channel_part = channel_text_ratio(lidar_name, station_name, channel_r, channel_t, 
                                      scc_channel_id_r, scc_channel_id_t)
    
    config_part = config_text(config_id, config_name)

    dateloc_part_cal = dateloc_text(start_date_cal, start_time_cal, stop_time_cal, laser_pointing_angle_cal)
   
    dateloc_part_ray = dateloc_text(start_date_ray, start_time_ray, stop_time_ray, laser_pointing_angle_ray)
    
    title = channel_part + ' - ' + sm_part + '\n'+\
            'Calibration: ' + dateloc_part_cal +' - '+\
            'Rayleigh: ' + dateloc_part_ray + '\n' +\
             config_part + ' - ' + mol_part + ' - ' + if_part_r + ' - ' + if_part_t

    return(title)

def intercomparison(channel_1, channel_2, metadata_1, metadata_2, args):
    
    station_name_1 = metadata_1['station_name']
    start_date_1 = metadata_1['start_date']
    start_time_1 = metadata_1['start_time'] 
    stop_time_1 = metadata_1['stop_time'] 
    lidar_name_1 = metadata_1['lidar_name'] 
    laser_pointing_angle_1 = metadata_1['laser_pointing_angle']

    station_name_2 = metadata_2['station_name']
    start_date_2 = metadata_2['start_date']
    start_time_2 = metadata_2['start_time'] 
    stop_time_2 = metadata_2['stop_time']     
    lidar_name_2 = metadata_2['lidar_name']
    laser_pointing_angle_2 = metadata_2['laser_pointing_angle']

    dwl_1 = metadata_1['dwl'].copy().loc[channel_1].values
    ewl_1 = metadata_1['ewl'].copy().loc[channel_1].values
    bdw_1 = metadata_1['bdw'].copy().loc[channel_1].values
    dwl_2 = metadata_2['dwl'].copy().loc[channel_2].values
    ewl_2 = metadata_2['ewl'].copy().loc[channel_2].values
    bdw_2 = metadata_2['bdw'].copy().loc[channel_2].values
    
    smooth = args['smooth']
    sm_lims = args['smoothing_range']
    sm_win = args['smoothing_window']
    sm_expo = args['smooth_exponential']
    
    mol_method_1 = metadata_1['mol_method']
    rs_station_name_1 = metadata_1['rs_station_name']
    rs_start_date_1 = metadata_1['rs_start_date']
    rs_start_time_1 = metadata_1['rs_start_time'] 
    wmo_id_1 = metadata_1['wmo_id']
    wban_id_1 = metadata_1['wban_id']
    
    mol_method_2 = metadata_2['mol_method']
    rs_station_name_2 = metadata_2['rs_station_name']
    rs_start_date_2 = metadata_2['rs_start_date']
    rs_start_time_2 = metadata_2['rs_start_time'] 
    wmo_id_2 = metadata_2['wmo_id']
    wban_id_2 = metadata_2['wban_id']
    
    sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)
        
    mol_part_1 = mol_text(mol_method_1, rs_station_name_1, wmo_id_1, wban_id_1, rs_start_date_1, rs_start_time_1)
    mol_part_2 = mol_text(mol_method_2, rs_station_name_2, wmo_id_2, wban_id_2, rs_start_date_2, rs_start_time_2)
            
    if_part_1 = if_text(ewl_1, dwl_1, bdw_1, label = channel_1)
    if_part_2 = if_text(ewl_2, dwl_2, bdw_2, label = channel_2)
    
    dateloc_part_1 = dateloc_text(start_date_1, start_time_1, stop_time_1, laser_pointing_angle_1)
    dateloc_part_2 = dateloc_text(start_date_2, start_time_2, stop_time_2, laser_pointing_angle_2)

    channel_part_1 = channel_text(lidar_name_1, station_name_1, channel_1, '')
    channel_part_2 = channel_text(lidar_name_2, station_name_2, channel_2, '')
        
    title = channel_part_1 + ' - ' + dateloc_part_1 + ' - ' + mol_part_1 + ' - ' + if_part_1 + '\n'+\
                channel_part_2 + ' - ' + dateloc_part_2 + ' - ' + mol_part_2 + ' - ' + if_part_2 + '\n'+\
                    sm_part
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

def channel_text(lidar_name, station, channel, scc_channel_id = ''):
    
    if scc_channel_id != '':
        channel_part = f'{lidar_name} {station} {channel} ({scc_channel_id})'.strip()
    else:
        channel_part = f'{lidar_name} {station} {channel}'.strip()

    return(channel_part)

def channel_text_ratio(lidar_name, station, channel_r, channel_t, scc_channel_id_r, scc_channel_id_t):
    
    channel_part = f'{lidar_name} {station} {channel_r} ({scc_channel_id_r}) to {channel_t} ({scc_channel_id_t})'.strip()

    return(channel_part)

def config_text(config_id, config_name):
    
    config_part = f'Config {config_id}: {config_name}'.strip()

    return(config_part)

def if_text(ewl, dwl, bdw, label = ''):
    
    ewl = np.round(float(ewl), decimals = 2)

    dwl = np.round(float(dwl), decimals = 2)
    
    bdw = np.round(float(bdw), decimals = 2)
    
    if len(label) > 0:
        
        if_part = f'EWL: {ewl} nm, DWL: {dwl} nm, BDW: {bdw} nm'.strip()
    else:
        if_part = f'{label} EWL: {ewl} nm, DWL: {dwl} nm, BDW: {bdw} nm'.strip()

    return(if_part)

def dateloc_text(start_date, start_time, stop_time, laser_pointing_angle):
    
    laser_pointing_angle = np.round(float(laser_pointing_angle), decimals = 1)

    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{stop_time[:2]}:{stop_time[2:4]}:{stop_time[4:6]}'
    
    dateloc_part = f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+f'{laser_pointing_angle}'+r'$^{o}$ ZA'.strip()
    
    return(dateloc_part)

def iter_text(iters, sampling_time_per_sector):
    
    iter_part = f'Iterations: {iters}, Sampling Time per Sector: {sampling_time_per_sector} s'.strip()
    
    return(iter_part)

def mol_text(mol_method, rs_station_name, wmo_id, wban_id, rs_start_date, rs_start_time):
    
    rs_date = f'{rs_start_date[6:]}.{rs_start_date[4:6]}.{rs_start_date[:4]}'
    
    rs_start = f'{rs_start_time[:2]}:{rs_start_time[2:4]}'    
    
    if mol_method == 'Radiosonde': 
        mol_part = f'{mol_method} {rs_station_name} {rs_date} {rs_start}UT {wmo_id} {wban_id}'.strip()
    else: mol_part = f'{mol_method}'
    
    return(mol_part)