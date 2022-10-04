#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 22:34:11 2022

@author: nick
"""

import numpy as np

def quicklook(start_time, end_time, lidar, channel, zan, lat, lon, elv):
    
    zan = np.round(float(zan), decimals = 1)
    
    lat = np.round(float(lat), decimals = 3)

    lon = np.round(float(lon), decimals = 3)

    elv = np.round(float(elv), decimals = 0)
    
    loc = f'lat: {lat}$^o$, lon: {lon}$^o$, alt: {elv} m'
        
    date = np.datetime64(start_time,'us').item().strftime('%d.%m.%Y')
    
    start = np.datetime64(start_time,'us').item().strftime('%H:%M:%S')
    
    end = np.datetime64(end_time,'us').item().strftime('%H:%M:%S')
    
    title = f'Channel {channel} : {lidar} at {loc}\n'+\
                f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+\
                    f' {zan}'+r'$^{o}$ off-zenith' 
    
    return(title)

def rayleigh(start_date, start_time, end_time, lidar, channel, 
             zan, lat, lon, elv):
    
    zan = np.round(float(zan), decimals = 1)
    
    lat = np.round(float(lat), decimals = 3)

    lon = np.round(float(lon), decimals = 3)

    elv = np.round(float(elv), decimals = 0)
    
    loc = f'lat: {lat}$^o$, lon: {lon}$^o$, alt: {elv} m'
        
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
        
    title = f'Channel {channel} : {lidar} at {loc}\n'+\
                f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+\
                    f' {zan}'+r'$^{o}$ off-zenith' 
    
    return(title)

def telecover(start_date, start_time, end_time, lidar, channel, 
             zan, lat, lon, elv):
    
    zan = np.round(float(zan), decimals = 1)
    
    lat = np.round(float(lat), decimals = 3)

    lon = np.round(float(lon), decimals = 3)

    elv = np.round(float(elv), decimals = 0)
    
    loc = f'lat: {lat}$^o$, lon: {lon}$^o$, alt: {elv} m'
        
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
        
    title = f'Channel {channel} : {lidar} at {loc}\n'+\
                f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+\
                    f' {zan}'+r'$^{o}$ off-zenith' 
    
    return(title)

def polarization_calibration(start_date_cal, start_time_cal, end_time_cal,
                             start_date_ray, start_time_ray, end_time_ray,  
                             lidar, channel_r, channel_t, zan, lat, lon, elv):
    
    zan = np.round(float(zan), decimals = 1)
    
    lat = np.round(float(lat), decimals = 3)

    lon = np.round(float(lon), decimals = 3)

    elv = np.round(float(elv), decimals = 0)
    
    loc = f'lat: {lat}$^o$, lon: {lon}$^o$, alt: {elv} m'
        
    date_cal = f'{start_date_cal[6:]}.{start_date_cal[4:6]}.{start_date_cal[:4]}'
    
    start_cal = f'{start_time_cal[:2]}:{start_time_cal[2:4]}:{start_time_cal[4:6]}'

    end_cal = f'{end_time_cal[:2]}:{end_time_cal[2:4]}:{end_time_cal[4:6]}'

    date_ray = f'{start_date_ray[6:]}.{start_date_ray[4:6]}.{start_date_ray[:4]}'
    
    start_ray = f'{start_time_ray[:2]}:{start_time_ray[2:4]}:{start_time_ray[4:6]}'

    end_ray = f'{end_time_ray[:2]}:{end_time_ray[2:4]}:{end_time_ray[4:6]}'
               
    title = f'Ratio {channel_r} to {channel_t}: {lidar} at {loc}\n '+\
                f'Calibration on {date_cal} from {start_cal} to {end_cal} UTC, '+\
                    r'$\nearrow$' + f' {zan}' + r'$^{o}$ off-zenith'+ f'\n'+\
                        f'Rayleigh on {date_ray} from {start_ray} to {end_ray} UTC, '+\
                            r'$\nearrow$' + f' {zan}' + r'$^{o}$ off-zenith'
    
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