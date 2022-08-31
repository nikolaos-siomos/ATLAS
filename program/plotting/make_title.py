#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 22:34:11 2022

@author: nick
"""

import numpy as np

def quicklook(start_time, end_time, lidar, zan, lat, lon, elv):
    
    loc = f'lat: {lat}$^o$, lon: {lon}$^o$, alt: {elv} m'
        
    date = np.datetime64(start_time,'us').item().strftime('%d.%m.%Y')
    
    start = np.datetime64(start_time,'us').item().strftime('%H:%M')
    
    end = np.datetime64(end_time,'us').item().strftime('%H:%M')
    
    title = f'{lidar} Time-Height cross sections at {loc}\n'+\
                f'On {date} from {start} to {end} UTC, '+r'$\nearrow$'+\
                    f' {zan}'+r'$^{o}$ off-zenith' 
    
    return(title)