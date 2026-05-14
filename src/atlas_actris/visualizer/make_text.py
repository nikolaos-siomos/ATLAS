#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 23:23:58 2026

@author: nikos
"""

import numpy as np
import xarray as xr
import pandas as pd
from typing import Any, Dict
from version import __version__
from dataclasses import dataclass

@dataclass(frozen=True)
class Libraries:
    system_info: xr.DataArray
    channel_info: xr.DataArray
    time_info: xr.DataArray
    settings: Dict[str, Any]
    
class GenerateText:
    
    def __init__(self, lib: Libraries, atlas_channel_id):
        
        self.system_info = lib.system_info
        self.channel_info = lib.channel_info
        self.time_info = lib.time_info
        
        self.settings = lib.settings
        
        start_timestamp = lib.time_info.isel({'time': 0}).sel({'parameters': 'start_time'}).item()
        stop_timestamp = lib.time_info.isel({'time': -1}).sel({'parameters': 'end_time'}).item()
        
        self.start_datetime = pd.to_datetime(start_timestamp)
        self.stop_datetime = pd.to_datetime(stop_timestamp)
        
        self.atlas_channel_id = atlas_channel_id
        
        self.scc_channel_ids = lib.channel_info.sel({'parameters':'scc_channel_id'})
    
        self.station_id = lib.system_info.loc["station_id"].item()
        self.lidar_id = lib.system_info.loc["lidar_id"].item()
        self.version_id = lib.system_info.loc["version_id"].item()
        self.configuration_id = lib.system_info.loc["configuration_id"].item()
        
    def make_filename(self, qa_test, extra_type = '', extra_channel = None):
        
        scc_channel_id = self.scc_channel_ids.sel({'channel':self.atlas_channel_id}).item()
        
        common_parts = [
            self.station_id, 
            self.configuration_id, 
            self.start_datetime.strftime("%Y%m%d"), 
            self.start_datetime.strftime("%H%M%S"), 
            qa_test, 
            self.atlas_channel_id, 
            scc_channel_id
            ]
        
        final_parts = [
            extra_type, 
            'ATLAS', 
            __version__
            ]
        
        if extra_channel is not None:
            extra_scc_channel = self.scc_channel_ids.sel({'channel':extra_channel}).item()

            extra_parts = [
                extra_channel,
                extra_scc_channel
                ]
                
        else:
            extra_parts = []
        
        parts = common_parts + extra_parts + final_parts
        filename = "_".join(
            str(part) for part in parts
            if part is not None and str(part) != ""
            )   
        
        return filename

    def make_quicklook_title(self):
            
        scc_channel_id = self.scc_channel_ids.sel({'channel':self.atlas_channel_id}).item()

        start_date = self.start_datetime.strftime("%Y%m%d")
        start_time = self.start_datetime.strftime("%H%M%S")
        stop_time = self.stop_datetime.strftime("%H%M%S")
        
        station_name = self.system_info.loc['station_name'].item()
        config_id = self.system_info.loc['configuration_id'].item()
        config_name = self.system_info.loc['configuration_name'].item()
        lidar_name = self.system_info.loc['lidar_name'].item()
        laser_pointing_angle = self.system_info.loc['zenith_angle'].item()
        
        smooth = self.settings['smooth']
        sm_lims = self.settings['smoothing_range']
        sm_win = self.settings['smoothing_constant_window']
        sm_expo = self.settings['smoothing_exponential']

        laser_pointing_angle = np.round(float(laser_pointing_angle), decimals = 1)
        
        sm_part = sm_text(smooth, sm_lims, sm_win, sm_expo)

        dateloc_part = dateloc_text(start_date, start_time, stop_time, laser_pointing_angle)
        
        channel_part = channel_text(lidar_name, station_name, self.atlas_channel_id, scc_channel_id)

        config_part =  config_text(config_id, config_name)

        title = channel_part + ' - ' + sm_part + '\n'+\
            config_part + ' - ' + dateloc_part
                            
        return title 
    

def sm_text(smooth, sm_lims, sm_win, sm_expo):

    if smooth != True:
        return "No Smoothing"

    sm_llim = np.round(float(sm_lims[0]), decimals=3)
    sm_ulim = np.round(float(sm_lims[-1]), decimals=3)

    if isinstance(sm_win, (list, tuple, np.ndarray)):
        sm_lwin = np.round(float(sm_win[0]), decimals=0)
        sm_uwin = np.round(float(sm_win[-1]), decimals=0)

        if sm_lwin > sm_uwin:
            change = "Decrease"
        else:
            change = "Increase"

        if sm_expo == True:
            sm_type = "Exp"
        else:
            sm_type = "Lin"

        sm_part = (
            f"Smoothing: {sm_llim} to {sm_ulim} km, "
            f"Win: {sm_lwin}m to {sm_uwin}m, "
            f"{change}: {sm_type}"
        )

    else:
        sm_win = np.round(float(sm_win), decimals=0)
        sm_part = f"Smoothing: {sm_llim} to {sm_ulim} km, Win: {sm_win}m"

    return sm_part

def channel_text(lidar_name, station, channel, scc_channel_id = ''):
    
    if scc_channel_id != '':
        channel_part = f'{lidar_name} {station} {channel} ({scc_channel_id})'.strip()
    else:
        channel_part = f'{lidar_name} {station} {channel}'.strip()

    return channel_part

def channel_text_ratio(lidar_name, station, channel_r, channel_t, scc_channel_id_r, scc_channel_id_t):
    
    channel_part = f'{lidar_name} {station} {channel_r} ({scc_channel_id_r}) to {channel_t} ({scc_channel_id_t})'.strip()

    return channel_part

def config_text(config_id, config_name):
    
    config_part = f'Config {config_id}: {config_name}'.strip()

    return(config_part)

def if_text(ewl, dwl, bdw, label=""):

    ewl = np.round(float(ewl), decimals=2)
    dwl = np.round(float(dwl), decimals=2)
    bdw = np.round(float(bdw), decimals=2)

    if len(label) > 0:
        if_part = f"{label} EWL: {ewl} nm, DWL: {dwl} nm, BDW: {bdw} nm"
    else:
        if_part = f"EWL: {ewl} nm, DWL: {dwl} nm, BDW: {bdw} nm"

    return if_part.strip()

def dateloc_text(start_date, start_time, stop_time, laser_pointing_angle):
    
    laser_pointing_angle = np.round(float(laser_pointing_angle), decimals = 1)

    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{stop_time[:2]}:{stop_time[2:4]}:{stop_time[4:6]}'
    
    dateloc_part = (
        f"On {date} from {start} to {end} UTC, "
        + r"$\nearrow$"
        + f"{laser_pointing_angle}"
        + r"$^{o}$ ZA"
    )  
    
    return dateloc_part

def iter_text(iters, sampling_time_per_sector):
    
    iter_part = f'Iterations: {iters}, Sampling Time per Sector: {sampling_time_per_sector} s'.strip()
    
    return(iter_part)

def mol_text(mol_method, rs_station_name, wmo_id, wban_id, rs_start_date, rs_start_time):
    
    rs_date = f'{rs_start_date[6:]}.{rs_start_date[4:6]}.{rs_start_date[:4]}'
    
    rs_start = f'{rs_start_time[:2]}:{rs_start_time[2:4]}'    
    
    if mol_method == 'Radiosonde': 
        mol_part = f'{mol_method} {rs_station_name} {rs_date} {rs_start}UT {wmo_id} {wban_id}'.strip()
    else: mol_part = f'{mol_method}'
    
    return mol_part