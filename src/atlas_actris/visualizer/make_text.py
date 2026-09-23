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

label = {
    'drk' : "Dark",
    'ray' : "Rayleigh",
    'tlc' : "Telecover",
    'pcb' : "Polarization Calibration",
    }

telescope_map = {
    'n' : 'Near Range 1',
    'm' : 'Near Range 2',
    'l' : 'Near Range 3',
    'f' : 'Far Range 1',
    'g' : 'Far Range 2',
    'h' : 'Far Range 3',
    'x' : 'Telescope 1',
    'y' : 'Telescope 2',
    'z' : 'Telescope 3',
    }

mode_map = {
    'a' : 'analog',
    'p' : 'photon'
    }

type_map = {
    'p' : 'Co-polar',
    'c' : 'Cross-polar',
    't' : 'Total',
    'v' : 'Vibrational Raman',
    'r' : 'Rotational Raman',
    'a' : 'Cabannes',
    'f' : 'Fluorescence'
    }

subtype_map = {
    'r' : 'Reflected',
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
    'x' : ''
    }

qa_dataset_labels = {
    'drk': 'Long dark',
    'ray': 'Rayleigh',
    'ray_pcb': 'Rayleigh pol. cal. mode',
    'tlc': 'Quadrant telecover',
    'tlc_rin': 'Ring telecover',
    'pcb': 'Pol. cal.',
    'pcb_aux': 'Auxialliary pol. cal.',
    'trg': 'Zero bin',
    'dtm': 'Dead time',
    "pcb_p45": 'Pol. cal.',
    "pcb_m45": 'Pol. cal.',
    "tlc_north": 'Quadr. telecover',
    "tlc_east": 'Quadr. telecover',
    "tlc_south": 'Quadr. telecover',
    "tlc_west": 'Quadrant telecover',
    "tlc_inner": 'Ring telecover',
    "tlc_outer": 'Ring telecover',
    "pcb_aux_p45": 'Auxialliary pol. cal.',
    "pcb_aux_m45": 'Auxialliary pol. cal.',
    "drk_ray": 'Rayleigh dark',
    "drk_pcb": 'Pol. cal. dark',
    "drk_tlc": 'Quadr. telecover dark',
    "drk_tlc_rin": 'Ring telecover dark',
    "drk_trg": 'Zero bin dark',
    "drk_dtm": 'Dead time dark measurent',
    "drk_ray_pcb": 'Rayleigh pol. cal. mode dark',
    "drk_pcb_aux": 'Auxilliary pol. cal. dark',
    }

qa_dataset_label_extensions = {
    'drk': 'test',
    'ray': 'measurement',
    'ray_pcb': 'measurement',
    'tlc': 'test',
    'tlc_rin': 'test',
    'pcb': 'measurement',
    'pcb_aux': 'measurement',
    'trg': 'test',
    'dtm': 'measurement',
    "pcb_p45": 'measurement | +45°',
    "pcb_m45": 'measurement | -45°',
    "tlc_north": 'test | North sector',
    "tlc_east": 'test | East sector',
    "tlc_south": 'test | South sector',
    "tlc_west": 'test | West sector',
    "tlc_inner": 'test (Inner sector)',
    "tlc_outer": 'test (Outer sector)',
    "pcb_aux_p45": 'measurement | -45°',
    "pcb_aux_m45": 'measurement | +45°',
    "drk_ray": 'measurent for Rayleigh',
    "drk_pcb": 'measurent for pol. cal.',
    "drk_tlc": 'measurent',
    "drk_tlc_rin": 'measurent',
    "drk_trg": 'measurent',
    "drk_dtm": 'measurent',
    "drk_ray_pcb": 'measurent',
    "drk_pcb_aux": 'measurent',
    }
@dataclass(frozen=True)
class Libraries:
    caller_info: Dict[str, Any]
    metadata: Dict[str, Any]
    extra_metadata: Dict[str, Any]
    settings: Dict[str, Any]
    qa_test_info: Dict[str, Any]
    
class GenerateText:
    
    def __init__(self, lib: Libraries):
        
        self.caller_info = lib.caller_info
        self.metadata = lib.metadata
        self.extra_metadata = lib.extra_metadata
        self.settings = lib.settings
        self.qa_test_info = lib.qa_test_info
        
        self.metadata['start_timestamp'] = pd.Timestamp(self.metadata['start_time_first'])
        self.metadata['stop_timestamp'] = pd.Timestamp(self.metadata['end_time_last'])
                
        if self.extra_metadata:
            self.extra_metadata['start_timestamp'] = pd.Timestamp(self.extra_metadata['start_time_first'])
            self.extra_metadata['stop_timestamp'] = pd.Timestamp(self.extra_metadata['end_time_last'])
        
            self.dateloc_part_extra = dateloc_text(
                start_timestamp = self.extra_metadata['start_timestamp'], 
                stop_timestamp = self.extra_metadata['stop_timestamp'], 
                laser_pointing_angle = self.extra_metadata['zenith_angle']
                )
            
        self.dateloc_part = dateloc_text(
            start_timestamp = self.metadata['start_timestamp'], 
            stop_timestamp = self.metadata['stop_timestamp'], 
            laser_pointing_angle = self.metadata['zenith_angle']
            )

        self.system_part = system_text(
            lidar_name = self.metadata['lidar_name'], 
            station_name = self.metadata['station_name'], 
            )   
        
        self.channel_part = channel_text(
            atlas_channel_id = self.metadata['atlas_channel_id'], 
            scc_channel_id = self.metadata['scc_channel_id']
            )

        self.config_part =  config_text(
            config_id = self.metadata['configuration_id'], 
            config_name = self.metadata['configuration_name']
            )
        
    def make_filename(self, qa_test, extra_type = '', extra_metadata = {}):
               
        start_date = self.metadata['start_timestamp'].strftime("%Y%m%d")
        start_time = self.metadata['start_timestamp'].strftime("%H%M%S")
        
        common_parts = [
            self.metadata['station_id'], 
            self.metadata['configuration_id'], 
            start_date, 
            start_time, 
            qa_test, 
            self.metadata['atlas_channel_id'], 
            self.metadata['scc_channel_id']
            ]
        
        final_parts = [
            extra_type, 
            'ATLAS', 
            __version__
            ]
        
        if extra_metadata:
            extra_channel = extra_metadata['atlas_channel_id']
            extra_scc_channel = extra_metadata['scc_channel_id']

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
    
    def make_filename_viewer(self, qa_test, extra_type = '', extra_metadata = {}):
               
        start_date = self.metadata['start_timestamp'].strftime("%Y%m%d")
        start_time = self.metadata['start_timestamp'].strftime("%H%M%S")
        
        common_parts = [
            self.metadata['station_id'], 
            self.metadata['configuration_id'], 
            start_date, 
            start_time, 
            qa_test, 
            self.metadata['atlas_channel_id'], 
            self.metadata['scc_channel_id']
            ]
        
        stage_parts = [
            self.qa_test_info['stage'],
            self.qa_test_info['db']
            ]
        
        final_parts = [
            extra_type, 
            'ATLAS', 
            __version__
            ]
        
        if extra_metadata:
            extra_channel = extra_metadata['atlas_channel_id']
            extra_scc_channel = extra_metadata['scc_channel_id']

            extra_parts = [
                extra_channel,
                extra_scc_channel
                ]
                
        else:
            extra_parts = []
        
        parts = common_parts + extra_parts + stage_parts + final_parts
        
        filename = "_".join(
            str(part) for part in parts
            if part is not None and str(part) != ""
            )   
        
        return filename
    
    def make_filename_pair(self, qa_test, extra_type = '', extra_metadata = {}):
               
        start_date = self.metadata['start_timestamp'].strftime("%Y%m%d")
        start_time = self.metadata['start_timestamp'].strftime("%H%M%S")
        
        common_parts = [
            self.metadata['station_id'], 
            self.metadata['configuration_id'], 
            start_date, 
            start_time, 
            qa_test, 
            self.qa_test_info['pair'], 
            ]
        
        final_parts = [
            extra_type, 
            'ATLAS', 
            __version__
            ]
        
        if extra_metadata:
            extra_channel = extra_metadata['atlas_channel_id']
            extra_scc_channel = extra_metadata['scc_channel_id']

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

        sm_part = sm_text(
            smooth = self.settings['smooth'], 
            sm_lims = self.settings['smoothing_range'], 
            sm_win = self.settings['smoothing_window'], 
            sm_expo = self.settings['smoothing_exponential']
            )
        
        qa_test = self.qa_test_info['qa_test']
        qck_text = f"Quicklook {qa_dataset_labels[qa_test]}"
        
        title = self.system_part + ' ' + self.channel_part + '\n'+\
            self.config_part  + ' - ' + sm_part + '\n'+\
                qck_text + ' - ' + self.dateloc_part
                            
        return title 
    
    def make_background_title(self):
        
        qa_test = self.qa_test_info['qa_test']
        bgd_text = f"Background {qa_dataset_labels[qa_test]}"

        title = self.system_part + ' ' + self.channel_part + '\n'+\
            self.config_part + '\n'+\
                bgd_text + ' - ' + self.dateloc_part
                            
        return title 
    
    def make_dark_title(self):

        sm_part = sm_text(
            smooth = self.settings['smooth'], 
            sm_lims = self.settings['smoothing_range'], 
            sm_win = self.settings['smoothing_window'], 
            sm_expo = False,
            )

        avg_rate_alias = self.settings['averaging_period']
        averaging_period = self.caller_info[f'{avg_rate_alias}_averaging_period']
        avg_rate_part = f'Averaging period ({avg_rate_alias}): {averaging_period}'
        
        qa_test = self.qa_test_info['qa_test']
        drk_text = f"{qa_dataset_labels[qa_test]}"
        
        title = self.system_part + ' ' + self.channel_part + '\n'+\
            self.config_part  + ' - ' + sm_part + ' - ' + avg_rate_part + '\n'+\
                drk_text + ' - ' + self.dateloc_part
                            
        return title 
    
    
    def make_vldr_title(self):

        sm_part = sm_text(
            smooth = self.settings['smooth'], 
            sm_lims = self.settings['smoothing_range'], 
            sm_win = self.settings['smoothing_window'], 
            sm_expo = self.settings['smoothing_exponential']
            )
        
        pair_part = f"VLDR ID: {self.qa_test_info['pair'].upper()}"
        channel_part = f"Ch R: {self.qa_test_info['ch_r']} - Ch T: {self.qa_test_info['ch_t']}"
        
        title = self.system_part + ' ' + pair_part + ' - ' + channel_part + '\n'+\
            self.config_part  + ' - ' + sm_part + '\n'+\
                'Quicklook VLDR - ' + self.dateloc_part
                            
        return title 
    
    
    def make_rayleigh_fit_title(self, ray_type = 'Rayleigh Fit'):
               
        radiosonde_timestamp = pd.Timestamp(self.metadata['radiosonde_time'])
         
        mol_part = mol_text(
            rs_format = self.metadata['radiosonde_format'], 
            rs_station_name = self.caller_info['rsonde_station_name'], 
            cloudnet_station_name = self.caller_info['cloudnet_station_name'], 
            wmo_id = self.caller_info['rsonde_station_wmo_id'], 
            rs_start_timestamp = radiosonde_timestamp
            )
        
        sm_part = sm_text(
            smooth = self.settings['smooth'], 
            sm_lims = self.settings['smoothing_range'], 
            sm_win = self.settings['smoothing_window'], 
            sm_expo = False
            )
        
        if_part = if_text(
            self.metadata['emitted_wavelength'], 
            self.metadata['detected_wavelength'], 
            self.metadata['channel_bandwidth']
            )
        
        title = self.system_part + ' ' + self.channel_part + ' - ' + self.dateloc_part + ' - ' + sm_part + '\n'+\
                    ray_type + ' - ' + self.config_part + ' - ' + mol_part + ' - ' + if_part 

        return title 


    def make_rayleigh_fit_mask_title(self, ray_type = 'Rayleigh Fit Mask'):
        
        radiosonde_timestamp = pd.Timestamp(self.metadata['radiosonde_time'])
         
        mol_part = mol_text(
            rs_format = self.metadata['radiosonde_format'], 
            rs_station_name = self.caller_info['rsonde_station_name'], 
            cloudnet_station_name = self.caller_info['cloudnet_station_name'], 
            wmo_id = self.caller_info['rsonde_station_wmo_id'], 
            rs_start_timestamp = radiosonde_timestamp
            )
        
        sm_part = sm_text(
            smooth = self.settings['smooth'], 
            sm_lims = self.settings['smoothing_range'], 
            sm_win = self.settings['smoothing_window'], 
            sm_expo = False
            )
        
        if_part = if_text(
            self.metadata['emitted_wavelength'], 
            self.metadata['detected_wavelength'], 
            self.metadata['channel_bandwidth']
            )
        
        title = self.system_part + ' ' + self.channel_part + ' - ' + sm_part + '\n'+\
                    self.config_part + ' - ' + self.dateloc_part + '\n'+\
                        ray_type + ' - ' + mol_part + ' - ' + if_part

        return title 
    
    
    def make_telecover_title(self, tlc_type = 'Qaudrant Telecover'):
        
        sm_part = sm_text_tlc(
            smooth = self.settings['smooth'], 
            sm_win = self.settings['smoothing_window'], 
            nr_ulim = self.settings['near_range_upper_limit']
            )
        
        if_part = if_text(
            self.metadata['emitted_wavelength'], 
            self.metadata['detected_wavelength'], 
            self.metadata['channel_bandwidth']
            )

        iter_part = iter_text(
                iters=self.qa_test_info.get("iters", np.nan),
                sampling_time_per_sector=self.qa_test_info.get("sampling_time_per_sector", np.nan),
            )
        
        title = (
            self.system_part + ' ' + self.channel_part + ' - ' + sm_part + '\n' +\
                iter_part + ' - ' + if_part + '\n' +\
                    tlc_type + ' - ' + self.config_part + ' - ' + self.dateloc_part
        )
                         
        return title 
    
    
    def make_polarization_calibration_title(self, metadata_r, metadata_t, 
                                            pcb_type = 'Pol. Calibration'):
        
        radiosonde_timestamp = pd.Timestamp(self.metadata['radiosonde_time'])
         
        mol_part = mol_text(
            rs_format = self.metadata['radiosonde_format'], 
            rs_station_name = self.caller_info['rsonde_station_name'], 
            cloudnet_station_name = self.caller_info['cloudnet_station_name'], 
            wmo_id = self.caller_info['rsonde_station_wmo_id'], 
            rs_start_timestamp = radiosonde_timestamp
            )
        
        sm_part = sm_text(
            smooth = self.settings['smooth'], 
            sm_lims = self.settings['smoothing_range'], 
            sm_win = self.settings['smoothing_window'], 
            sm_expo = False
            )
        
        if_part_r = if_text(
            metadata_r['emitted_wavelength'], 
            metadata_r['detected_wavelength'], 
            metadata_r['channel_bandwidth']
            )
        
        if_part_t = if_text(
            metadata_t['emitted_wavelength'], 
            metadata_t['detected_wavelength'], 
            metadata_t['channel_bandwidth']
            )
        
        
        ch_r_text = channel_text(
            atlas_channel_id=metadata_r['atlas_channel_id'], 
            scc_channel_id=metadata_r['scc_channel_id']
            )
        
        ch_t_text = channel_text(
            atlas_channel_id=metadata_t['atlas_channel_id'], 
            scc_channel_id=metadata_t['scc_channel_id']
            )
        
        title = self.system_part + ' - ' + sm_part + '\n'+\
            f"Rayleigh: {self.dateloc_part}" + ' - ' + f"Calibration: {self.dateloc_part_extra}" + '\n'+\
                    pcb_type + ' - ' + self.config_part + ' - ' + mol_part + '\n'+\
                        'Ch R: ' + ch_r_text + ', ' + if_part_r + ' - ' +'Ch T: ' + ch_t_text + ', ' + if_part_t

        return title 
    
    
    def make_sig_mlines_title(self):
        
        stage_part = stage_text(self.qa_test_info['stage'])

        if self.qa_test_info['db'] == 'profile':
            sig_type = self.qa_test_info['sig_type'].capitalize()
        else:
            sig_type = f"Averaged {self.qa_test_info['sig_type']}"
            
        qa_dataset_text = qa_dataset_labels[self.qa_test_info['qa_test']]
        alias = self.qa_test_info['qa_test']
        
        signal_part = f"{sig_type} - {qa_dataset_text} ({alias})"
            
            
        
        title = self.system_part + ' ' + self.channel_part + '\n'+\
            self.config_part + ' - ' + stage_part + '\n'+\
                signal_part + ': ' + self.dateloc_part
                            
        return title 
    
    def make_header_dark(self):
                
        parts = []
        
        parts.append(
            header_system_text(
                station_id = self.metadata['station_id'],
                station_name = self.metadata['station_name'],
                lidar_name = self.metadata['lidar_name']
                )
            )

        parts.append(
            header_signal_text(self.metadata['atlas_channel_id'])
            )
        
        parts.append(
            header_time_text(
                start_timestamp = self.metadata['start_timestamp'], 
                stop_timestamp = self.metadata['stop_timestamp'], 
                label = label['drk']
                )
            )
        
        parts.append(
            header_dark_text()
            )
        
        header = '\n'.join(parts)
        
        return header
    
    def make_header_rayleigh_fit(self):
                
        parts = []
        
        radiosonde_timestamp = pd.Timestamp(self.metadata['radiosonde_time'])

        parts.append(
            header_system_text(
                station_id = self.metadata['station_id'],
                station_name = self.metadata['station_name'],
                lidar_name = self.metadata['lidar_name']
                )
            )

        parts.append(
            header_signal_text(self.metadata['atlas_channel_id'])
            )
        
        parts.append(
            header_time_text(
                start_timestamp = self.metadata['start_timestamp'], 
                stop_timestamp = self.metadata['stop_timestamp'], 
                label = label['ray']
                )
            )
        
        parts.append(
            header_radiosonde_text(
                rs_station_name = self.caller_info['rsonde_station_name'], 
                wmo_id = self.caller_info['rsonde_station_wmo_id'], 
                rs_start_timestamp = radiosonde_timestamp, 
                )
            )
        
        parts.append(
            header_rayleigh_fit_text(self.qa_test_info['norm_region'])
            )
        
        header = '\n'.join(parts)
        
        return header
    
    def make_header_telecover(self):
                
        parts = []
        
        parts.append(
            header_system_text(
                station_id = self.metadata['station_id'],
                station_name = self.metadata['station_name'],
                lidar_name = self.metadata['lidar_name']
                )
            )

        parts.append(
            header_signal_text(self.metadata['atlas_channel_id'])
            )
        
        parts.append(
            header_time_text(
                start_timestamp = self.metadata['start_timestamp'], 
                stop_timestamp = self.metadata['stop_timestamp'], 
                label = label['tlc']
                )
            )
        
        parts.append(
            header_telecover_text(
                iters = self.qa_test_info['iters'], 
                secs = self.qa_test_info['available_sectors'], 
                extra_sec = self.qa_test_info['extra_sec']
                )
            )
        
        header = '\n'.join(parts)
        
        return header

    def make_header_polcal(self):
                
        parts = []
        
        parts.append(
            header_system_text(
                station_id = self.metadata['station_id'],
                station_name = self.metadata['station_name'],
                lidar_name = self.metadata['lidar_name']
                )
            )

        parts.append(
            header_signal_text_polcal(
                atlas_channel_id_r = self.metadata['atlas_channel_id'], 
                atlas_channel_id_t = self.extra_metadata['atlas_channel_id'])
            )
        
        parts.append(
            header_time_text(
                start_timestamp = self.metadata['start_timestamp'], 
                stop_timestamp = self.metadata['stop_timestamp'], 
                label = label['pcb']
                )
            )
        
        parts.append(
            header_time_text(
                start_timestamp = self.extra_metadata['start_timestamp'], 
                stop_timestamp = self.extra_metadata['stop_timestamp'], 
                label = label['ray']
                )
            )
        
        parts.append(
            header_polcal_text(
                G_R = self.qa_test_info['G_R'], 
                G_T = self.qa_test_info['G_T'], 
                H_R = self.qa_test_info['G_R'], 
                H_T = self.qa_test_info['G_T'], 
                K = self.qa_test_info['K']
                )
            )
        
        header = '\n'.join(parts)
    
        return header                     
    

def header_system_text(station_id, station_name, lidar_name):
    
    if not station_id:
        station_id = ""
        
    if not station_name:
        station_name = ""
    
    if not lidar_name:
        lidar_name = ""
        
    line_1 = f"station ID = {station_id}"
    
    line_2 = f"system = {lidar_name} - {station_name}"
    
    text = f"{line_1}\n{line_2}"

    return text 

def header_signal_text(atlas_channel_id):
         
    if telescope_map:
        telescope_text = telescope_map[atlas_channel_id[4]]
    else:
        telescope_text = ''
        
    expression = f"{telescope_text} {subtype_map[atlas_channel_id[7]]} {type_map[atlas_channel_id[5]]}"

    text = f"signal = {atlas_channel_id[:4].lstrip('0')}, {expression}, {mode_map[atlas_channel_id[6]]}, dark-subtracted"
    
    return text

def header_signal_text_polcal(atlas_channel_id_r, atlas_channel_id_t):
    
    if telescope_map:
        telescope_text_r = telescope_map[atlas_channel_id_r[4]]
        telescope_text_t = telescope_map[atlas_channel_id_t[4]]
    
    else:
        telescope_text_r = ''
        telescope_text_t = ''
        
    expression_r = f'R: {telescope_text_r} {type_map[atlas_channel_id_r[5]]} {mode_map[atlas_channel_id_r[6]]}'
    
    expression_t = f'T: {telescope_text_t} {type_map[atlas_channel_id_t[5]]} {mode_map[atlas_channel_id_t[6]]}'
            
    text = f"signal = {atlas_channel_id_r[:4].lstrip('0')}, {expression_r}, {expression_t}, {mode_map[atlas_channel_id_r[6]]}, dark-subtracted"
    
    return text

def header_time_text(start_timestamp, stop_timestamp, label = ""):
    
    start_date = start_timestamp.strftime("%d.%m.%Y")
    start_time = start_timestamp.strftime("%H:%M:%S")
    
    duration = (stop_timestamp - start_timestamp).total_seconds

    text = f"date of {label} measurement, time, duration of measurement = {start_date}, {start_time}UTC, {duration} s"

    return text

def header_radiosonde_text(rs_station_name, wmo_id, rs_start_timestamp):

    start_date = rs_start_timestamp.strftime("%d.%m.%Y")
    start_time = rs_start_timestamp.strftime("%H:%M:%S")
    
    text = f"location, WMO radiosonde station ID, date of radiosonde = {rs_station_name}, {wmo_id}, {start_date}, {start_time}UT"

    return text

def header_dark_text():
    
    line_1 = "range_RawSignal, RawSignal, range_BackgrCorrectedSignal, BackgrCorrectedSignal, range_RangeCorrectedSignal, RangeCorrectedSignal, range_RayleighRangeCorrectedSignal, RayleighRangeCorrectedSignal, RayleighRangeCorrectedSignalError, range_RayleighDarkRangeCorrectedSignal, RayleighDarkRangeCorrectedSignal, RayleighDarkRangeCorrectedSignalError"
    
    text = f"{line_1}"
    
    return text

def header_rayleigh_fit_text(norm_region):
    
    ray_l = np.round(norm_region[0], decimals = 1)
    ray_u = np.round(norm_region[1], decimals = 1)
    
    line_1 = f"lower and upper Rayleigh height limits = {ray_l}, {ray_u}"
    
    line_2 = "range, attnRayleighBSC, RangeCorrectedSignal"
    
    text = f"{line_1}\n{line_2}"
    
    return text

def header_telecover_text(iters, secs, extra_sec):
    
    extra = [f'{key}{iters+1}' for key in extra_sec.keys() if extra_sec[key]]
    
    iter_num = np.arange(1,iters+1,1)
    
    combo = []
    for num in iter_num:
        for sec in secs:
            combo.append(f'{sec}{num}')
    
    sec_text = ', '.join(combo+extra) 
    
    text = f"range, {sec_text}"
    
    return text

def header_polcal_text(G_R, G_T, H_R, H_T, K):
    
    line_1 = f"GR, GT, HR, HT, K = {G_R} {G_T} {H_R} {H_T} {K}"
    
    line_2 = "range, ITplus45, IRplus45, ITminus45, IRminus45, ITRayleigh, IRRayleigh"
    
    text = f"{line_1}\n{line_2}"

    return text
    
def sm_text(smooth, sm_lims, sm_win, sm_expo, flavour = ''):

    if flavour:
        caption = 'Smoothing'
    else:
        caption = f'Smoothing {flavour}'
        
    if smooth != True:
        return f"No Smoothing {flavour}"
    
    if isinstance(sm_win, (list, tuple, np.ndarray)):
        if len(sm_win):
            
            sm_lwin = np.round(float(sm_win[0]), decimals=3)
            sm_uwin = np.round(float(sm_win[-1]), decimals=3)
    
            if sm_lwin > sm_uwin:
                change = "Decrease"
            else:
                change = "Increase"
    
            if sm_expo == True:
                sm_type = "Exp"
            else:
                sm_type = "Lin"
                
            if len(sm_lims) == 0:
                sm_part = (
                    f"{caption}: All range, "
                    f"Win: {sm_lwin} km to {sm_uwin} km, "
                    f"{change}: {sm_type}"
                )
                
            else:
                sm_llim = np.round(float(sm_lims[0]), decimals=3)
                sm_ulim = np.round(float(sm_lims[-1]), decimals=3)
                sm_part = (
                    f"{caption}: {sm_llim} to {sm_ulim} km, "
                    f"Win: {sm_lwin} km to {sm_uwin} km, "
                    f"{change}: {sm_type}"
                )
        
        else:
            sm_part = "No {caption}"

    else:
        if sm_win:
            sm_win = np.round(float(sm_win), decimals=3)
            
            if len(sm_lims) == 0:
                sm_part = (
                    f"{caption}: All range, "
                    f"Win: {sm_win} km "
                )
                
            else:
                sm_llim = np.round(float(sm_lims[0]), decimals=3)
                sm_ulim = np.round(float(sm_lims[-1]), decimals=3)
                sm_part = (
                    f"{caption}: {sm_llim} to {sm_ulim} km, "
                    f"Win: {sm_win} km"
                )
       
        else:
            sm_part = "No {caption}"

    return sm_part

def sm_text_tlc(smooth, sm_win, nr_ulim):

    if smooth != True:
        return "No Smoothing"

    if sm_win:
        sm_win = np.round(float(sm_win), decimals=3)
        sm_part = (
            f"Smoothing: 0 to {nr_ulim} km, "
            f"Win: {sm_win} km - Win above: 0.5 km"
            )
   
    else:
        sm_part = "No Smoothing"

    return sm_part

def channel_text(atlas_channel_id, scc_channel_id=""):

    channel_part = str(atlas_channel_id)

    if scc_channel_id:
        channel_part = f"{channel_part} ({scc_channel_id})"

    return channel_part

def q_text(qa_test):

    return f'Quicklook {qa_dataset_labels[qa_test]}'

def stage_text(stage):
    
    return f'ATLAS Processing Stage: {stage}'

def system_text(lidar_name, station_name):
    
    parts = []

    if lidar_name:
        parts.append(str(lidar_name))

    if station_name:
        parts.append(str(station_name))

    channel_part = " ".join(parts)

    return channel_part


def channel_text_ratio(
    lidar_name,
    station,
    channel_r,
    channel_t,
    scc_channel_id_r="",
    scc_channel_id_t="",
):
    parts = []

    if lidar_name:
        parts.append(str(lidar_name))

    if station:
        parts.append(str(station))

    channel_part = " ".join(parts)

    channel_r_text = str(channel_r)
    channel_t_text = str(channel_t)

    if scc_channel_id_r:
        channel_r_text = f"{channel_r_text} ({scc_channel_id_r})"

    if scc_channel_id_t:
        channel_t_text = f"{channel_t_text} ({scc_channel_id_t})"

    ratio_part = f"{channel_r_text} to {channel_t_text}"

    if channel_part:
        return f"{channel_part} {ratio_part}"

    return ratio_part

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

def dateloc_text(start_timestamp, stop_timestamp, laser_pointing_angle):
    
    start_date = start_timestamp.strftime("%d.%m.%Y")
    start_time = start_timestamp.strftime("%H:%M:%S")
    stop_time = stop_timestamp.strftime("%H:%M:%S")
    
    laser_pointing_angle = np.round(float(laser_pointing_angle), decimals = 1)

    dateloc_part = (
        f"On {start_date} from {start_time} to {stop_time} UTC, "
        + r"$\nearrow$"
        + f"{laser_pointing_angle}"
        + r"$^{o}$ ZA"
    )  
    
    return dateloc_part

def iter_text(iters, sampling_time_per_sector):
    
    iter_part = f'Iterations: {iters}, Total sampling time per sector: {np.round(sampling_time_per_sector)} s'.strip()
    
    return(iter_part)

def mol_text(rs_format, rs_station_name, cloudnet_station_name, wmo_id, rs_start_timestamp):
        
    start_date = rs_start_timestamp.strftime("%d.%m.%Y")
    start_time = rs_start_timestamp.strftime("%H:%M:%S")
    
    if rs_format == 'ecmwf':
        station_name = cloudnet_station_name
    else:
        station_name = rs_station_name
        
    mol_part = f'{rs_format.capitalize()} {station_name} {start_date} {start_time}UT {wmo_id}'.strip()
    
    return mol_part
# -----------------------------------------------------------------------------
# Intercomparison text helpers
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class IntercomparisonLibraries:
    """Inputs required for intercomparison plot text generation.

    This is deliberately separate from ``Libraries`` above because an
    intercomparison plot represents several named entries, which may include
    multiple channels or pairs from the same dataset.
    """

    intercomparison_info: Dict[str, Any]
    group_id: str
    group_kind: str
    group: Dict[str, Any]


class GenerateIntercomparisonText:
    """Generate titles and filenames for intercomparison figures.

    Existing ``GenerateText`` methods are intentionally left untouched.
    """

    def __init__(self, lib: IntercomparisonLibraries):
        self.intercomparison_info = lib.intercomparison_info
        self.group_id = str(lib.group_id)
        self.group_kind = str(lib.group_kind)
        self.group = lib.group

    @staticmethod
    def _safe_token(value):
        text = str(value).strip().lower()
        chars = []
        previous_underscore = False
        for char in text:
            if char.isalnum() or char in {"-", "_"}:
                chars.append(char)
                previous_underscore = False
            elif not previous_underscore:
                chars.append("_")
                previous_underscore = True
        return "".join(chars).strip("_") or "group"

    def _reference_label(self):
        reference_entry = self.group.get("reference_entry")
        entry = self.group.get("entries", {}).get(reference_entry, {})
        if entry.get("entry_label") or entry.get("label"):
            return entry.get("entry_label") or entry.get("label")

        dataset_id = entry.get("dataset_id") or self.intercomparison_info.get("reference_dataset", "reference")
        dataset = self.intercomparison_info.get("datasets", {}).get(dataset_id, {})
        dataset_label = (
            dataset.get("dataset_label")
            or dataset.get("system_label")
            or dataset_id
        )
        product_id = entry.get("atlas_channel_id") or entry.get("atlas_pair_id")
        if product_id:
            return f"{dataset_label} - {product_id}"
        return dataset_label

    def make_title(self):
        general = self.intercomparison_info.get("general", {})
        label = self.group.get("label") or self.group_id
        reference_label = self._reference_label()
        vertical_scale = general.get("vertical_scale", "height_asl")
        plot_native = bool(general.get("plot_native_scale", False))

        if self.group_kind == "channel_group":
            plot_type = "Channel Intercomparison"
        elif self.group_kind == "pair_group":
            plot_type = "Pair Intercomparison"
        else:
            plot_type = "Intercomparison"

        if plot_native:
            scale_part = f"native {vertical_scale} scales"
        else:
            method = general.get("vertical_method", "")
            scale_part = f"{vertical_scale} | {method}"

        return (
            f"ATLAS {plot_type} - {label}\n"
            f"Reference: {reference_label} - Vertical: {scale_part}"
        )

    def make_filename(self):
        if self.group_kind == "channel_group":
            kind = "channel"
        elif self.group_kind == "pair_group":
            kind = "pair"
        else:
            kind = "group"

        return "_".join(
            [
                "intercomparison",
                kind,
                self._safe_token(self.group_id),
                "ATLAS",
                str(__version__),
            ]
        )
