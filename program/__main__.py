#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:48:11 2022

@author: nick
"""

import warnings, os, sys
from readers.parse_config import parse_config
from readers import read_files
from lidar_processing import short_prepro
from molecular import atmosphere

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the command line argument information
args = parse_config()

print('-----------------------------------------')
print('Start reading the QA file(s)...')
print('-----------------------------------------')

meas_info, channel_info, time_info, time_info_d, \
    sig_raw, sig_raw_d, shots, shots_d = \
        read_files.short_reader(args['input_file'])

meas_type = meas_info.Measurement_type
meas_label = {'ray' : 'Rayleigh', 
              'tlc' : 'Telecover', 
              'pcl' : 'Polarization Calibration'}
print(f'-- {meas_label[meas_type]} QA file succesfully parsed!')

if 'Rayleigh_File_Name' in meas_info.keys():
    ray_path = os.path.join(os.path.dirname(args['input_file']),  
                            meas_info['Rayleigh_File_Name'])
    
    meas_info_r, channel_info_r, time_info_r, time_info_dr, \
        sig_raw_r, sig_raw_dr, shots_r, shots_dr = \
            read_files.short_reader(ray_path)

    print(f"-- {meas_label['ray']} QA file succesfully parsed!")

print('-----------------------------------------')
print("")

        
if meas_type == 'drk' or not isinstance(sig_raw_d, list):
    drk, drk_pack, time_info_d = \
        short_prepro.dark(sig_raw = sig_raw_d, 
                          shots = shots_d, 
                          meas_info = meas_info, 
                          channel_info = channel_info, 
                          time_info = time_info_d,
                          external_info = args)

if 'Rayleigh_File_Name' in meas_info.keys() and not isinstance(sig_raw_d, list):
        drk_r, drk_pack_r, time_info_dr = \
            short_prepro.dark(sig_raw = sig_raw_dr, 
                              shots = shots_dr, 
                              meas_info = meas_info_r, 
                              channel_info = channel_info_r, 
                              time_info = time_info_dr,
                              external_info = args)
    
if meas_type == 'ray':
            
    sig, sig_pack, time_info = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = meas_type,
                              sig_drk = drk_pack['sig_flt'])

if meas_type == 'tlc':
    
    sig, sig_pack, time_info = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = meas_type,
                              sig_drk = drk_pack['sig_flt'])

if meas_type == 'pcl':
            
    sig, sig_pack, time_info = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = meas_type,
                              sig_drk = drk_pack['sig_flt'])

    sig_r, sig_pack_r, time_info_r = \
        short_prepro.standard(sig_raw = sig_raw_r, 
                              shots = shots_r, 
                              meas_info = meas_info_r, 
                              channel_info = channel_info_r, 
                              time_info = time_info_r,
                              external_info = args,
                              meas_type = 'ray',
                              sig_drk = drk_pack_r['sig_flt'])

if args['quicklook']:
            
    qck, qck_pack, time_info = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = 'qck',
                              sig_drk = drk_pack['sig_flt'])


if 'Sounding_File_Name' in meas_info.keys() or \
    ('Pressure_at_Lidar_Station' in meas_info.keys() and \
        'Temperature_at_Lidar_Station' in meas_info.keys()):
            molec, molec_info = \
                atmosphere.short_molec(heights = sig_pack['heights'],
                                       meas_info = meas_info, 
                                       channel_info = channel_info, 
                                       time_info = time_info,
                                       external_info = args)


    
    

# Quicklook --> 2 time scales (timeframes or file number) - be able to trim time dimension as argument
# Automate pretrigger bin flagging --> sequential removal of outliers
# Automate background calculations when a pretrigger does not exists
# Heymann check noise paper!