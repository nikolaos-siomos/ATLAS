#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:48:11 2022

@author: nick
"""

import warnings, os, sys
from readers.parse_config import main_parser
from readers import read_files
from lidar_processing import short_prepro
from lidar_molecular import atmosphere
from export import nc_dataset
from version import __version__

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the command line argument information
args = main_parser()

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
    
    if all([ch in channel_info_r.index.values 
            for ch in channel_info.index.values]):
        print(f"-- {meas_label['ray']} QA file used for calibration "+
              "succesfully parsed!")
    else:
        sys.exit("-- Error! The Rayleigh measurement accompaning the " +
                 "calibration measurement {meas_info['Rayleigh_File_Name']} " +
                 "lacks at least of the calibration measurement channels! " +
                 "Please make sure that the Rayleigh measuremen contains " +
                 "at least the channels channels necessary for the calibration")

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
    sig_drk = drk_pack['sig_flt']
else:
    sig_drk = []

if 'Rayleigh_File_Name' in meas_info.keys() and not isinstance(sig_raw_d, list):
    drk_r, drk_pack_r, time_info_dr = \
        short_prepro.dark(sig_raw = sig_raw_dr, 
                          shots = shots_dr, 
                          meas_info = meas_info_r, 
                          channel_info = channel_info_r, 
                          time_info = time_info_dr,
                          external_info = args)
    sig_drk_r = drk_pack_r['sig_flt']
else:
    sig_drk_r = []
        
if meas_type == 'ray':
            
    ray, ray_pack, time_info_ray = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = meas_type,
                              sig_drk = sig_drk)
    
    molec, molec_info, meteo = \
        atmosphere.short_molec(heights = ray_pack['heights'],
                               ranges = ray_pack['ranges'],
                               meas_info = meas_info, 
                               channel_info = channel_info, 
                               time_info = time_info_ray,
                               external_info = args)
    
    nc_ray = nc_dataset.rayleigh(sig = ray, 
                                 molec = molec,
                                 meteo = meteo,
                                 meas_info = meas_info, 
                                 channel_info = channel_info, 
                                 time_info = time_info_ray, 
                                 molec_info = molec_info, 
                                 heights = ray_pack['heights'], 
                                 ranges = ray_pack['ranges'],
                                 version = __version__,
                                 dir_out = args['output_folder'])

if meas_type == 'tlc':
    
    tlc, tlc_pack, time_info_tlc = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = meas_type,
                              sig_drk = sig_drk)

    nc_tlc = nc_dataset.telecover(sig = tlc, 
                                  meas_info = meas_info, 
                                  channel_info = channel_info, 
                                  time_info = time_info_tlc, 
                                  heights = tlc_pack['heights'], 
                                  ranges = tlc_pack['ranges'],
                                  version = __version__,
                                  dir_out = args['output_folder'])

if meas_type == 'pcl':
            
    pcl, pcl_pack, time_info_pcl = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = meas_type,
                              sig_drk = sig_drk)
    
    ray, ray_pack, time_info_ray = \
        short_prepro.standard(sig_raw = sig_raw_r, 
                              shots = shots_r, 
                              meas_info = meas_info_r, 
                              channel_info = channel_info_r, 
                              time_info = time_info_r,
                              external_info = args,
                              meas_type = 'ray',
                              sig_drk = sig_drk_r)
    
    molec, molec_info, meteo = \
        atmosphere.short_molec(heights = ray_pack['heights'],
                               ranges = ray_pack['ranges'],
                               meas_info = meas_info_r, 
                               channel_info = channel_info_r, 
                               time_info = time_info_ray,
                               external_info = args)
        
    pcl_ch = pcl.channel.values
    
    nc_pcl = nc_dataset.calibration(sig = pcl, 
                                    sig_ray = ray.loc[dict(channel = pcl_ch)],
                                    molec = molec.loc[dict(channel = pcl_ch)],
                                    meteo = meteo.loc[dict(channel = pcl_ch)],
                                    meas_info = meas_info, 
                                    meas_info_ray = meas_info_r, 
                                    channel_info = channel_info, 
                                    channel_info_ray = channel_info.loc[pcl_ch], 
                                    time_info = time_info_pcl, 
                                    time_info_ray = time_info_ray, 
                                    molec_info = molec_info, 
                                    heights = pcl_pack['heights'], 
                                    ranges = pcl_pack['ranges'],
                                    ranges_ray = ray_pack['ranges'],
                                    heights_ray = ray_pack['heights'], 
                                    version = __version__,
                                    dir_out = args['output_folder'])

if args['quicklook']:
            
    qck, qck_pack, time_info_qck = \
        short_prepro.standard(sig_raw = sig_raw, 
                              shots = shots, 
                              meas_info = meas_info, 
                              channel_info = channel_info, 
                              time_info = time_info,
                              external_info = args,
                              meas_type = 'qck',
                              sig_drk = sig_drk)
        
    nc_qck = nc_dataset.quicklook(sig = qck, 
                                 meas_info = meas_info, 
                                 channel_info = channel_info, 
                                 time_info = time_info_qck, 
                                 heights = qck_pack['heights'], 
                                 ranges = qck_pack['ranges'],
                                 version = __version__,
                                 dir_out = args['output_folder'])

    
    

# Quicklook --> 2 time scales (timeframes or file number) - be able to trim time dimension as argument
# Automate pretrigger bin flagging --> sequential removal of outliers
# Automate background calculations when a pretrigger does not exists
# Heymann check noise paper!