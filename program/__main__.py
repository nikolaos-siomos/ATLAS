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

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the command line argument information
args = parse_config()

meas_info, channel_info, time_info, time_info_d, \
    sig_raw, sig_raw_d, shots, shots_d = \
        read_files.short_reader(args['input_file'])

if meas_info.Measurement_type == 'drk' or not isinstance(sig_raw_d, list):
    drk, drk_pack, time_info_d = short_prepro.dark(sig_raw = sig_raw_d, 
                                                   shots = shots_d, 
                                                   meas_info = meas_info, 
                                                   channel_info = channel_info, 
                                                   time_info = time_info_d,
                                                   external_info = args)
    
if meas_info.Measurement_type == 'ray':
    sig, sig_pack, time_info = short_prepro.rayleigh(sig_raw = sig_raw, 
                                                     shots = shots, 
                                                     meas_info = meas_info, 
                                                     channel_info = channel_info, 
                                                     time_info = time_info,
                                                     external_info = args,
                                                     sig_drk = drk_pack['sig_flt'])

if meas_info.Measurement_type == 'tlc':
    sig, sig_pack, time_info = \
        short_prepro.telecover(sig_raw = sig_raw, 
                               shots = shots, 
                               meas_info = meas_info, 
                               channel_info = channel_info, 
                               time_info = time_info,
                               external_info = args,
                               sig_drk = drk_pack['sig_flt'])

if meas_info.Measurement_type == 'pcl':
    sig, sig_pack, time_info = \
        short_prepro.calibration(sig_raw = sig_raw, 
                                 shots = shots, 
                                 meas_info = meas_info, 
                                 channel_info = channel_info, 
                                 time_info = time_info,
                                 external_info = args,
                                 sig_drk = drk_pack['sig_flt'])

if args['quicklook']:
    qck, qck_pack = short_prepro.quicklook(sig_raw = sig_raw, 
                                           shots = shots, 
                                           meas_info = meas_info, 
                                           channel_info = channel_info, 
                                           external_info = args,
                                           sig_drk = drk_pack['sig_flt'])
    
    

 
# if iscr and external_info['debug']: pack_out['sig_ovf'] = sig.copy()

# Make the timeframes starting time exactly as they appear in the header
# Use a time scale separately of the xarray
# What happens when a square root channel is provided??? --> make a filter in scc_converter

#Quicklook --> 2 time scales (timeframes or file number) - be able to trim time dimension as argument
#Automate pretrigger bin flagging --> sequential removal of outliers
# Automate background calculations when a pretrigger does not exists
# The Q-switch trigger lasts 20μs --> average the dark above!
# Dark subtraction, doing it in the rangecorrected signal maybe adds some noise but probably it's ok
# Heymann check noise paper!