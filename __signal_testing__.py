#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:23:31 2024

@author: nikos
"""

import os, glob, warnings
from helper_functions import testing_utils
import matplotlib

warnings.filterwarnings('ignore')

#------------------------------------------------------------------------------
# A) Inputs
#------------------------------------------------------------------------------
# Path to the parent folder fs ATLAS 
# parent_folder = '/home/nikos/Nextcloud3/ARS-stations/the/179_199_664_20240221'
parent_folder = '/home/nikos/Nextcloud3/ARS-stations/aky/138_146_437_20240427/'
# parent_folder = os.path.join('/home/nikos/Nextcloud/ACTRIS-CARS-LMU/Instruments (POLIS-6, POLIS-1064)/POLIS-1064/Laser/Interspersion_measurements_APD_flex/',
#                              '20250219','50us_pretrigger')

timescale = '10min' # Set to None to skip averaging. Use e.g. 10s for 10 second averages, 30min for 30 minute averages, or 1h for 1 hour averages. Use 'all' to average all files. Use None to apply no averaging

smoothing_window = None#1000.#1000. #1. # in meters. Set to None to skip smoothing

smoothing_range = None#[0., 60.] #[2, 60] # in km. Set to None to smooth the whole profile. It will be applied only if the smoothing window is not None

normalization_range = None#[10., 20.] # in km. Set to None to skip normalization. All profiles will be normalized to their temporal mean

background_range = None#[9., 60.] # in km. Set to None to skip background correction of each profile. It will be ignored if the signal type is set to 'rangecor'

statistics_range = [5.7, 7.3] # in km. Select a range to calculate statistics on signals to be displayed on the plots per channel. Statistics are calculated after averaging. Set to None to skip

mtype = 'ray' # set to either 'ray', 'drk', 'tlc', 'pcb' to plot the signals of the corresponding QA test 

channels = None# ['1064xpat'] # use the ATLAS IDs to plot only specific channels, can be a list or scalar, set to None to plot all channels
    
#------------------------------------------------------------------------------
# B) Calculations 
#------------------------------------------------------------------------------
# Path to the 'netcdf' folder that contains the files exported from ATLAS (more than 1 paths can be provided)
netcdf_folder = os.path.join(parent_folder, 'netcdf')

# Path to the folder where the plots will be placed. If set to None
plot_folder = os.path.join(parent_folder, 'html')

options = dict(netcdf_folder = netcdf_folder,
               plot_folder = plot_folder,
               timescale = timescale,
               smoothing_window = smoothing_window,
               smoothing_range = smoothing_range,
               normalization_range = normalization_range,
               background_range = background_range,
               statistics_range = statistics_range,
               mtype = mtype,
               channels = channels)

testing_utils.check_options(options)

#------------------------------------------------------------------------------
# C) Reading files 
#------------------------------------------------------------------------------
fpath_list_raw = testing_utils.get_fpaths(netcdf_folder = options['netcdf_folder'], 
                                          signal_type = 'raw',
                                          mtype = options['mtype'])

fpath_list_rangecor = testing_utils.get_fpaths(netcdf_folder = options['netcdf_folder'], 
                                               signal_type = 'rangecor',
                                               mtype = options['mtype'])
    
for fpath in fpath_list_raw:
    sig, ranges, date_info, stats, system_info, channel_info, time_info, shots = \
        testing_utils.get_converter_signals(fpath = fpath,
                                           options = options)
    
for fpath in fpath_list_rangecor:
    sig_rc, atb, ranges_rc, date_info_rc, stats_rc = \
        testing_utils.get_prepro_signals(fpath = fpath,
                                        options = options)

