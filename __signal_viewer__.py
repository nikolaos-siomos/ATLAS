#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:23:31 2024

@author: nikos
"""

import os, glob, warnings
from helper_functions import viewer_utils

warnings.filterwarnings('ignore')

#------------------------------------------------------------------------------
# A) Inputs
#------------------------------------------------------------------------------

# parent_folders = ['/home/nikos/Nextcloud3/ARS-stations/the/179_199_664_20231221/netcdf',
#                   '/home/nikos/Nextcloud3/ARS-stations/the/179_199_664_20240116/netcdf',
#                   '/home/nikos/Nextcloud3/ARS-stations/the/179_199_665_20230930/netcdf',
#                   '/home/nikos/Nextcloud3/ARS-stations/the/179_199_665_20231221/netcdf',
#                   '/home/nikos/Nextcloud3/ARS-stations/the/179_199_664_20240221/netcdf',
#                   '/home/nikos/Nextcloud3/ARS-stations/the/179_199_665_20240224/netcdf']                                 

parent_folders = ['/home/nikos/Nextcloud3/ARS-stations/tst/8_8_11_20240705/01/netcdf/']                                 
# parent_folders = ['/home/nikos/Nextcloud3/ARS-stations/the/archive/179_199_664_20240221/netcdf/']                                 
# parent_folders = ['/home/nikos/Nextcloud3/ARS-stations/brc/123_202_812_20240205/netcdf/']                                 

timescale = ''
smoothing_window = 0.

raw_pattern = '_ray_ATLAS_'
prp_pattern = '_qck_ray_ATLAS_'

raw_mtype = 'Ray'
prp_mtype = 'Ray'

# channels = ['0354xrax', '0354xrpx', '0355xtat', '0355xtpt', '1064xtax', '0530xrax',
#         '0530xrpx', '0532xtat', '0532xtpt', '0532xcar', '0532xcpr', '0355xcar',
#         '0355xcpr'] #can be a list or scalar

channels = ['0355xpar', '0355xppr', '0355xcat', '0355xcpt', '0387xvan', '0387xvpn'] #can be a list or scalar

# xlims in bins or km
raw_mline_xlims = [0, 3000]
prp_mline_xlims = [0, 3000]

raw_tseries_xlims = [2000, 3000]
prp_tseries_xlims = [5.5, 7.5]

raw_ylims = [0,5]
prp_ylims = [-2E6,10E6]

raw_plot_mline = True
raw_plot_tseries = False

prp_plot_mline = False
prp_plot_tseries = False

exp_dir = '/home/nikos/Nextcloud3/ARS-stations/tst/raw_sig_test'
# exp_dir = ''

#------------------------------------------------------------------------------
# B) Calculations 
#------------------------------------------------------------------------------
fpath_raw, fpath_prp = viewer_utils.get_fpaths(parent_folders, 
                                               raw_pattern = raw_pattern,
                                               prp_pattern = prp_pattern)

for fpath in fpath_raw:
    system_info, channel_info, time_info, time_info_d,\
        sig, sig_d, shots, shots_d = viewer_utils.get_converter_signals(fpath)

for fpath in fpath_prp:
    sig_rc, std_rc, ranges, \
        start_date, start_time, end_time = \
            viewer_utils.get_prepro_signals(fpath, timescale = timescale, 
                                            smoothing_window = smoothing_window)

for fpath in fpath_raw:
    if raw_plot_mline:
        viewer_utils.multi_line_plot(signal = sig, 
                                     xlims = raw_mline_xlims,
                                     ylims = raw_ylims,
                                     channels = channels, 
                                     zero_bins = channel_info.DAQ_Trigger_Offset,
                                     resol = channel_info.Raw_Data_Range_Resolution,
                                     mtype = raw_mtype,
                                     stype = 'Raw',
                                     start_date = system_info.RawData_Start_Date,
                                     start_time = system_info.RawData_Start_Time_UT,
                                     end_time = system_info.RawData_Stop_Time_UT,
                                     exp_dir = exp_dir)
        
    if raw_plot_tseries:
        viewer_utils.time_series_plot(signal = sig, 
                                      xlims = raw_tseries_xlims,
                                      ylims = [],
                                      channels = channels, 
                                      mtype = raw_mtype,
                                      stype = 'Raw',
                                      start_date = system_info.RawData_Start_Date,
                                      start_time = system_info.RawData_Start_Time_UT,
                                      end_time = system_info.RawData_Stop_Time_UT,
                                      exp_dir = exp_dir)
        

for fpath in fpath_prp:
    if prp_plot_mline:
        viewer_utils.multi_line_plot(signal = sig_rc, 
                                     ranges = ranges, 
                                     xlims = prp_mline_xlims,
                                     ylims = prp_ylims,
                                     channels = channels, 
                                     mtype = prp_mtype,
                                     stype = 'RC',
                                     timescale = timescale,
                                     smoothing_window = smoothing_window,
                                     start_date = start_date,
                                     start_time = start_time,
                                     end_time = end_time,
                                     exp_dir = exp_dir)

    if prp_plot_tseries:
        viewer_utils.time_series_plot(signal = sig_rc, 
                                      xlims = prp_tseries_xlims,
                                      ylims = [],
                                      ranges = ranges, 
                                      channels = channels, 
                                      mtype = raw_mtype,
                                      stype = 'RC',
                                      start_date = system_info.RawData_Start_Date,
                                      start_time = system_info.RawData_Start_Time_UT,
                                      end_time = system_info.RawData_Stop_Time_UT,
                                      exp_dir = exp_dir)