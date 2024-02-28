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

parent_folders = ['/home/nikos/Nextcloud3/ARS-stations/the/179_199_664_20231221/netcdf',
                  '/home/nikos/Nextcloud3/ARS-stations/the/179_199_664_20240116/netcdf',
                  '/home/nikos/Nextcloud3/ARS-stations/the/179_199_665_20230930/netcdf',
                  '/home/nikos/Nextcloud3/ARS-stations/the/179_199_665_20231221/netcdf']                                 

#parent_folders = ['/home/nikos/Nextcloud3/ARS-stations/gar/203_246_823_20240222/netcdf/']                                 

timescale = ''
smoothing_window = 4000.

raw_pattern = '_ray_ATLAS_'
prp_pattern = '_qck_ray_ATLAS_'

raw_mtype = 'Raw Ray'
prp_mtype = 'RC Ray'

raw_channels = '0532xcpr' #can be a list or scalar
prp_channels = '1064xtax' #can be a list or scalar

raw_xlims = [0,9000] #in bins
prp_xlims = [0, 20] #in kms

raw_ylims = []
prp_ylims = [-0.4E7,0.4E7]

# exp_dir = '/home/nikos/Nextcloud3/ARS-stations/the/dark_distortions'
exp_dir = ''

#------------------------------------------------------------------------------
# B) Calculations 
#------------------------------------------------------------------------------
fpath_raw, fpath_prp = viewer_utils.get_fpaths(parent_folders, 
                                               raw_pattern = raw_pattern,
                                               prp_pattern = prp_pattern)

for fpath in fpath_raw:
    system_info, channel_info, time_info, time_info_d,\
        sig, sig_d, shots, shots_d = viewer_utils.get_converter_signals(fpath)

    viewer_utils.multi_line_plot(signal = sig, 
                                 xlims = raw_xlims,
                                 ylims = raw_ylims,
                                 channels = raw_channels, 
                                 mtype = raw_mtype,
                                 start_date = system_info.RawData_Start_Date,
                                 start_time = system_info.RawData_Start_Time_UT,
                                 end_time = system_info.RawData_Stop_Time_UT,
                                 exp_dir = exp_dir)

for fpath in fpath_prp:
    sig_rc, std_rc, ranges, \
        start_date, start_time, end_time = \
            viewer_utils.get_prepro_signals(fpath, timescale = timescale, 
                                            smoothing_window = smoothing_window)

    viewer_utils.multi_line_plot(signal = sig_rc, 
                                 x_vals = ranges, 
                                 xlims = prp_xlims,
                                 ylims = prp_ylims,
                                 channels = prp_channels, 
                                 mtype = prp_mtype,
                                 timescale = timescale,
                                 smoothing_window = smoothing_window,
                                 start_date = start_date,
                                 start_time = start_time,
                                 end_time = end_time,
                                 exp_dir = exp_dir)
