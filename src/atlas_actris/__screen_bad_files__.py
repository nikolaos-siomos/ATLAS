#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:23:31 2024

@author: nikos
"""

import os, glob, warnings
from helper_functions import viewer_utils
import matplotlib
from processor.lidar_processing.signal import dark_correction
import numpy as np
import shutil
from pathlib import Path

warnings.filterwarnings('ignore')

#------------------------------------------------------------------------------
# A) Inputs
#------------------------------------------------------------------------------
# Path to the parent folder fs ATLAS 
output_folder = '/home/nikos/Big_Data/Intercomparison_Garmisch_2025/Analysis/mun/POLIS_9/9_9_9_20251107_20250412_002417'

input_folder = '/home/nikos/Big_Data/Network_drive_clone/mun/POLIS_9/9_9_9_20251107/nrm_02_screened'
# parent_folder = os.path.join('/home/nikos/Nextcloud/ACTRIS-CARS-LMU/Instruments (POLIS-6, POLIS-1064)/POLIS-1064/Laser/Interspersion_measurements_APD_flex/',
#                              '20250219','50us_pretrigger')

timescale = None # Set to None to skip averaging. Use e.g. 10s for 10 second averages, 30min for 30 minute averages, or 1h for 1 hour averages. Use 'all' to average all files. Use None to apply no averaging


smoothing_window = None # in meters. Set to None to skip smoothing

smoothing_range = None # in km. Set to None to smooth the whole profile. It will be applied only if the smoothing window is not None

normalization_range = None#[10., 20.] # in km. Set to None to skip normalization. All profiles will be normalized to their temporal mean

background_range = None # in km. Set to None to skip background correction of each profile. It will be ignored if the signal type is set to 'rangecor'

statistics_range = None # in km. Select a range to calculate statistics on signals to be displayed on the plots per channel. Statistics are calculated after averaging. Set to None to skip

dark_subtraction = True # Set to True to enable the dark signal subtraction. It can be applied only if the signal_type is 'raw', the mtype is not 'drk', and there is a drk file available in the netcdf folder

channels = None # use the ATLAS IDs to plot only specific channels, can be a list or scalar, set to None to plot all channels
    
colorscale = 'sequential' # select the colorscale ('sequential' or 'discrete') - If discrete is used 4 different colors will be applied iteratively for each measurement. If sequential is used, a colorscale of up to 256 colors will be used with the first/last measurement getting the first/last color repsectively. If there are more measurements than colors then the same colors might be used more than once for adjacent measurements   

custom_label = '' # provide a label to be added to the plot title. Use '' to skip

x_range = None # Add the default range in bins for the x axis, applied to all channels  Example: x_range = [0, 1000]

y_range_analog = None #[3.6, 4.1] # Add the default range in mv, applied to all analog channels. Example: y_range_analog = [0, 3]

y_range_photon = None # Add the default range in mv, applied to all analog channels.  Example: y_range_analog = [0, 3]

#------------------------------------------------------------------------------
# B) Calculations 
#------------------------------------------------------------------------------
# Path to the 'netcdf' folder that contains the files exported from ATLAS (more than 1 paths can be provided)
netcdf_folder = os.path.join(output_folder, 'netcdf')

# Path to the folder where the plots will be placed. If set to None
plot_folder = os.path.join(output_folder, 'html')

options = dict(netcdf_folder = netcdf_folder,
               plot_folder = plot_folder,
               timescale = timescale,
               smoothing_window = smoothing_window,
               smoothing_range = smoothing_range,
               normalization_range = normalization_range,
               background_range = background_range,
               statistics_range = statistics_range,
               dark_subtraction = dark_subtraction,
               mtype = 'ray',
               signal_type = 'raw',
               channels = channels,
               colorscale = colorscale,
               custom_label = custom_label)

options_rc = dict(netcdf_folder = netcdf_folder,
               plot_folder = plot_folder,
               timescale = timescale,
               smoothing_window = smoothing_window,
               smoothing_range = smoothing_range,
               normalization_range = normalization_range,
               background_range = background_range,
               statistics_range = statistics_range,
               dark_subtraction = dark_subtraction,
               mtype = 'ray',
               signal_type = 'rangecor',
               channels = channels,
               colorscale = colorscale,
               custom_label = custom_label)

viewer_utils.check_options(options)


#------------------------------------------------------------------------------
# C) Reading files 
#------------------------------------------------------------------------------
fpath_list = viewer_utils.get_fpaths(netcdf_folder = options['netcdf_folder'], 
                                     signal_type = options['signal_type'],
                                     mtype = options['mtype'])
for fpath in fpath_list:
    sig, ranges, date_info, stats, system_info, channel_info, time_info, shots = \
        viewer_utils.get_converter_signals(fpath = fpath,
                                           options = options)
    
fpath_list = viewer_utils.get_fpaths(netcdf_folder = options_rc['netcdf_folder'], 
                                     signal_type = options_rc['signal_type'],
                                     mtype = options_rc['mtype'])

for fpath in fpath_list:
    sig_rc, ranges_rc, date_info_rc, stats_rc = \
        viewer_utils.get_prepro_signals(fpath = fpath,
                                        options = options_rc)

#------------------------------------------------------------------------------
# E) Screening 
#------------------------------------------------------------------------------
def move_files(file_paths, dest_folder):
    dest_folder = Path(dest_folder)
    dest_folder.mkdir(parents=True, exist_ok=True)

    for src in file_paths:
        src = Path(src)
        if src.is_file():
            dest = dest_folder / src.name
            shutil.move(str(src), str(dest))  # move the file
            print(f"Moved: {src} → {dest}\n")
        else:
            print(f"WARNING: File not found: {src}")
            
bad_folder = os.path.join(os.path.dirname(input_folder), 'bad_files')

for j in range(sig_rc.channel.size):
    
    bad_indices = np.where(sig_rc[:,j,100:200].mean('bins')<0)[0]

    bad_fnames = time_info.Filename[bad_indices].values
    
    bad_filepaths = [os.path.join(input_folder,file) for file in bad_fnames]
    
    move_files(bad_filepaths, dest_folder = bad_folder)
    
        
    