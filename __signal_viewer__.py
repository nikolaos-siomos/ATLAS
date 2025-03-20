#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:23:31 2024

@author: nikos
"""

import os, glob, warnings
from helper_functions import viewer_utils
import matplotlib

warnings.filterwarnings('ignore')

#------------------------------------------------------------------------------
# A) Inputs
#------------------------------------------------------------------------------
# Path to the parent folder fs ATLAS 
parent_folder = '/home/nikos/Nextcloud3/ARS-stations/the/179_199_664_20240221'
# parent_folder = os.path.join('/home/nikos/Nextcloud/ACTRIS-CARS-LMU/Instruments (POLIS-6, POLIS-1064)/POLIS-1064/Laser/Interspersion_measurements_APD_flex/',
#                              '20250219','50us_pretrigger')

timescale = '10min' # Set to None to skip averaging. Use e.g. 10s for 10 second averages, 30min for 30 minute averages, or 1h for 1 hour averages. Use 'all' to average all files. Use None to apply no averaging

smoothing_window = None#1000.#1000. #1. # in meters. Set to None to skip smoothing

smoothing_range = None#[0., 60.] #[2, 60] # in km. Set to None to smooth the whole profile. It will be applied only if the smoothing window is not None

normalization_range = None#[10., 20.] # in km. Set to None to skip normalization. All profiles will be normalized to their temporal mean

background_range = None#[9., 60.] # in km. Set to None to skip background correction of each profile. It will be ignored if the signal type is set to 'rangecor'

statistics_range = [5.7, 7.3] # in km. Select a range to calculate statistics on signals to be displayed on the plots per channel. Statistics are calculated after averaging. Set to None to skip

mtype = 'drk' # set to either 'ray', 'drk', 'tlc', 'pcb' to plot the signals of the corresponding QA test 

signal_type = 'raw' # select either 'raw' or 'rangecor' to plots the raw or the rangecorrected signals, respectively

channels = None# ['1064xpat'] # use the ATLAS IDs to plot only specific channels, can be a list or scalar, set to None to plot all channels
    
colorscale = 'discrete' # select the colorscale ('sequential' or 'discrete') - If discrete is used 4 different colors will be applied iteratively for each measurement. If sequential is used, a colorscale of up to 256 colors will be used with the first/last measurement getting the first/last color repsectively. If there are more measurements than colors then the same colors might be used more than once for adjacent measurements   

custom_label = 'pretrigger_50us' # provide a label to be added to the plot title. Use '' to skip

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
               signal_type = signal_type,
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
    
if options['signal_type'] == 'raw':
    for fpath in fpath_list:
        sig, ranges, date_info, stats, system_info, channel_info, time_info, shots = \
            viewer_utils.get_converter_signals(fpath = fpath,
                                               options = options)
    
if options['signal_type'] == 'rangecor':
    for fpath in fpath_list:
        sig, ranges, date_info, stats = \
            viewer_utils.get_prepro_signals(fpath = fpath,
                                            options = options)

#------------------------------------------------------------------------------
# E) Plotting 
#------------------------------------------------------------------------------
for k in range(len(fpath_list)):
    options['channels'] = channels

    titles = viewer_utils.make_title_mline(channels = sig.channel.values, 
                                           options = options,
                                           date_info = date_info) 

bname = \
    viewer_utils.make_plot_filename(fpath_list = fpath_list,
                                    options = options)
    
viewer_utils.multi_line_plot_bokeh(signal = sig, 
                                   signal_type = options['signal_type'],
                                   ranges = ranges,
                                   channels = sig.channel.values, 
                                   colorscale = options['colorscale'],
                                   titles = titles,
                                   plot_folder = plot_folder,
                                   filename = bname[k],
                                   stats = stats,
                                   statistics_range = options['statistics_range'])
