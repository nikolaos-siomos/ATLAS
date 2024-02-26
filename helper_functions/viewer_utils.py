#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:41:41 2024

@author: nikos
"""

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from processor.readers.read_files import short_reader
import netCDF4 as nc
import os, glob
from visualizer.tools.smoothing import sliding_average_2D_fast

def get_fpaths(parent_folders, raw_pattern = '_ray_ATLAS_', 
               prp_pattern = '_qck_drk_ATLAS_'): 
    
    fpath_raw = []
    fpath_prp = []
    
    for folder in parent_folders:
        raw_files = glob.glob(os.path.join(folder, 'converter', f'*{raw_pattern}*.nc'))
        fpath_raw.extend(raw_files)

        prp_files = glob.glob(os.path.join(folder, 'preprocessor', f'*{prp_pattern}*.nc'))
        fpath_prp.extend(prp_files)
    
    return(fpath_raw, fpath_prp)

def get_converter_signals(fpath_raw):
    if len(fpath_raw) > 0 and os.path.exists(fpath_raw):
        system_info, channel_info, time_info, time_info_d,\
            sig, sig_d, shots, shots_d = \
                   short_reader(fpath = fpath_raw, 
                                exclude_telescope_type = '', 
                                exclude_channel_type = '', 
                                exclude_acquisition_mode = '', 
                                exclude_channel_subtype = '', 
                                use_channels = None)
    else:
        print("-- Warning: Raw file not found!")
        
    return(system_info, channel_info, time_info, time_info_d,\
           sig, sig_d, shots, shots_d)

def get_prepro_signals(fpath_prp, timescale = '',
                       smoothing_window = [], smoothing_llim = None, 
                       smoothing_ulim = None):
            
    if len(fpath_prp) > 0 and os.path.exists(fpath_prp):
        data = xr.open_dataset(fpath_prp)
        
        # heights = data.Height_levels / 1000.
        ranges = data.Range_levels / 1000.
        
        # Extract signal
        sig_rc = data.Range_Corrected_Signals
        sig_rc = sig_rc.copy().where(sig_rc != nc.default_fillvals['f8'])
        
        start_date = data.RawData_Start_Date
        start_time = data.RawData_Start_Time_UT
        end_time = data.RawData_Stop_Time_UT

        std_rc = np.nan * sig_rc.copy()
        
        if len(timescale) > 0:
            sig_rc = sig_rc.copy().resample(time = timescale).mean()
            std_rc = std_rc.copy().resample(time = timescale).mean()
    
        if not isinstance(smoothing_window, list):
            for j in range(sig_rc.channel.size):
                if smoothing_llim == None:
                    smoothing_llim = ranges[j,0].copy().values
                if smoothing_ulim == None:
                    smoothing_ulim = ranges[j,-1].copy().values
                        
                sig_sm, sig_std = \
                    sliding_average_2D_fast(z_vals = sig_rc[:,j,:].copy().values, 
                                            y_vals = ranges[j,:].copy().values,
                                            y_sm_lims = [smoothing_llim, smoothing_ulim],
                                            y_sm_win = smoothing_window,
                                            err_type = 'std')
                sig_rc[:,j,:] = sig_sm
                std_rc[:,j,:] = sig_std
        else:
            std_rc = np.nan * sig_rc.copy()
    else:
        print("-- Warning: Preprocessed file not found!")

    return(sig_rc, std_rc, ranges, start_date, start_time, end_time)

def date_text(start_date, start_time, end_time):
    
    date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
    
    start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'

    end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
    
    date_part = f'On {date} from {start} to {end} UTC'.strip()
    
    return(date_part)

def make_title(mtype, timescale, smoothing_window, channel, date_part):
    
    if timescale != '': time_label = f'Avg: {timescale}'
    else: time_label = ''

    if not isinstance(smoothing_window, list): 
        sm_label = f'Sm Win: {int(smoothing_window)} m'
    else: sm_label = ''    
    
    ch_label = f'Ch: {channel}'
    
    label = ''
    for part in [mtype, time_label, sm_label, ch_label]:
        if len(part) > 0:
            if len(label) == 0:
                label = label + part
            else:
                label = label + ', ' + part
    
    label = str.strip(label)
    
    if len(date_part) > 0:
        title = str.strip(f'{label}\n {date_part}')
        
    return(title)

def make_filename(mtype, timescale, smoothing_window, channel, 
                  start_date, start_time, end_time):
    
    mtype_label = mtype.replace(' ','_')

    if not isinstance(smoothing_window, list): 
        sm_label = f'{int(smoothing_window)}m'
    else: sm_label = ''    
        
    fname = ''
    for part in [start_date, start_time, end_time, mtype_label, 
                 timescale, sm_label, channel]:
        if len(part) > 0:
            if len(fname) == 0:
                fname = fname + part
            else:
                fname = fname + '_' + part

    return(fname)

def multi_line_plot(signal, x_vals = [], channels = None,
                    xlims = [], ylims = [], 
                    mtype = '', timescale = '', smoothing_window = [],
                    start_date = '', start_time = '', end_time = '',
                    exp_dir = '', dpi_val = 300):
    
    if channels == None:
        channels = signal.channel.values
    else:
        channels = np.array(channels)
    
    if len(x_vals) == 0:
        x_vals = signal.bins.values
        x_units = 'bins'
    else:
        x_units = 'range' 
    
    if len(start_date) > 0 and len(start_time) > 0 and len(end_time) > 0:
        date_part = date_text(start_date = start_date, 
                              start_time = start_time,
                              end_time = end_time) 
    else:
        date_part = ''
        
    for j in range(signal.channel.size):
        ch = signal.channel.values[j]

        if ch in channels:
            title = make_title(mtype = mtype, timescale = timescale, 
                               smoothing_window = smoothing_window, 
                               channel = ch, date_part = date_part)
            plt.title(title)

            if len(xlims) > 0 and x_units == 'range':
                if len(np.where(x_vals[j,:] < xlims[0])[0]) == 0:
                    lbin = 0
                else:
                    lbin = np.where(x_vals[j,:] < xlims[0])[0][-1]

                if len(np.where(x_vals[j,:] > xlims[1])[0]) == 0:
                    ubin = 0
                else:
                    ubin = np.where(x_vals[j,:] > xlims[1])[0][0]
                                        
                xbins = [lbin, ubin]
                
            elif len(xlims) > 0 and x_units == 'bins':
                xbins = xlims
                         
            for i in range(signal.time.size):
                
                if x_units == 'bins': 
                    if len(xlims) > 0:
                        plt.plot(x_vals[xbins[0]:xbins[1] + 1], 
                                 signal[i, j, xbins[0]:xbins[1] + 1])
                        plt.xlim([xlims[0], xlims[1]])
                    else:
                        plt.plot(x_vals, signal[i, j, :])
                    plt.xlabel('Bins')

                else:
                    if len(xlims) > 0:
                        plt.plot(x_vals[j, xbins[0]:xbins[1] + 1], 
                                 signal[i, j, xbins[0]:xbins[1] + 1])
                        plt.xlim([xlims[0], xlims[1]])
                    else:
                        plt.plot(x_vals[j, :], signal[i, j, :])
                    plt.xlabel('Range [km]')        

                if len(ylims) > 0:
                    plt.ylim([ylims[0],ylims[1]])
                plt.ylabel('Signal [AU]')                    
            plt.grid()
            
            if len(exp_dir) > 0 and os.path.exists(exp_dir):
                fname = make_filename(mtype, timescale = timescale, 
                                      smoothing_window = smoothing_window, 
                                      channel = ch, 
                                      start_date = start_date, 
                                      start_time = start_time, 
                                      end_time = end_time)

                fpath = os.path.join(exp_dir, fname)
                
                plt.tight_layout()
                
                plt.savefig(fpath, dpi = dpi_val)
            
            plt.show()
            plt.clf()
            
            plt.close()
    
    return()