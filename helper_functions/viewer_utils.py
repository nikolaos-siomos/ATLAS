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
import matplotlib.dates as mdates
from scipy.stats import linregress

def get_fpaths(parent_folders, raw_pattern = '_ray_ATLAS_', 
               prp_pattern = '_qck_drk_ATLAS_'): 
    
    fpath_raw = []
    fpath_prp = []
    
    for folder in parent_folders:
        raw_files = glob.glob(os.path.join(folder, 'converter', f'*{raw_pattern}*.nc'))
        fpath_raw.extend(raw_files)

        prp_files = glob.glob(os.path.join(folder, 'preprocessor', f'*{prp_pattern}*.nc'))
        fpath_prp.extend(prp_files)
    
    if len(fpath_raw) == 0:
        print(f"--Warning: No raw files matching the given patern ({raw_pattern}) were detected!")
    
    if len(fpath_raw) == 0:
        print(f"--Warning: No precessed files matching the given patern ({prp_pattern}) were detected!")
    
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

def make_title_mline(mtype, stype, timescale, smoothing_window, channel, date_part):
    
    if timescale != '': time_label = f'Avg: {timescale}'
    else: time_label = ''

    if not isinstance(smoothing_window, list): 
        sm_label = f'Sm Win: {int(smoothing_window)} m'
    else: sm_label = ''    
    
    ch_label = f'Ch: {channel}'
    
    label = ''
    for part in [mtype, stype, time_label, sm_label, ch_label]:
        if len(part) > 0:
            if len(label) == 0:
                label = label + part
            else:
                label = label + ', ' + part
    
    label = str.strip(label)
    
    if len(date_part) > 0:
        title = str.strip(f'{label}\n {date_part}')
        
    return(title)

def make_title_tseries(mtype, stype, timescale, smoothing_window, 
                       channel, date_part, xlims, x_units):
    
    if timescale != '': time_label = f'Avg: {timescale}'
    else: time_label = ''

    if not isinstance(smoothing_window, list): 
        sm_label = f'Sm Win: {int(smoothing_window)} m'
    else: sm_label = ''   
    
    if x_units == 'bins':
        extension = 'bins'
    else:
        extension = 'km'
        
    if len(xlims) > 0: 
        sum_label = f'Sum: {xlims[0]} - {xlims[1]} {extension}'
    else: sum_label = ''   
    
    ch_label = f'Ch: {channel}'
    
    label = ''
    for part in [mtype, stype, sum_label, time_label, sm_label, ch_label]:
        if len(part) > 0:
            if len(label) == 0:
                label = label + part
            else:
                label = label + ', ' + part
    
    label = str.strip(label)
    
    if len(date_part) > 0:
        title = str.strip(f'{label}\n {date_part}')
        
    return(title)

def make_filename(mtype, stype, timescale, smoothing_window, channel, 
                  start_date, start_time, end_time, ptype):
    
    if not isinstance(smoothing_window, list): 
        sm_label = f'{int(smoothing_window)}m'
    else: sm_label = ''    
        
    fname = ''
    for part in [start_date, start_time, end_time, mtype, stype, ptype,
                 timescale, sm_label, channel]:
        if len(part) > 0:
            if len(fname) == 0:
                fname = fname + part
            else:
                fname = fname + '_' + part

    return(fname)

def multi_line_plot(signal, ranges = [], channels = None, zero_bins = [],
                    resol = [], xlims = [], ylims = [],
                    mtype = '', stype = 'Raw', 
                    timescale = '', smoothing_window = [],
                    start_date = '', start_time = '', end_time = '',
                    exp_dir = '', dpi_val = 300):
    
    if channels == None:
        channels = signal.channel.values
    else:
        channels = np.array(channels)
    
    if len(ranges) == 0:
        x_units = 'bins'
    else:
        x_units = 'range' 
    
    if stype not in ['RC','Raw']:
        raise Exception(f"-- Error: The provided stype ({stype}) was not recognize. Please select one of: ['RC','Raw']")
            
    bins = signal.bins.values
    
    if len(start_date) > 0 and len(start_time) > 0 and len(end_time) > 0:
        date_part = date_text(start_date = start_date, 
                              start_time = start_time,
                              end_time = end_time) 
    else:
        date_part = ''
        
    for j in range(signal.channel.size):
        ch = signal.channel.values[j]

        if ch in channels:
            title = make_title_mline(mtype = mtype, stype = stype,
                                     timescale = timescale, 
                                     smoothing_window = smoothing_window, 
                                     channel = ch, date_part = date_part)
            
            fig = plt.figure(figsize=(6. , 3.))
            ax = fig.add_axes([0.1,0.15,0.86,0.70])
            
            ax.set_title(title)
            
            if x_units == 'range':
                x_vals = ranges[j,:].values
            else:
                x_vals = bins

            xbins = get_range_lims(x_vals = x_vals, xlims = xlims, 
                                   x_units = x_units)
                        
            if stype == 'Raw' and ch[6] == 'a':
                y_units = 'mV'
            elif stype == 'Raw' and ch[6] == 'p':
                y_units = 'Counts'
            elif stype == 'RC':
                y_units = 'AU'

            for i in range(signal.time.size):
                
                if x_units == 'bins': 
                    if len(xlims) > 0:
                        ax.plot(x_vals[xbins[0]:xbins[1] + 1], 
                                 signal[i, j, xbins[0]:xbins[1] + 1])
                        ax.set_xlim([xlims[0], xlims[1]])
                    else:
                        ax.set_plot(x_vals, signal[i, j, :])

                else:
                    if len(xlims) > 0:
                        ax.plot(x_vals[xbins[0]:xbins[1] + 1], 
                                 signal[i, j, xbins[0]:xbins[1] + 1])
                        ax.set_xlim([xlims[0], xlims[1]])
                    else:
                        ax.plot(x_vals, signal[i, j, :])

                if len(ylims) > 0:
                    ax.set_ylim([ylims[0],ylims[1]])
            
            if x_units == 'bins':
                ax.set_xlabel('Bins')                
                if stype == 'Raw' and not isinstance(zero_bins, list) and \
                    not isinstance(resol, list):
                    x_tick_bins = ax.get_xticks()
                    
                    x_tick_bins = x_tick_bins[x_tick_bins >= 0]
                    ax1 = ax.twiny()
                    
                    zero_bin = zero_bins.loc[ch]
                    resol_ch = resol.loc[ch]
                    
                    x1_ticks_l = np.round((x_tick_bins + zero_bin) * resol_ch * 1E-3,1)

                    ax1.set_xticks(x_tick_bins-x_tick_bins[0],x1_ticks_l)
                    ax1.tick_params(direction = "in", pad = -15)
                    ax1.set_xlabel('Range [km]',labelpad = -23)    
            else:
                ax.set_xlabel('Range [km]')    
                
            ax.set_ylabel(f'Signal [{y_units}]')                    
            ax.grid()
            ax.ticklabel_format(axis = 'y', useMathText = True, style='sci', scilimits=(-1,1))
                        
            if len(exp_dir) > 0 and os.path.exists(exp_dir):
                fname = make_filename(mtype = mtype, stype = stype,
                                      timescale = timescale, 
                                      smoothing_window = smoothing_window, 
                                      channel = ch, 
                                      start_date = start_date, 
                                      start_time = start_time, 
                                      end_time = end_time,
                                      ptype = 'mline')

                fpath = os.path.join(exp_dir, fname)
                
                # plt.tight_layout()
                
                fig.savefig(fpath, dpi = dpi_val)
            
            plt.show()
            fig.clf()
            
            plt.close()
    
    return()

def time_series_plot(signal, channels = None, ranges = [], 
                     xlims = [], ylims = [], mtype = '', stype = 'Raw', 
                     timescale = '', smoothing_window = [],
                     start_date = '', start_time = '', end_time = '',
                     exp_dir = '', dpi_val = 300):
    
    if channels == None:
        channels = signal.channel.values
    else:
        channels = np.array(channels)
    
    if len(xlims) == 0:
        x_units = 'bins'
    elif max(xlims) > 100:
        x_units = 'bins'
    else:
        x_units = 'range' 
    
    t_vals = signal.time.values
    
    if len(start_date) > 0 and len(start_time) > 0 and len(end_time) > 0:
        date_part = date_text(start_date = start_date, 
                              start_time = start_time,
                              end_time = end_time) 
    else:
        date_part = ''
        
    for j in range(signal.channel.size):
        ch = signal.channel.values[j]

        if ch in channels:
            title = make_title_tseries(mtype = mtype, stype = stype,
                                       timescale = timescale, 
                                       smoothing_window = smoothing_window, 
                                       channel = ch, date_part = date_part,
                                       xlims = xlims, x_units = x_units)
            
            # if y_units == None and 
            if x_units == 'range':
                x_vals = ranges[j,:].values
            else:
                x_vals = signal.bins.values
                
            if stype == 'Raw' and ch[6] == 'a':
                y_units = 'mV'
            elif stype == 'Raw' and ch[6] == 'p':
                y_units = 'Counts per shot'
            elif stype == 'RC':
                y_units = 'AU'
                
            xbins = get_range_lims(x_vals = x_vals, xlims = xlims, 
                                   x_units = x_units)
            
            signal_m = signal[:,j,xbins].copy().mean(dim = 'bins')
            
            sigma = signal_m.copy().std(dim = 'time').values
            mean = signal_m.copy().mean(dim = 'time').values

            fig = plt.figure(figsize=(15. , 3.))
            
            fig.suptitle(title)

            ax1 = fig.add_axes([0.045,0.15,0.64,0.70])
            ax2 = fig.add_axes([0.735,0.15,0.25,0.70])

            ax1.plot(t_vals, signal_m)
            
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

            ax1.set_xlabel('Time')        
            ax1.set_ylabel(f'Signal [{y_units}]')                    
            
            if len(ylims) > 0:
                ax1.set_ylim([ylims[0],ylims[1]])
                                 
            ax1.grid()
            ax1.ticklabel_format(axis = 'y', useMathText = True, style='sci', scilimits=(-1,1))
            
            t_vals_sec = (t_vals - t_vals[0])/ np.timedelta64(1, 's')
            slope, _, _, _, stderr = linregress(t_vals_sec, signal_m)
            
            if slope > 2. * stderr:
                slope_color = 'tab:red'
            else:
                slope_color = 'tab:green'
                
            ax1.text(0.1, 0.9, f"slope: {round_it(slope,3)} per s",
                     bbox=dict(facecolor=slope_color, alpha=0.1, zorder = 3),
                     transform = ax1.transAxes)

            ax2.hist(signal_m)
            
            ax2.set_xlabel(f'Signal [{y_units}]')                    
            ax2.set_ylabel('Counts')        
            
            if len(ylims) > 0:
                ax2.set_xlim([xlims[0],xlims[1]])
                     
            ax2.grid()
            ax2.ticklabel_format(axis = 'y', useMathText = True, style='sci', scilimits=(-1,1))
            ax2.ticklabel_format(axis = 'x', useMathText = True, style='sci', scilimits=(-1,1))
            
            ax2.text(0.1, 0.9, f"μ: {'{:.2E}'.format(mean)}",
                     bbox=dict(facecolor='tab:green', alpha=0.1, zorder = 3),
                     transform = ax2.transAxes)
            
            ax2.text(0.1, 0.8, f"σ: {'{:.2E}'.format(sigma)}",
                     bbox=dict(facecolor='tab:green', alpha=0.1, zorder = 3),
                     transform = ax2.transAxes)
            
            if len(exp_dir) > 0 and os.path.exists(exp_dir):
                fname = make_filename(mtype = mtype, stype = stype,
                                      timescale = timescale, 
                                      smoothing_window = smoothing_window, 
                                      channel = ch, 
                                      start_date = start_date, 
                                      start_time = start_time, 
                                      end_time = end_time,
                                      ptype = 'tseries')

                fpath = os.path.join(exp_dir, fname)
                
                # plt.tight_layout()
                
                fig.savefig(fpath, dpi = dpi_val)
            
            plt.show()
            fig.clf()
            
            plt.close()
    
    return()

def get_range_lims(x_vals, xlims, x_units):

    if len(xlims) > 0 and x_units == 'range':
        if len(np.where(x_vals < xlims[0])[0]) == 0:
            lbin = 0
        else:
            lbin = np.where(x_vals < xlims[0])[0][-1]

        if len(np.where(x_vals > xlims[1])[0]) == 0:
            ubin = 0
        else:
            ubin = np.where(x_vals > xlims[1])[0][0]
                                
        xbins = [lbin, ubin]
        
    elif len(xlims) > 0 and x_units == 'bins':
        xbins = xlims
    
    else:
        xbins = [0, -1]
        
    return(xbins)

def round_it(x, sig):
    
    if not np.isfinite(x) or np.isnan(x):
        x = -999.
        sig = 3
        
    if x != 0.:
        x_out = np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)
    else:
        x_out = 0.
        
    return x_out