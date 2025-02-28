#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:41:41 2024

@author: nikos
"""

import xarray as xr
import numpy as np
import pandas as pd
from processor.readers.read_files import short_reader
import netCDF4 as nc
import os, glob
from visualizer.tools.smoothing import sliding_average_2D_bin_fast, sliding_average_2D_fast
from scipy.stats import linregress, shapiro
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import TabPanel, Tabs, LinearAxis, Range1d, BoxAnnotation, Label
from bokeh.palettes import Category10, turbo
import webbrowser 

# https://docs.bokeh.org/en/latest/docs/reference/palettes.html#bokeh-palettes

def check_options(options, ranges = None):
    
    none_allowed = dict(netcdf_folder = False,
                        plot_folder = False,
                        timescale = True,
                        smoothing_window = True,
                        background_range = True,
                        normalization_range = True,
                        mtype = False,
                        signal_type = False,
                        channels = True,
                        colorscale = False,
                        statistics_range = True,
                        custom_label = False)
    
    option_type = dict(netcdf_folder = [str],
                       plot_folder = [str],
                       timescale = [str],
                       smoothing_window = [int, float],
                       background_range =[int, float],
                       normalization_range = [int, float],
                       mtype = [str],
                       signal_type = [str],
                       channels = [str],
                       colorscale = [str],
                       statistics_range = [int, float],
                       custom_label = [str])
    
    option_shape = dict(netcdf_folder = 'scalar',
                        plot_folder = 'scalar',
                        timescale = 'scalar',
                        smoothing_window = 'scalar',
                        background_range = 'list',
                        normalization_range = 'list',
                        mtype = 'scalar',
                        signal_type = 'scalar',
                        channels = 'list',
                        colorscale = 'scalar',
                        statistics_range = 'list',
                        custom_label = 'scalar')
    
    
    check_none(options = options, none_allowed = none_allowed)
    
    check_list(options = options, option_shape = option_shape)
    
    check_option_type(options = options, option_type = option_type)    
    
    check_path(options  = options, key = 'netcdf_folder')
    check_parent_directory(options  = options, key = 'plot_folder')
            
    check_channel_format(options)
        
    allowed_signal_types = ['raw', 'rangecor']
    allowed_mtypes = ['ray', 'tlc', 'pcb', 'drk']
    colorbar = ['sequential', 'discreet']
    check_allowed_values(options = options,
                         keys = ['signal_type', 'mtype', 'colorscale'],
                         allowed_values = [allowed_signal_types, allowed_mtypes, colorbar])
                         
    check_conflicts(options = options)

    if not type(ranges) == type(None):
        check_smoothing_window(options = options, 
                               max_range = ranges[-1])
        
        for key in ['smoothing_range', 'normalization_range', 
                    'background_range', 'statistics_range']:
            check_range_limits(options = options, 
                               key = key, 
                               min_range = ranges[0], 
                               max_range = ranges[-1])
    
    return()

def check_none(options, none_allowed):
    
    for key in none_allowed.keys():
        if options[key] == None and none_allowed[key] == False:
            raise Exception(f"--Error: The {key} option cannot be None.")

    return()  

def check_list(options, option_shape):
    
    for key in option_shape.keys():
        if not options[key] == None:
            if option_shape[key] == 'list' and not isinstance(options[key], list) :
                raise Exception(f"--Error: The {key} option should be a list")
            if option_shape[key] == 'scalar' and not np.isscalar(options[key]):
                raise Exception(f"--Error: The {key} option should be a scalar")

    return()    

def check_option_type(options, option_type):
    
    for key in option_type.keys():
        if not options[key] == None:
            if isinstance(options[key], list):
                for opt in options[key]:
                    if type(opt) not in option_type[key]:
                        raise Exception(f"--Error: All of the the {key} option elements should be {option_type[key]}. Type detected: {type(opt)}")
            if np.isscalar(options[key]):
                if type(options[key]) not in option_type[key]:
                    raise Exception(f"--Error: The {key} option should be an {option_type[key]}. Type detected: {type(opt)}")

    return()           
    

def check_scalar_type(option, key):
    
    if not np.isscalar(option[key]):
        raise Exception(f"--Error: The {key} option must be a scalar, not a list or array")
    
    return()
    
def check_channel_format(options):
    
    if options['channels'] != None:
        if not isinstance(options['channels'], list):
            raise Exception("--Error: The channels option must always be provided as either list or None. Please provide a list even if a single channel is provided, e.g. ['1064xcar']")
        
        for ch in options['channels']:
            if len(ch) != 8:
                raise Exception(f"--Error: The provided channel {ch} has a wrong ID. Please use a string with lenght 8 where:\n-- the first 4 letters are wavelength e.g. '0355'\n--the 5th letter is the telescope\n--the 6th letter is the channel type\n--the 7th letter is the channel mode ('a' or 'p')\n--the 8th letter is the channel subtype\nPlease refer to the ATLAS manual for more information about the ATLAS channel ID")

    return()

def check_conflicts(options):
   
    if not options['normalization_range'] == None and options['signal_type'] == 'raw':
        print("--Warning: Normalization is not supported for raw signals. The provided normalization_range will be ignored")

    if not options['background_range'] == None and options['signal_type'] == 'rangecor':
        print("--Warning: Background correction is not supported for rangecorrected signals. The provided background_range will be ignored")
    
    return()
    

def check_allowed_values(options, keys, allowed_values):

    for key, values in zip(keys, allowed_values):
        if options[key] not in values:
            raise Exception(f"--Error: The provided {key} parameter is not supported:\n{options[key]}\nPlease select one of: {values}")
    
    return()

def check_path(options, key):
    
    if not os.path.exists(options[key]):
        raise Exception(f"-- Warning: The provided {key} does not exist:\n{options[key]}") 

    return()

def check_parent_directory(options, key):
    
    if not options[key] == None and not os.path.exists(os.path.dirname(os.path.normpath(options[key]))):
        raise Exception(f"-- The provided {key} is not placed in a valid directory:\n{options[key]}\nPlease use a valid directory to create the {key} or use a valid existing path for the ")
    
    return()

def check_smoothing_window(options, max_range):
    
    if not options['smoothing_window'] == None:

        if options['smoothing_window'] < 0:
            raise Exception("--Error: The smoothing window must be positive")

        if options['smoothing_window'] < 7.:
            raise Exception("--Error: The smoothing window must be provided in meters")
                
        if options['smoothing_window'] >= 1e3 * max_range:
            raise Exception(f"--Error: The smoothing window cannot be higher than the maximum range of the channel ({1e-3 * max_range} km). ")      
        
    return()

def check_range_limits(options, key, min_range, max_range):
    
    if not options[key] == None:
        
        if len(options[key]) != 2:
            raise Exception(f"--Error: The key option must be a list of exaclty 2 elements, the lower and the upper limit. A {len(options[key])} element list was provided")
            
        llim = options[key][0]
        ulim = options[key][1]
            
        if llim < min_range or ulim > max_range or \
            llim > max_range or ulim < min_range:
            raise Exception(f'--Error: The provided {key} is out of the signal limits:\n{options[key]}\nPlease provide values:\n    between {min_range} and {max_range}\n    or None to smooth over the whole profile')
    
        if llim > ulim:
            raise Exception(f'--Error: The lower limit of the provided {options[key]} is higher than the upper limit: {options[key]}')
            
    return()

def get_fpaths(netcdf_folder, signal_type, mtype): 
        
    if signal_type == 'raw':   
        pattern = f"_{mtype}_ATLAS_"
        files = glob.glob(os.path.join(netcdf_folder, 'converter', f'*{pattern}*.nc'))
        
        if len(files) == 0:
            raise Exception(f"--Warning: No raw files matching the given patern {pattern} were detected!")

    elif signal_type == 'rangecor':   
        prp_pattern = f"_qck_{mtype}_ATLAS_"
        files = glob.glob(os.path.join(netcdf_folder, 'preprocessor', f'*{prp_pattern}*.nc'))
    
        if len(files) == 0:
            raise Exception(f"--Warning: No raw files matching the given patern {pattern} were detected!")

    else:
        raise Exception(f"--Error: The provided mtype parameter ({signal_type}) is not supported. Please select one of: 'raw', 'rangecor'")


    return(files)


def make_plot_filename(fpath_list, options):
    
    smoothing_window = options['smoothing_window']
    timescale = options['timescale']
    normalization_range = options['normalization_range']
    background_range = options['background_range']
    smoothing_range = options['smoothing_range']
    signal_type = options['signal_type']
    
    if len(fpath_list) > 0:
        if timescale == None: avg = ''
        else:
            avg = f"avg_{timescale}"
            
        if smoothing_range == None: sm_range = ''
        else:
            sm_range = f"sm_{int(smoothing_range[0])}_{smoothing_range[1]}"
        
        if smoothing_window == None: sm_win = ''
        else:
            sm_win = f"win_{int(smoothing_window)}"
        
        if normalization_range == None or signal_type == 'raw': norm = ''
        else:
            norm = f"nrm_{int(normalization_range[0])}_{normalization_range[1]}"

        if background_range == None or signal_type == 'rangecor': bline = ''
        else:
            bline = f"bgr_{int(background_range[0])}_{background_range[1]}"

        bname = ['_'.join(filter(None, [f"{os.path.basename(os.path.splitext(path)[0])}", avg, sm_range, sm_win, norm, bline])) for path in fpath_list]
    else:
        bname = []
        
    return(bname)

def get_converter_signals(fpath, options):
    
    channels = options['channels']
    timescale = options['timescale']
    smoothing_window = options['smoothing_window']
    smoothing_range = options['smoothing_range']
    background_range = options['background_range']
    statistics_range = options['statistics_range']
  
    if os.path.exists(fpath):
        system_info, channel_info, time_info, _,\
            sig, _, shots, _ = \
                   short_reader(fpath = fpath, 
                                exclude_telescope_type = '', 
                                exclude_channel_type = '', 
                                exclude_acquisition_mode = '', 
                                exclude_channel_subtype = '', 
                                use_channels = None)

        date_info = dict()
        date_info['start_date'] = system_info.RawData_Start_Date
        date_info['start_time'] = system_info.RawData_Start_Time_UT
        date_info['end_time'] = system_info.RawData_Stop_Time_UT
        
        all_channels = sig.copy().channel.values
        
        channels = check_channels(channels, all_channels)
        
        sig = sig.copy().sel({'channel' : channels})
        channel_info = channel_info.copy().loc[channels,:] 
        
        ranges = bin_to_range(sig = sig,
                              zero_bin = channel_info.DAQ_Trigger_Offset,
                              range_resolution = channel_info.Raw_Data_Range_Resolution)
        
        for ch in channels:
            check_options(options = options, 
                          ranges = ranges.loc[ch,:].copy().values)

        stats = calculate_statistics(sig = sig, 
                                     ranges = ranges, 
                                     statistics_range = statistics_range)
       
        sig = time_averaging(sig = sig, 
                             timescale = timescale)
        
        sig = smoothing(sig = sig, 
                        ranges = ranges, 
                        smoothing_window = smoothing_window, 
                        smoothing_range = smoothing_range)
        
        sig = background_correction(sig = sig, 
                                    ranges = ranges,
                                    background_range = background_range)
                     
    else:
        raise Exception("-- Warning: No raw files found") 
        
    return(sig, ranges, date_info, stats, system_info, channel_info, time_info, shots)

def get_prepro_signals(fpath, options):
            
    channels = options['channels']
    timescale = options['timescale']
    smoothing_window = options['smoothing_window']
    smoothing_range = options['smoothing_range']
    normalization_range = options['normalization_range']
    statistics_range = options['statistics_range']
  
    if os.path.exists(fpath):
        data = xr.open_dataset(fpath)
        
        ranges = data.Range_levels / 1000.
        
        # Extract signal
        sig = data.Range_Corrected_Signals
        sig = sig.copy().where(sig != nc.default_fillvals['f8'])
        
        date_info = dict()
        date_info['start_date'] = data.RawData_Start_Date
        date_info['start_time'] = data.RawData_Start_Time_UT
        date_info['end_time'] = data.RawData_Stop_Time_UT
        
        all_channels = sig.copy().channel.values
        
        channels = check_channels(channels, all_channels)
            
        sig = sig.copy().sel({'channel' : channels})
        ranges = ranges.copy().sel({'channel' : channels})
        
        for ch in channels:
            check_options(options = options, 
                          ranges = ranges.loc[ch,:].copy().values)

        stats = calculate_statistics(sig = sig, 
                                     ranges = ranges, 
                                     statistics_range = statistics_range)
       
        sig = time_averaging(sig = sig, 
                             timescale = timescale)
        
        sig = smoothing(sig = sig, 
                        ranges = ranges, 
                        smoothing_window = smoothing_window, 
                        smoothing_range = smoothing_range)
        
        
        sig = normalize(sig = sig, 
                        ranges = ranges, 
                        normalization_range = normalization_range)
        
                    
    else:
        raise Exception("-- Error: Preprocessed file not found")

    return(sig, ranges, date_info, stats)

def bin_to_range(sig, zero_bin, range_resolution):
    
    channels = sig.copy().channel.values
    bins = sig.copy().bins.values.astype(int)
    
    ranges = xr.DataArray(dims = sig.dims[1:3], coords = [sig.coords[sig.dims[1]], sig.coords[sig.dims[2]]])

    for ch in channels:
        ranges.loc[ch,:] = (bins + zero_bin.loc[ch] + 0.5) * range_resolution.loc[ch] * 1E-3 

    return(ranges)

def time_averaging(sig, timescale):
    
    sig_avg = sig.copy()
    
    if not timescale == None:
        
        if timescale == 'all':
           sig_avg  = sig.copy().mean(dim = 'time', keepdims = True)
        else:
            count = sig.copy().resample(time = timescale).count()
            sig_avg = sig.copy().resample(time = timescale).mean()\
                .where(count >= 0.7 * count.max()).dropna(dim = 'time', how = 'all')
        
    return(sig_avg)

def smoothing(sig, ranges, smoothing_window, smoothing_range):
    
    sig_sm = sig.copy()
    channels = sig.copy().channel.values
    
    if not smoothing_window == None:
        
        for ch in channels:
        
            if smoothing_range == None:
                smoothing_llim = ranges.loc[ch,:].copy()[0].values
            else: smoothing_llim = smoothing_range[0]

            if smoothing_range == None:
                smoothing_ulim =  ranges.loc[ch,:].copy()[-1].values
            else: smoothing_ulim = smoothing_range[1]
                  
            sig_sm_ch, _ = \
                sliding_average_2D_fast(z_vals = sig.loc[:,ch,:].copy().values, 
                                        y_vals = ranges.loc[ch,:].copy().values,
                                        y_sm_lims = [smoothing_llim, smoothing_ulim],
                                        y_sm_win = smoothing_window,
                                        err_type = 'std')
            
            sig_sm.loc[:, ch, :] = sig_sm_ch
    
    return(sig_sm)

def normalize(sig, ranges, normalization_range):
    
    sig_nrm = sig.copy()
    channels = sig.copy().channel.values
    
    if not normalization_range == None: 
        
        for ch in channels:

            normalization_llim = np.where(ranges.loc[ch,:].values >= normalization_range[0])[0][0]
            normalization_ulim = np.where(ranges.loc[ch,:].values <= normalization_range[1])[0][-1]

            sig_div = sig.copy().loc[:, ch, :][:, normalization_llim:normalization_ulim].mean(dim = 'bins')
            sig_mlt = sig_div.copy().mean(dim = 'time')

            sig_nrm.loc[:,ch,:] = (sig.copy().loc[:,ch,:] / sig_div.copy()) * sig_mlt.copy()
    
    return(sig_nrm)

def background_correction(sig, ranges, background_range):
    
    sig_bgc = sig.copy()
    channels = sig.copy().channel.values

    if background_range != None:

        for ch in channels:
            
            background_llim = np.where(ranges.loc[ch,:].values >= background_range[0])[0][0]
            background_ulim = np.where(ranges.loc[ch,:].values <= background_range[1])[0][-1]

            if ch[6] == 'a':
                sig_bgr = sig.copy().loc[:, ch, :][:, background_llim:background_ulim].mean(dim = 'bins')
                sig_bgc.loc[:, ch, :] = sig.copy().loc[:, ch, :] - sig_bgr.copy()
            else:
                sig_bgc.loc[:, ch, :] = sig.copy().loc[:, ch, :]
    
    return(sig_bgc)

def calculate_statistics(sig, ranges, statistics_range):
    
    time = (sig.copy().time - sig.copy().time[0]).dt.seconds.values + \
        1e-6 * (sig.copy().time - sig.copy().time[0]).dt.microseconds.values
    
    channels = sig.copy().channel.values
    
    if not statistics_range == None:
        stats = pd.DataFrame(index = channels, columns = ['mean', 'sdev', 'sem', 'vert_slope', 'temp_slope', 'gaussian_noise', 'profiles', 'bins', 'points'])

        for ch in channels:

            statistics_llim = np.where(ranges.loc[ch,:].values >= statistics_range[0])[0][0]
            statistics_ulim = np.where(ranges.loc[ch,:].values <= statistics_range[1])[0][-1]
                    
            x_vals = ranges.loc[ch,:][statistics_llim:statistics_ulim]
            y_vals = sig.loc[:,ch,:][:, statistics_llim:statistics_ulim].copy().mean(dim = 'time').values
            
            vert_fit = linregress(x = x_vals, y = y_vals)

            xt_vals = time
            yt_vals = sig.loc[:,ch,:][:, statistics_llim:statistics_ulim].copy().mean(dim = 'bins').values
            
            temp_fit = linregress(x = xt_vals, y = yt_vals)
            
            stats.loc[ch, 'profiles'] = time.size
            stats.loc[ch, 'bins'] = x_vals.size
            stats.loc[ch, 'points'] = sig.loc[:,ch,:][:, statistics_llim:statistics_ulim].copy().count().values
            
            stats.loc[ch, 'mean'] = round_it(sig.loc[:,ch,:][:, statistics_llim:statistics_ulim].copy().mean().values, 5)
            stats.loc[ch, 'sdev'] = round_it(sig.loc[:,ch,:][:, statistics_llim:statistics_ulim].copy().std().values, 3)
            stats.loc[ch, 'sem'] = round_it(stats.loc[ch, 'sdev'] / np.sqrt(stats.loc[ch, 'points']), 3)

            stats.loc[ch, 'vert_slope'] = vert_fit[3] <= 0.05
            stats.loc[ch, 'temp_slope'] = temp_fit[3] <= 0.05
            stats.loc[ch, 'gaussian_noise'] = shapiro(y_vals - np.mean(y_vals))[1] > 0.05
            
    else: stats = None
    
    return(stats)

def date_text(start_date, start_time, end_time):
    
    if len(start_date) > 0 and len(start_time) > 0 and len(end_time) > 0:
        
        date = f'{start_date[6:]}.{start_date[4:6]}.{start_date[:4]}'
        
        start = f'{start_time[:2]}:{start_time[2:4]}:{start_time[4:6]}'
    
        end = f'{end_time[:2]}:{end_time[2:4]}:{end_time[4:6]}'
        
        date_part = f'On {date} from {start} to {end} UTC'.strip()

    else:
        date_part = ''
        
    return(date_part)

def make_title_line(iterables):
    
    label = ''
    for part in iterables:
        if len(part) > 0:
            if len(label) == 0:
                label = label + part
            else:
                label = label + ', ' + part
    
    label = str.strip(label)
    
    return(label)

def make_title_mline(channels, options, date_info):

    smoothing_window = options['smoothing_window']
    smoothing_range = options['smoothing_range']
    timescale = options['timescale']
    normalization_range = options['normalization_range']
    background_range = options['background_range']
    signal_type = options['signal_type'].capitalize()
    mtype = options['mtype'].capitalize()
    custom_label = options['custom_label'].capitalize()
    
    date_part = date_text(start_date = date_info['start_date'], 
                          start_time = date_info['start_time'],
                          end_time = date_info['end_time']) 
    
    if not timescale == None: time_label = f'Avg: {timescale}'
    else: time_label = ''

    units = 'km'
        
    if not smoothing_window == None:
        sm_win_label = f"Smoothing Window: {int(smoothing_window)} m"
        
        if not smoothing_range == None:
            sm_label = f"Smoothing Range: {smoothing_range[0]} - {smoothing_range[1]} {units}"
        else: sm_label = ''
    
    else: 
        sm_win_label = ''
        sm_label = ''
        
    if not normalization_range == None:
        nr_label = f"Normalization Range: {normalization_range[0]} - {normalization_range[1]} {units}"
    else: nr_label = ''

    if background_range == None or signal_type == 'rangecor': bl_label = ''
    else: 
            bl_label = f"Background Range: {background_range[0]} - {background_range[1]} {units}"
    
    options_part_1 = make_title_line(iterables = [time_label, sm_win_label, sm_label])
    options_part_2 = make_title_line(iterables = [nr_label, bl_label])

    titles = dict()
    
    for channel in channels:
        ch_label = f'Ch: {channel}'
        
        ch_part = make_title_line(iterables = [custom_label, mtype, signal_type, ch_label])
        
        if len(date_part) > 0:
            titles[channel] = str.strip(f'{ch_part}\n{date_part}\n{options_part_1}\n{options_part_2}')
        
    return(titles)

def multi_line_plot_bokeh(signal, signal_type, ranges = [], channels = None, 
                          zero_bins = [], range_resolution = [], 
                          mtype = '', colorscale = 'sequential',
                          titles = None, plot_folder = None, filename = None,
                          stats = None, statistics_range = None):
    
    if plot_folder != None and os.path.exists(os.path.dirname(os.path.normpath(plot_folder))) == False:
        raise Exception(f"-- The provided plot_folder is not placed in a valid directory:\n{plot_folder}\nPlease use a valid directory to create the plot_folder or use a valid existing path for the plot_folder")
    
    n_time = signal.copy().time.size
    n_channel = signal.copy().channel.size
    
    bins = signal.copy().bins.values
    channels = signal.copy().channel.values
    tabs_obj = np.nan * np.zeros(n_channel, dtype = object)
    
    colors = make_colorscale(colorscale = colorscale, 
                             n_time = n_time) 
    
    
    for j in range(n_channel):
        ch = channels[j]

        ranges_ch = ranges[j,:].copy().values
        
        x_label_1, x_label_2, y_label, y_units = \
            make_labels(signal_type = signal_type, 
                        channel_mode = ch[6])
            
        x_vals_1, x_range_1, x_range_2 = \
            find_x_limits(signal_type = signal_type, 
                          bins = bins, 
                          ranges = ranges_ch)
        
        
        p = figure(width = 1600, 
                   height = 800, 
                   title = titles[ch],
                   x_axis_label = x_label_1, 
                   y_axis_label = y_label,
                   x_range = x_range_1)   

        for i in range(n_time):
            p.line(x_vals_1, signal[i, j, :].values, line_width = 2, color = colors[i])

        if signal_type == 'raw':
            p.extra_x_ranges = {"x2": Range1d(start = x_range_2[0], end = x_range_2[1])}
            p.add_layout(LinearAxis(x_range_name = "x2", axis_label = x_label_2), 'above')

        if not statistics_range == None:
            if signal_type == 'raw':
                box_llim = np.where(ranges_ch >= statistics_range[0])[0][0]
                box_ulim = np.where(ranges_ch <= statistics_range[1])[0][-1]
            else:
                box_llim = statistics_range[0]
                box_ulim = statistics_range[1]               
            
            box = BoxAnnotation(left = box_llim, right = box_ulim, fill_color = 'gray', fill_alpha = 0.3)
            p.add_layout(box)
            
            mean = Label(x = 960, y = 430, text = f"Mean: {stats.loc[ch,'mean']} {y_units}", x_units = 'screen', y_units = 'screen')
            sdev = Label(x = 960, y = 410, text = f"SDev: {stats.loc[ch,'sdev']} {y_units}", x_units = 'screen', y_units = 'screen')
            sem = Label(x = 960, y = 390, text = f"SEM: {stats.loc[ch,'sem']} {y_units}", x_units = 'screen', y_units = 'screen')
            points = Label(x = 960, y = 370, text = f"Points: {stats.loc[ch,'points']}", x_units = 'screen', y_units = 'screen')

            vert_slope = Label(x = 1130, y = 430, text = f"Sign. slope (range): {stats.loc[ch,'vert_slope']}", x_units = 'screen', y_units = 'screen')
            temp_slope = Label(x = 1130, y = 410, text = f"Sign. slope (time): {stats.loc[ch,'temp_slope']}", x_units = 'screen', y_units = 'screen')
            noise = Label(x = 1130, y = 390, text = f"Gaussian noise: {stats.loc[ch,'gaussian_noise']}", x_units = 'screen', y_units = 'screen')
            nums = Label(x = 1130, y = 370, text = f"Bins: {stats.loc[ch,'bins']}, Profiles: {stats.loc[ch,'profiles']}", x_units = 'screen', y_units = 'screen')

            p.add_layout(mean)
            p.add_layout(sdev)
            p.add_layout(sem)
            p.add_layout(points)

            p.add_layout(vert_slope)
            p.add_layout(temp_slope)
            p.add_layout(noise)
            p.add_layout(nums)
        
        p.title.text_font_size = '22pt'    
        p.xaxis.axis_label_text_font_size = '22pt'    
        p.yaxis.axis_label_text_font_size = '22pt'    
        p.xaxis.major_label_text_font_size = '22pt'    
        p.yaxis.major_label_text_font_size = '22pt'  
        
        p.ygrid.minor_grid_line_color = 'black'
        p.ygrid.minor_grid_line_alpha = 0.1
        p.xgrid.minor_grid_line_color = 'black'
        p.xgrid.minor_grid_line_alpha = 0.1
        
        tabs_obj[j] = TabPanel(child = p, title = ch)

    tabs = Tabs(tabs=list(tabs_obj[tabs_obj == tabs_obj]))
    
    if plot_folder == None:
        show(tabs)
    else:
        os.makedirs(plot_folder, exist_ok = True)
        fpath = os.path.join(plot_folder, f'{filename}.html')
        output_file(fpath, mode = 'inline')
        
        save(tabs)
        webbrowser.open_new_tab(fpath) 

    return(tabs)

def check_channels(channels, all_channels):
    
    if channels == None:
        channels = all_channels
    else:
        ch_exists = [ch in all_channels for ch in channels]
        if not np.all(ch_exists):
            raise Exception(f"--Error: The provided channels ({channels}) do not exist. Please use a subset of: {all_channels}")
        else:
            channels = np.array(channels)
        
    return(channels)

def make_colorscale(colorscale, n_time):
    
    if colorscale == 'discreet':
        colors = int(np.ceil(n_time / 10)) * Category10[10]
    elif colorscale == 'sequential':
        if n_time <= 256:
            colors = turbo(n_time)
        else:
            max_colors = turbo(256)
            indexes = np.arange(0, n_time)
            color_index = (np.round(255 * indexes / n_time, decimals = 0)).astype(int)
            colors = tuple([max_colors[ind] for ind in color_index])
                
    return(colors)

def make_labels(signal_type, channel_mode):
    
    if signal_type == 'raw' and channel_mode == 'a':
        y_units = 'mV'
        y_label = f'Raw Signal [{y_units}]'
        x_label_1 = 'Bins'
        x_label_2 = 'Range [km]'
    elif signal_type == 'raw' and channel_mode == 'p':
        y_units = 'Counts'
        y_label = f'Raw Signal [{y_units}]'
        x_label_1 = 'Bins'
        x_label_2 = 'Range [km]'
    elif signal_type == 'rangecor':
        y_units = 'AU'
        y_label = 'RC Signal'
        x_label_1 = 'Range [km]'
        x_label_2 = ''
    else:
        raise Exception(f"-- Error: The provided stype ({signal_type}) was not recognize. Please select one of: ['rangecor','raw']")
    
    return(x_label_1, x_label_2, y_label, y_units)

def find_x_limits(signal_type, bins, ranges):
    
    n_bins = bins.size

    buffer_bins = np.round(0.05 * n_bins, 0)
    range_resolution = ranges[1] - ranges[0]

    if signal_type == 'raw':
        x_vals_1 = bins
        x_llim_1 = x_vals_1[0] - buffer_bins # in bins
        x_ulim_1 = x_vals_1[-1] + buffer_bins # in bins
        x_vals_2 = ranges
        x_llim_2 = x_vals_2[0] - buffer_bins * range_resolution # in km
        x_ulim_2 = x_vals_2[-1] + buffer_bins * range_resolution # in km
        x_range_1 = (x_llim_1, x_ulim_1)
        x_range_2 = (x_llim_2, x_ulim_2)
    elif signal_type == 'rangecor':
        x_vals_1 = ranges
        x_llim_1 = x_vals_1[0] - buffer_bins * range_resolution # in km
        x_ulim_1 = x_vals_1[-1] + buffer_bins * range_resolution # in km
        x_range_1 = (x_llim_1, x_ulim_1)
        x_range_2 = []
    else:
        raise Exception(f"-- Error: The provided stype ({signal_type}) was not recognize. Please select one of: ['rangecor','raw']")
       
    return(x_vals_1, x_range_1, x_range_2)

def round_it(x, sig):
    
    if not np.isfinite(x) or np.isnan(x):
        x = -999.
        sig = 3
        
    if x != 0.:
        x_out = np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)
    else:
        x_out = 0.
        
    return x_out
