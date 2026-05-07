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
from visualizer.tools.smoothing import sliding_average_2D_bin_fast, sliding_average_2D_fast, sliding_average_1D_fast
from scipy.stats import linregress, shapiro
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import TabPanel, Tabs, LinearAxis, Range1d, BoxAnnotation, Label
from bokeh.palettes import Category10, turbo
import webbrowser 
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from PIL import Image
import numpy as np
from pathlib import Path
from visualizer.readers.read_prepro import unpack
from helper_functions.time_conversions import iso_to_datetimes

# https://docs.bokeh.org/en/latest/docs/reference/palettes.html#bokeh-palettes

# allowed_mtypes = ['ray', 'drk']

def check_options(options, heights = None):
    
    options["channels"] = None
    
    none_allowed = dict(plot_folder = False,
                        timescale = True,
                        smoothing_window = True,
                        normalization_range = True,
                        # mtype = False,
                        colorscale = False,
                        statistics_range = True,
                        custom_label = False,
                        use_molecular = False,
                        start_date = True,
                        slice_time = True)
    
    option_type = dict(plot_folder = [str],
                       timescale = [str],
                       smoothing_window = [int, float],
                       normalization_range = [int, float],
                       # mtype = [str],
                       colorscale = [str],
                       statistics_range = [int, float],
                       custom_label = [str],
                       use_molecular = [bool],
                       start_date = [str],
                       slice_time = [str])
    
    option_shape = dict(plot_folder = 'scalar',
                        timescale = 'scalar',
                        smoothing_window = 'scalar',
                        normalization_range = 'list',
                        # mtype = 'scalar',
                        colorscale = 'scalar',
                        statistics_range = 'list',
                        custom_label = 'scalar',
                        use_molecular = 'scalar',
                        start_date = 'scalar',
                        slice_time = 'list')
    
    
    check_none(options = options, none_allowed = none_allowed)
    
    check_list(options = options, option_shape = option_shape)
    
    check_option_type(options = options, option_type = option_type)    
    
    check_parent_directory(options  = options, key = 'plot_folder')
                    
    colorbar = ['sequential', 'discrete']
    check_allowed_values(options = options,
                         # keys = ['mtype', 'colorscale'],
                         keys = ['colorscale'],
                         allowed_values = [colorbar])
                         
    if not type(heights) == type(None):
        check_smoothing_window(options = options, 
                               max_range = heights[-1])
        
        for key in ['smoothing_range', 'normalization_range', 
                    'statistics_range']:
            check_range_limits(options = options, 
                               key = key, 
                               min_range = heights[0], 
                               max_range = heights[-1])
    
    return()


def hhmm_to_datetime(hhmm: str, base: pd.Timestamp) -> pd.Timestamp:
    """Convert hhmm string to datetime on the date of base timestamp."""
    return pd.Timestamp(year=base.year,
                        month=base.month,
                        day=base.day,
                        hour=int(hhmm[:2]),
                        minute=int(hhmm[2:]))

def make_interval(start_str: str, stop_str: str, base: pd.Timestamp):
    start_dt = hhmm_to_datetime(start_str, base)
    stop_dt = hhmm_to_datetime(stop_str, base)

    # If stop is "before" start → assume it's the next day
    if stop_dt <= start_dt:
        stop_dt += pd.Timedelta(days=1)

    return start_dt, stop_dt

def check_channel_lists(channels_1, channels_2):

    channel_list_checks(channels_1, "channels_1")
    channel_list_checks(channels_2, "channels_2")
    
    if len(channels_1) != len(channels_2):
        raise Exception(f"Lists channels_1 and channels_2 don't have the same length.\nChannels_1: {channels_1}\nChannels_2: {channels_2}")
        

def channel_list_checks(channels, identifier):
    
    if not isinstance(channels, list):
        raise Exception(f"Parameter {identifier} is not a list. Please provide a list of valid ATLAS channel IDS")
    else:
        if len(channels) == 0:
            raise Exception(f"Parameter {identifier} is empty. Please provide a list of valid ATLAS channel IDS")
    
    for ch in channels:
        if type(ch) != str:
            raise Exception(f"Wrong value format in {identifier} list. ATLAS channel ID must be a string. Provided: {ch}")
    
    
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
    

def check_allowed_values(options, keys, allowed_values):

    for key, values in zip(keys, allowed_values):
        if options[key] not in values:
            raise Exception(f"--Error: The provided {key} parameter is not supported:\n{options[key]}\nPlease select one of: {values}")
    
    return()

def check_path(options, key):
    
    p = Path(options[key])
    
    if p.parent.parent.parent.exists() and p.parent.parent.parent.is_dir():
        p.mkdir(exist_ok=True)
        
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

# def get_fpaths(netcdf_folder, mtype): 
    
#     if mtype not in allowed_mtypes:
#         raise Exception(f"--Error: Provided mtype {mtype} not supported. Allowed mtypes: {allowed_mtypes}")
        
#     prp_pattern = f"_qck_{mtype}_ATLAS_"
#     files = glob.glob(os.path.join(netcdf_folder, 'preprocessor', f'*{prp_pattern}*.nc'))
#     # exclude_pattern = f"_qck_{mtype}_ATLAS_"
#     # exclude_files = glob.glob(os.path.join(netcdf_folder, 'preprocessor', f'*{exclude_pattern}*.nc'))
#     # files = [file for file in files if file not in exclude_files]
    
#     if len(files) == 0:
#         raise Exception(f"--Warning: No preprocessed files matching the given patern {prp_pattern} were detected!")

#     return(files)


def make_plot_filename(fpath_list, options):
    
    smoothing_window = options['smoothing_window']
    timescale = options['timescale']
    normalization_range = options['normalization_range']
    smoothing_range = options['smoothing_range']
    
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
        
        if normalization_range == None: norm = ''
        else:
            norm = f"nrm_{int(normalization_range[0])}_{normalization_range[1]}"

        bname = ['_'.join(filter(None, [f"{os.path.basename(os.path.splitext(path)[0])}", avg, sm_range, sm_win, norm])) for path in fpath_list]
    else:
        bname = []
        
    return(bname)

def get_converter_signals(fpath, options, channels):
    
    timescale = options['timescale']
    smoothing_window = options['smoothing_window']
    smoothing_range = options['smoothing_range']
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
        
        heights = bin_to_range(sig = sig,
                              zero_bin = channel_info.DAQ_Trigger_Offset,
                              range_resolution = channel_info.Raw_Data_Range_Resolution)
        
        for ch in channels:
            check_options(options = options, 
                          heights = heights.loc[ch,:].copy().values)

        stats = calculate_statistics(sig = sig, 
                                     heights = heights, 
                                     statistics_range = statistics_range)
       
        sig = time_averaging(sig = sig, 
                             timescale = timescale)
        
        sig = smoothing(sig = sig, 
                        heights = heights, 
                        smoothing_window = smoothing_window, 
                        smoothing_range = smoothing_range)

    else:
        raise Exception("-- Warning: No raw files found") 
        
    return(sig, heights, date_info, stats, system_info, channel_info, time_info, shots)

def get_prepro_signals(fpath, options, channels):
            
    smoothing_window = options['smoothing_window']
    smoothing_range = options['smoothing_range']
    normalization_range = options['normalization_range']
    statistics_range = options['statistics_range']
    use_molecular = options['use_molecular']
    start_date = options['start_date']
    slice_time = options['slice_time']
  
    if os.path.exists(fpath):
        profiles, metadata = unpack(fpath)
        # data = xr.open_dataset(fpath)
        
        # if "heights" in profiles.keys():
        heights = profiles["heights"] / 1000.
        # elif "heights_cal" in profiles.keys():
        #     heights = profiles["heights_cal"] / 1000.

        # sig = []
        # atb = []
        
        # if "sig" in profiles.keys():
        # Extract signal
        sig = profiles["sig"]
        sig = sig.copy().where(sig != nc.default_fillvals['f8'])
        
        atb = profiles["atb"]
        atb = atb.copy().where(atb != nc.default_fillvals['f8'])
        # elif "sig_pcb" in profiles.keys():
        #     sig_p45 = profiles["sig_p45"]
        #     sig_m45 = profiles["sig_m45"]
        #     sig_p45 = sig_p45.copy().where(sig != nc.default_fillvals['f8'])            
        #     sig_m45 = sig_p45.copy().where(sig != nc.default_fillvals['f8']) 
        #     sig = xr.concat([sig_p45, sig_m45], dim = "time").sortby("time")
        # elif 
        
        
        # date_info = dict()
        # date_info['start_date'] = data.RawData_Start_Date
        # date_info['start_time'] = data.RawData_Start_Time_UT
        # date_info['end_time'] = data.RawData_Stop_Time_UT
        
        all_channels = sig.copy().channel.values
        channels = check_channels(channels, all_channels)
            
        sig = sig.copy().sel({'channel' : channels})
        atb = atb.copy().sel({'channel' : channels})
        heights = heights.copy().sel({'channel' : channels})
        
        for ch in channels:
            check_options(options = options, 
                          heights = heights.loc[ch,:].copy().values)

        stats = calculate_statistics(sig = sig, 
                                     heights = heights, 
                                     statistics_range = statistics_range)
       
        # start_times = iso_to_datetimes(profiles['sig'].time)

        if slice_time != None:
            mask_time = None
    
            if start_date == None:    
                base = pd.to_datetime(metadata['start_date'] + metadata['start_time'], utc=True)
            
            else:
                base = pd.to_datetime(start_date, utc=True)
                
    
            for i in range(0, len(slice_time), 2):
                start_t = slice_time[i]
                stop_t = slice_time[i+1]
                
                time = sig.time
                
                start_dt, stop_dt =  make_interval(
                    start_str = start_t, 
                    stop_str = stop_t, 
                    base = base
                    )
                
                mask_interval = (time >= start_dt) & (time <= stop_dt)
                
                if mask_time is None:
                    mask_time = mask_interval
                else:
                    mask_time = (mask_time | mask_interval)                    
            
            sig = sig.where(mask_time)
            times = sig.time[mask_time].values
        else:

            times = sig.time.values

        # Extract the datetime value
        t = pd.to_datetime(times)

        # Format into strings
        dates_str = t.strftime("%Y%m%d") #t.strftime("%d.%m.%Y")
        times_str = t.strftime("%H%M") #t.strftime("%H:%M")
        
        metadata["start_date"] = dates_str[0]
        metadata["start_time"] = times_str[0]
        
        metadata["stop_date"] = dates_str[-1]
        metadata["stop_time"] = times_str[-1]

        
        # dates_str_fname = t.strftime("%Y%m%d")
        # times_str_fname = t.strftime("%H%M")
        
        # start_date = dates_str[0]
        # start_time = times_str[0]
        
        # end_date = dates_str[-1]
        # end_time = times_str[-1]
        
        sig = sig.mean("time")

        
        sig, sig_err = smoothing(sig = sig, 
                                 heights = heights, 
                                 smoothing_window = smoothing_window, 
                                 smoothing_range = smoothing_range)

        atb,_ = smoothing(sig = atb, 
                          heights = heights, 
                          smoothing_window = smoothing_window, 
                          smoothing_range = smoothing_range)
        
        sig, sig_norm_val = normalize(sig = sig, 
                                      heights = heights, 
                                      normalization_range = normalization_range)
        
        sig_err = sig_err.copy() / sig_norm_val
        
        if not use_molecular:
        #     _, atb_norm_val = normalize(sig = atb, 
        #                                 heights = heights, 
        #                                 normalization_range = normalization_range)
            
        #     sig = sig * atb_norm_val
        
        # else:
            atb = xr.full_like(atb, np.nan)     
        
                    
    else:
        raise Exception("-- Error: Preprocessed file not found")

    return(sig, sig_err, atb, heights, metadata, stats)

def bin_to_range(sig, zero_bin, range_resolution):
    
    channels = sig.copy().channel.values
    bins = sig.copy().bins.values.astype(int)
    
    heights = xr.DataArray(dims = sig.dims[1:3], coords = [sig.coords[sig.dims[1]], sig.coords[sig.dims[2]]])

    for ch in channels:
        heights.loc[ch,:] = (bins + zero_bin.loc[ch] + 0.5) * range_resolution.loc[ch] * 1E-3 

    return(heights)

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

def smoothing(sig, heights, smoothing_window, smoothing_range):
    
    sig_sm = sig.copy()
    sig_err = np.nan * sig.copy()
    channels = sig.copy().channel.values
    
    if not smoothing_window == None:
        
        for ch in channels:
            ch_d = {"channel" : ch}
            sig_ch = sig.sel(ch_d).values
            heights_ch = heights.sel(ch_d).values
            
            if smoothing_range == None:
                smoothing_llim = heights.loc[ch,:].copy()[0].values
            else: smoothing_llim = smoothing_range[0]

            if smoothing_range == None:
                smoothing_ulim =  heights.loc[ch,:].copy()[-1].values
            else: smoothing_ulim = smoothing_range[1]
                  
            if "time" in sig.dims:
                sig_sm_ch, sig_err_ch = \
                    sliding_average_2D_fast(z_vals = sig_ch, 
                                            y_vals = heights_ch,
                                            y_sm_lims = [smoothing_llim, smoothing_ulim],
                                            y_sm_win = smoothing_window,
                                            err_type = 'std')
            else:
                sig_sm_ch, sig_err_ch = sliding_average_1D_fast(y_vals = sig_ch,
                                                       x_vals = heights_ch,
                                                       x_sm_lims = [smoothing_llim, smoothing_ulim],
                                                       x_sm_win = smoothing_window,
                                                       err_type = 'std')
                
            
            sig_sm.loc[ch_d] = sig_sm_ch
            sig_err.loc[ch_d] = sig_err_ch
    
    return(sig_sm, sig_err)

def normalize(sig, heights, normalization_range):
    
    sig_nrm = sig.copy()
    
    nrm = xr.ones_like(sig.channel, dtype = float)
    # channels = sig.copy().channel.values
    
    if normalization_range != None: 
                
        mask_heights = (heights >= normalization_range[0]) & (heights <= normalization_range[1])
        
        if "time" in sig.dims:
            nrm = sig.where(mask_heights).mean("time").mean("bins")
        else:
            nrm = sig.where(mask_heights).mean("bins")

        sig_nrm = sig / nrm

    return(sig_nrm, nrm)

def calculate_statistics(sig, heights, statistics_range):
    
    time = (sig.copy().time - sig.copy().time[0]).dt.seconds.values + \
        1e-6 * (sig.copy().time - sig.copy().time[0]).dt.microseconds.values
    
    channels = sig.copy().channel.values
    
    if not statistics_range == None:
        stats = pd.DataFrame(index = channels, columns = ['mean', 'sdev', 'sem', 'vert_slope', 'temp_slope', 'gaussian_noise', 'profiles', 'bins', 'points'])

        for ch in channels:

            statistics_llim = np.where(heights.loc[ch,:].values >= statistics_range[0])[0][0]
            statistics_ulim = np.where(heights.loc[ch,:].values <= statistics_range[1])[0][-1]
                    
            x_vals = heights.loc[ch,:][statistics_llim:statistics_ulim]
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

def round_it(x, sig):
    
    if not np.isfinite(x) or np.isnan(x):
        x = -999.
        sig = 3
        
    if x != 0.:
        x_out = np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)
    else:
        x_out = 0.
        
    return x_out


def plot_intercomparison(
        dir_out, fname, title, dpi_val, color_reduction, use_lin, norm_region,
        X1, Y1, X2, Y2, Y1E, Y2E,
        XC, YC1, YC2, YC1E, YC2E, YCM,
        xlims, x_axis_label, y_axis_label, label_1, label_2
        ):

    
    if  np.isnan(YCM).all() == False:
        norm = np.nanmean(YCM[(XC >= norm_region[0]) & (XC <= norm_region[1])])

        Y1 = norm * Y1
        Y2 = norm * Y2
        Y1E = norm * Y1E
        Y2E = norm * Y2E
        
    else:
        norm = 1.
    
    N1 = np.mean(Y1[(X1 >= norm_region[0]) & (X1 <= norm_region[1])])
    N2 = np.mean(Y2[(X2 >= norm_region[0]) & (X2 <= norm_region[1])])
        
    export_pack = {"height": XC,
                   "signal_1": YC1,
                   "signal_2": YC2,
                   "signal_1_sdev": YC1E,
                   "signal_2_sdev": YC2E,
                   "molecular": YCM}
    
    # Create the figure
    fig = plt.figure(figsize=(15 , 3.5))

    fig.suptitle(title)

    ax = fig.add_axes([0.05,0.145,0.52,0.63])
        
    ax.plot(X1, Y1, color = 'tab:blue', label = label_1)
    ax.plot(X2, Y2, color = 'tab:red', label = label_2)

    if  np.isnan(YCM).all() == False:
        ax.plot(XC, YCM, color = 'tab:green', label = "molecular")
        
    if np.isnan(Y1E).all() == False:
        ax.fill_between(X1, Y1 - Y1E, Y1 + Y1E, color = 'tab:blue', alpha = 0.3)

    if np.isnan(Y2E).all() == False:
        ax.fill_between(X2, Y2 - Y2E, Y2 + Y2E, color = 'tab:red', alpha = 0.3)

    x_tick = choose_tick_step(xlims[1])
    
    x_ticks = np.arange(x_tick * np.ceil(xlims[0] / x_tick), 
                        x_tick * (np.floor(xlims[1] / x_tick) + 1.), 
                        x_tick)

        
    if np.abs(xlims[0] - x_ticks[0]) < x_tick * 0.25:
        x_ticks[0] = xlims[0]
    else:
        x_ticks = np.hstack((xlims[0], x_ticks))

    if np.abs(xlims[1] - x_ticks[-1]) < x_tick * 0.25:
        x_ticks[-1] = xlims[1]
    else:
        x_ticks = np.hstack((x_ticks, xlims[1]))

    x_ticks = np.round(x_ticks, decimals = 3)
    
    max_y = 2 * max(np.nanmax(Y1), np.nanmax(Y2))
    min_y = min(N1,N2) * np.exp(-0.2*(xlims[1]-np.mean(norm_region)))
    ylims = [min_y, max_y]

    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([xlims[0], xlims[1]])
    ax.set_xlabel(x_axis_label)
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.))

    ax.set_ylim(ylims)
    ax.set_ylabel(y_axis_label)
    if use_lin == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    
    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc = 'upper right')
        
    ax.axvspan(norm_region[0], 
               norm_region[1], 
               alpha = 0.2, facecolor = 'tab:grey')

    ax.text(0.45, 0.9, 
            f'norm region: {norm_region[0]} - {norm_region[1]} km',
            transform = ax.transAxes,
            bbox = dict(facecolor = 'tab:blue', alpha = 0.22, zorder = 3))
    # ax.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')


    ax2 = fig.add_axes([0.625,0.145,0.36,0.63])
    
    YCE = np.sqrt(np.power(YC1E,2) + np.power(YC2E,2))
    
    # if np.isnan(YC1E).all() == False and np.isnan(YC2E).all() == True:
    #     ax2.fill_between(XC, (YC1 - YC1E - YC2) / YC2, 
    #                      (YC1 + YC1E - YC2) / YC2, color = 'tab:blue', 
    #                      alpha = 0.3, label = 'sem')
    #     ax2.plot(XC, (YC1 - YC2) / YC2, color = 'tab:blue',label = 'mean')

    # if np.isnan(YC1E).all() == False and np.isnan(YC2E).all() == False:
    #     ax2.fill_between(XC, (YC1 - YC1E - YC2) / YC2, 
    #                      (YC1 + YC1E - YC2) / YC2, color = 'tab:blue', 
    #                      alpha = 0.3, label = 'sem1')
    #     ax2.plot(XC, (YC1 - YC2) / YC2, color = 'tab:blue',label = 'mean')
    #     ax2.fill_between(XC, (YC1 - YC2E - YC2) / YC2, 
    #                      (YC1 + YC2E - YC2) / YC2, color = 'tab:red', 
    #                      alpha = 0.3, label = 'sem2')

    if np.isnan(YCE).all() == False and np.isnan(YCE).all() == True:
        ax2.fill_between(XC, (YC1 - YC2 - YCE) / YC1, 
                         (YC1 - YC2 + YCE) / YC1, color = 'tab:blue', 
                         alpha = 0.3, label = 'stdev')
        ax2.plot(XC, (YC1 - YC2) / YC2, color = 'tab:blue',label = 'mean')

    if np.isnan(YCE).all() == False and np.isnan(YCE).all() == False:
        ax2.fill_between(XC, (YC1 - YC2 - YCE) / YC1, 
                         (YC1 - YC2 + YCE) / YC1, color = 'tab:blue', 
                         alpha = 0.3, label = 'stdev')
        ax2.plot(XC, (YC1 - YC2) / YC1, color = 'tab:blue',label = 'mean')
        
        if ax2.get_legend_handles_labels() != ([], []):
            ax2.legend(loc = 'lower left')

    else:
        ax2.plot(XC, (YC1 - YC2) / YC2, color = 'tab:blue')
        
    
    ax2.axhline(c = 'k')

    x_tick_2 = 2. * x_tick 
    x_ticks_2 = np.arange(x_tick_2 * np.floor(xlims[0] / x_tick_2), 
                          x_tick_2 * (np.ceil(xlims[1] / x_tick_2) + 1.), 
                          x_tick_2)
    
    ax2.set_xticks(x_ticks_2, labels = x_ticks_2)
    ax2.set_xlim([xlims[0], xlims[1]])
    ax2.set_xlabel(x_axis_label)
    ax2.xaxis.set_minor_locator(MultipleLocator(x_tick_2 / 2.))

    y_ticks = np.round(np.arange(-0.10, 0.10 + 0.025, 0.025), decimals = 2)
    ax2.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax2.set_ylim([y_ticks[0], y_ticks[-1]])
    ax2.set_ylabel('Relative Diff. ')
    
    ax2.grid(which = 'both')
    
    # ax2.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')
    
    fpath = os.path.join(dir_out, fname)
            
    os.makedirs(os.path.dirname(fpath), exist_ok = True)
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    perform_color_reduction(color_reduction, fpath)
    
    return(fpath, export_pack)

def perform_color_reduction(color_reduction, plot_path):
    
    if color_reduction == True:
        im = Image.open(plot_path)
        im = im.convert('P', palette = Image.ADAPTIVE, colors = 255) 
        im.save(plot_path)

    return()

def choose_tick_step(xmax_km, target=10):
    """
    Return a 'nice' major tick step in km for a 0..xmax_km axis.
    Prefers round km steps; for small ranges, allows 250, 100 m steps.
    """
    # Add finer increments (0.1 km = 100 m)
    candidates_km = [30, 20, 10, 5, 2, 1, 0.5, 0.25, 0.1]

    best = None
    best_score = None

    for step in candidates_km:
        # Skip overly fine steps for large ranges (>1 km)
        if xmax_km > 1 and step < 0.25:
            continue

        n_ticks = int(np.floor(xmax_km / step)) + 1
        diff = abs(n_ticks - target)
        penalty = 0 if 8 <= n_ticks <= 12 else 1
        score = (penalty, diff, -step)  # smaller = better

        if best_score is None or score < best_score:
            best_score = score
            best = step

    return best