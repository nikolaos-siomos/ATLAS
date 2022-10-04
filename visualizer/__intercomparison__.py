#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, os, glob
import xarray as xr
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from .readers.parse_cmp_args import call_parser, check_parser
from .readers.check import check_channels_no_exclude as check_channels
from .plotting import make_axis, make_title, make_plot
from .tools.smoothing import sliding_average_1D
from .tools import normalize

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print('Initializing the Lidar Intercomparison...')
    print('-----------------------------------------')
    
    # Read the quicklook file
    data1 = xr.open_dataset(args['input_files'][0])
    data2 = xr.open_dataset(args['input_files'][1])
    
    # Extract signal
    sig1 = data1.Range_Corrected_Signals
    sig1 = sig1.copy().where(sig1 != nc.default_fillvals['f8'])
    
    sig2 = data2.Range_Corrected_Signals
    sig2 = sig2.copy().where(sig2 != nc.default_fillvals['f8'])
    
    # Extract Attenuated Backscatter
    atb1 = data1.Attenuated_Backscatter
    atb2 = data2.Attenuated_Backscatter
    
    # Extract signal time, channels, and bins
    time1 = sig1.time.values
    bins1 = sig1.bins.values
    time2 = sig2.time.values
    bins2 = sig2.bins.values
    
    # Check if the parsed channels exist
    channels1 = check_channels(sel_channels = args['channels_1'], 
                               all_channels = data1.channel.values)
    channels2 = check_channels(sel_channels = args['channels_2'], 
                               all_channels = data2.channel.values)
    
    # iterate over the channels
    for j in range(channels1.size):
        print(f"-- ch: {channels1[j]} & {channels2[j]}")

        
        ch1_d = dict(channel = channels1[j])
        ch2_d = dict(channel = channels2[j])
    
        sig1_ch = sig1.loc[ch1_d]
        sig2_ch = sig2.loc[ch2_d]
        
        atb1_ch = atb1.loc[ch1_d]
        atb2_ch = atb2.loc[ch2_d]
        
        alt1 = data1.Height_levels.loc[ch1_d].values
        alt2 = data2.Height_levels.loc[ch2_d].values
        
        
        alt_com = np.sort(np.unique(np.hstack([alt1,alt2])))
    
        rng1 = data1.Range_levels.loc[ch1_d].values
        rng2 = data2.Range_levels.loc[ch2_d].values

        rng_com = np.sort(np.unique(np.hstack([rng1,rng2])))
    
        if args['use_distance'] == False:
            sig1_ch = xr.DataArray(sig1_ch, dims = ['altitude'],  coords = [alt1])
            sig2_ch = xr.DataArray(sig2_ch, dims = ['altitude'],  coords = [alt2])
    
            atb1_ch = xr.DataArray(atb1_ch, dims = ['altitude'],  coords = [alt1])
            atb2_ch = xr.DataArray(atb2_ch, dims = ['altitude'],  coords = [alt2])
                    
            sig1_ch = sig1_ch.interp(altitude = alt_com)
            sig2_ch = sig2_ch.interp(altitude = alt_com)
    
            atb1_ch = atb1_ch.interp(altitude = alt_com)
            atb2_ch = atb2_ch.interp(altitude = alt_com)
            
        else:
            sig1_ch = xr.DataArray(sig1_ch, dims = ['range'],  coords = [alt1])
            sig2_ch = xr.DataArray(sig2_ch, dims = ['range'],  coords = [alt2])
    
            atb1_ch = xr.DataArray(atb1_ch, dims = ['range'],  coords = [alt1])
            atb2_ch = xr.DataArray(atb2_ch, dims = ['range'],  coords = [alt2])
                    
            sig1_ch = sig1_ch.interp(range = rng_com)
            sig2_ch = sig2_ch.interp(range = rng_com)
    
            atb1_ch = atb1_ch.interp(range = rng_com)
            atb2_ch = atb2_ch.interp(range = rng_com)            
        
        sig1_ch = sig1_ch.values
        sig2_ch = sig2_ch.values

        atb1_ch = atb1_ch.values
        atb2_ch = atb2_ch.values  
        
        # Create the y axis (height/range)
        x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
            make_axis.rayleigh_x(heights = alt_com, 
                                 ranges = rng_com,
                                 x_lims = args['x_lims'], 
                                 use_dis = args['use_distance'])
    
        # Smoothing
        if args['smooth']:
            y1_vals_sm, y1_errs = \
                sliding_average_1D(y_vals = sig1_ch.copy(), 
                                   x_vals = x_vals,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])
    
            y2_vals_sm, y2_errs = \
                sliding_average_1D(y_vals = sig2_ch.copy(), 
                                   x_vals = x_vals,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])
        
        else:
            y1_vals_sm = sig1_ch.copy()
            y2_vals_sm = sig2_ch.copy()
            y1_errs = np.nan * y1_vals_sm
            y2_errs = np.nan * y2_vals_sm
    
        # Normalization of the signals        
        coef1 = normalize.to_a_point(sig = y1_vals_sm, 
                                     sig_b = atb1_ch, 
                                     x_vals = x_vals,
                                     norm = args['normalization_height'],
                                     hwin = args['half_normalization_window'],
                                     axis = 0)
        
        coef2 = normalize.to_a_point(sig = y2_vals_sm, 
                                     sig_b = atb2_ch, 
                                     x_vals = x_vals,
                                     norm = args['normalization_height'],
                                     hwin = args['half_normalization_window'],
                                     axis = 0)        
    
        # Create the y axis (signal)
        y_llim, y_ulim, y_label = \
            make_axis.intercomparison_y(sig1 = coef1 * y1_vals_sm[slice(x_lbin,x_ubin+1)],
                                        sig2 = coef2 * y2_vals_sm[slice(x_lbin,x_ubin+1)] , 
                                        y_lims = args['y_lims'] , 
                                        use_lin = args['use_lin_scale'])    
                
        # Make title
        title = \
            make_title.intercomparison(start_date = data2.RawData_Start_Date,
                                       start_time = data2.RawData_Start_Time_UT, 
                                       end_time = data2.RawData_Stop_Time_UT, 
                                       lidar_1 = data1.Lidar_Name,
                                       lidar_2 = data2.Lidar_Name, 
                                       channel_1 = channels1[j],
                                       channel_2 = channels2[j], 
                                       zan = data2.Laser_Pointing_Angle,
                                       lat = data2.Latitude_degrees_north, 
                                       lon = data2.Longitude_degrees_east, 
                                       elv = data2.Altitude_meter_asl)
    
        # Make filename
        fname = f'cmp_{data1.Measurement_ID}_{data2.Measurement_ID}_{channels1[j]}_{channels2[j]}.png'

        # Make the plot
        fpath = make_plot.intercomparison(dir_out = args['output_folder'], 
                                          fname = fname, title = title,
                                          dpi_val = args['dpi'],
                                          use_lin = args['use_lin_scale'],
                                          x_refr = args['normalization_height'],
                                          refr_hwin = args['half_normalization_window'],
                                          x_vals = x_vals, 
                                          y1_vals = y1_vals_sm, 
                                          y2_vals = y2_vals_sm, 
                                          y1_errs = y1_errs, 
                                          y2_errs = y2_errs,
                                          y3_vals = atb1_ch,
                                          coef1 = coef1,
                                          coef2 = coef2,
                                          use_molecular = args['use_molecular'],
                                          x_lbin = x_lbin, x_ubin = x_ubin,
                                          x_llim = x_llim, x_ulim = x_ulim, 
                                          y_llim = y_llim, y_ulim = y_ulim, 
                                          x_label = x_label, y_label = y_label,
                                          x_tick = args['x_tick'],
                                          lidars = [data1.Lidar_Name, data2.Lidar_Name])  

    print('-----------------------------------------')
    print(' ')
    return()

if __name__ == '__main__':
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args)
    
    # sys.exit()
    # # Add metadata to the quicklook plot
    # from PIL import Image
    # from PIL import PngImagePlugin
   
    # METADATA = {"processing_software" : f"ATLAS_{data.version}",
    #             "measurement_id" : f"{data.Measurement_ID}",
    #             "channel" : f"{ch}",
    #             "smooth" : f"{args['smooth']}",
    #             "smoothing_exponential" : f"{args['smooth_exponential']}",
    #             "smoothing_range (lower)" : f"{args['smoothing_range'][0]}",
    #             "smoothing_range (upper)" : f"{args['smoothing_range'][-1]}",
    #             "half_window (lower)": f"{args['half_window'][0]}",
    #             "half_window (upper)": f"{args['half_window'][-1]}",
    #             "dpi" : f"{args['dpi']}",
    #             "use_log_scale" : f"{args['use_log_scale']}",
    #             "use_distance" : f"{args['use_distance']}",
    #             "x_lims (lower)" : f"{x_llim}",
    #             "x_lims (upper)" : f"{x_ulim}",
    #             "y_lims (lower)" : f"{y_vals[y_llim]}",
    #             "y_lims (upper)" : f"{y_vals[y_ulim]}",
    #             "z_lims (lower)" : f"{z_llim}",
    #             "z_lims (upper)" : f"{z_ulim}",
    #             "x_tick" : f"{x_tick}",
    #             "y_tick" : f"{args['y_tick']}"}
            
    # im = Image.open(fpath)
    # meta = PngImagePlugin.PngInfo()

    # for x in METADATA.keys():
    #     meta.add_text(x, METADATA[x])
        
    # im.save(fpath, "png", pnginfo = meta)
