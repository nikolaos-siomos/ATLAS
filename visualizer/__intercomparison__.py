#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, os, glob, sys
import xarray as xr
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from .readers.parse_cmp_args import call_parser, check_parser
from .readers.check import check_channels_no_exclude as check_channels
from .plotting import make_axis, make_title, make_plot
from .tools.smoothing import sliding_average_1D
from .tools import normalize
from .tools import normalize, average, differentiate, curve_fit

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
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
    
    # Relative Stnadard error of the mean upper limit for the auto_fit
    rsem_lim = 0.01
    
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

        alt1 = data1.Height_levels.loc[ch1_d].values
        alt2 = data2.Height_levels.loc[ch2_d].values
        
        
        alt_com = np.sort(np.unique(np.hstack([alt1,alt2])))
    
        rng1 = data1.Range_levels.loc[ch1_d].values
        rng2 = data2.Range_levels.loc[ch2_d].values

        rng_com = np.sort(np.unique(np.hstack([rng1,rng2])))
    
        if args['use_range'] == False:
            sig1_ch = xr.DataArray(sig1_ch, dims = ['altitude'],  coords = [alt1])
            sig2_ch = xr.DataArray(sig2_ch, dims = ['altitude'],  coords = [alt2])

            sig1_ch = sig1_ch.interp(altitude = alt_com)
            sig2_ch = sig2_ch.interp(altitude = alt_com)
            
        else:
            sig1_ch = xr.DataArray(sig1_ch, dims = ['range'],  coords = [alt1])
            sig2_ch = xr.DataArray(sig2_ch, dims = ['range'],  coords = [alt2])
           
            sig1_ch = sig1_ch.interp(range = rng_com)
            sig2_ch = sig2_ch.interp(range = rng_com)
    
        sig1_ch = sig1_ch.values
        sig2_ch = sig2_ch.values

        # Create the y axis (height/range)
        x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
            make_axis.rayleigh_x(heights = alt_com, 
                                 ranges = rng_com,
                                 x_lims = args['x_lims'], 
                                 use_dis = args['use_range'])
    
        # Smoothing
        if args['smooth']:
            if not isinstance(args['smoothing_window'],list):
                from .tools.smoothing import sliding_average_1D_fast as smooth_1D
            else:
                from .tools.smoothing import sliding_average_1D as smooth_1D

            y1_vals_sm, y1_errs = \
                smooth_1D(y_vals = sig1_ch.copy(), 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')
                
            y2_vals_sm, y2_errs = \
                smooth_1D(y_vals = sig2_ch.copy(), 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')      
        
        else:
            y1_vals_sm = sig1_ch.copy()
            y2_vals_sm = sig2_ch.copy()
            y1_errs = np.nan * y1_vals_sm
            y2_errs = np.nan * y2_vals_sm
    
        rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef = \
            curve_fit.stats(y1 = sig1_ch.copy(),
                            y2 = sig2_ch.copy(), 
                            x  = x_vals,
                            min_win = 1.,
                            max_win = 4.,
                            step = 0.1,
                            llim = 2.,
                            ulim = 16.,
                            rsem_lim = rsem_lim,
                            cross_check_all_points = False,
                            cross_check_type = 'all')

            
        # Negate the cross-check mask for the signal comparison
        # mneg = np.ones(mder.shape, dtype = bool)
        
        norm_region, idx, ismol = \
            curve_fit.scan(mfit = mfit,
                           dflt_region = args['normalization_region'],
                           auto_fit = args['auto_fit'])
                
        coef_c = float(coef[idx].values)
        rsem_c = float(rsem[idx].values)
        nder_c = float(nder[idx].values)
        mder_c = float(mder[idx].values)
    
        # Create the y axis (signal)
        y_llim, y_ulim, y_label = \
            make_axis.intercomparison_y(sig1 = coef_c * y1_vals_sm[slice(x_lbin,x_ubin+1)], 
                                        sig2 = y2_vals_sm[slice(x_lbin,x_ubin+1)], 
                                        y_lims = args['y_lims'],
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
        fname = f'{data1.Measurement_ID}_{data2.Measurement_ID}_{data1.Lidar_Name}_{data2.Lidar_Name}_cmp_{channels1[j]}_{channels2[j]}_ATLAS_{__version__}.png'

        # Make the plot
        fpath = make_plot.rayleigh(dir_out = args['output_folder'], 
                                   fname = fname, title = title,
                                   dpi_val = args['dpi'],
                                   color_reduction = args['color_reduction'],
                                   use_lin = args['use_lin_scale'],
                                   norm_region = norm_region,
                                   x_vals = x_vals, 
                                   y1_vals = y1_vals_sm,
                                   y2_vals = y2_vals_sm,
                                   y1_errs = y1_errs,
                                   y2_errs = y2_errs,
                                   coef = coef_c,
                                   rsem = rsem_c,
                                   rslope = nder_c,
                                   pval = mder_c,
                                   rsem_lim = rsem_lim,
                                   x_lbin = x_lbin, x_ubin = x_ubin,
                                   x_llim = x_llim, x_ulim = x_ulim, 
                                   y_llim = y_llim, y_ulim = y_ulim, 
                                   x_label = x_label, y_label = y_label,
                                   x_tick = args['x_tick'],
                                   label_1 = channels1[j],
                                   label_2 = channels2[j])  
        
        if args['auto_fit']:
            fname_mask = f'{data1.Measurement_ID}_{data2.Measurement_ID}_{data1.Lidar_Name}_{data2.Lidar_Name}_cmp_{channels1[j]}_{channels2[j]}_mask_ATLAS_{__version__}.png'
    
            make_plot.rayleigh_mask(dir_out = args['output_folder'], 
                                    fname = fname_mask, title = title,
                                    dpi_val = args['dpi'],
                                    color_reduction = args['color_reduction'],
                                    mfit = mfit, mder = mder, msec = msec, 
                                    mshp = mshp, mcrc = mcrc, rsem = rsem,
                                    rsem_lim = rsem_lim)
        
        # Add metadata to the quicklook plot
        from PIL import Image
        from PIL import PngImagePlugin
       
        METADATA = {"processing_software_1" : f"ATLAS_{data1.version}",
                    "processing_software_2" : f"ATLAS_{data2.version}",
                    "measurement_id_1" : f"{data1.Measurement_ID}",
                    "measurement_id_2" : f"{data2.Measurement_ID}",
                    "channel_1" : f"{ch1_d}",
                    "channel_2" : f"{ch2_d}",
                    "smooth" : f"{args['smooth']}",
                    "smoothing_exponential" : f"{args['smooth_exponential']}",
                    "smoothing_range" : f"{args['smoothing_range']}",
                    "smoothing_window": f"{args['smoothing_window']}",
                    "dpi" : f"{args['dpi']}",
                    "color_reduction" : f"{args['color_reduction']}",
                    "normalization_region" : f"{args['normalization_region']}",
                    "use_lin_scale" : f"{args['use_lin_scale']}",
                    "use_range" : f"{args['use_range']}",
                    "x_lims" : f"{args['x_lims']}",
                    "y_lims" : f"{args['y_lims']}",
                    "x_tick" : f"{args['x_tick']}"}

                
        im = Image.open(fpath)
        meta = PngImagePlugin.PngInfo()
    
        for x in METADATA.keys():
            meta.add_text(x, METADATA[x])
            
        im.save(fpath, "png", pnginfo = meta)
        

    print('-----------------------------------------')
    print(' ')
    return()

if __name__ == '__main__':
    
    sys.path.append('../')
    
    from version import __version__
    
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args, __version__)
    
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
    #             "use_range" : f"{args['use_range']}",
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
