#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, glob, os
import xarray as xr
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from .readers.parse_pcl_args import call_parser, check_parser
from .readers.check import check_channels
from .tools import normalize
from .plotting import make_axis, make_title, make_plot
from .tools.smoothing import sliding_average_1D
from .tools import average

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print('Initializing the Pol. Calibration...')
    print('-----------------------------------------')
    
    # Read the pol cal file
    data = xr.open_dataset(args['input_file'])
    
    # Extract signal
    sig_ray = data.Range_Corrected_Signals_Rayleigh
    sig_ray = sig_ray.copy().where(sig_ray != nc.default_fillvals['f8'])
    
    sig_m45 = data.Range_Corrected_Signals_minus_45
    sig_m45 = sig_m45.copy().where(sig_m45 != nc.default_fillvals['f8'])\
        .mean(dim = 'time_m45')
    
    sig_p45 = data.Range_Corrected_Signals_plus_45
    sig_p45 = sig_p45.copy().where(sig_p45 != nc.default_fillvals['f8'])\
        .mean(dim = 'time_p45')
    
    # Extract signal time, channels, and bins
    channels = data.channel.values
    channels_r = args['ch_r']
    channels_t = args['ch_t']
    
    # Extract the K value
    if args['K'] == None:
        K = len(channels_r) * [1.]
    else:
        K = args['K']
    
    # Check if the parsed channels exist
    channels_r = check_channels(sel_channels = channels_r, all_channels = channels)
    channels_t = check_channels(sel_channels = channels_t, all_channels = channels)
    
    # Extract Molecular Depolarization Ratio and Calucalte the Atm. Parameter alpha
    mldr = data.Molecular_Linear_Depolarization_Ratio
    a_m = (1. - mldr) / (1. + mldr)
    
    # Iterate over the channels
    for ch_r, ch_t, K_ch in zip(channels_r, channels_t, K):
        
        print(f"-- channels: {ch_r} & ch_t")

        ch_r_d = dict(channel = ch_r)
        ch_t_d = dict(channel = ch_t)
        
        sig_r_p45_ch = sig_p45.loc[ch_r_d].values
        sig_t_p45_ch = sig_p45.loc[ch_t_d].values
        sig_r_m45_ch = sig_m45.loc[ch_r_d].values
        sig_t_m45_ch = sig_m45.loc[ch_t_d].values
        sig_r_ray_ch = sig_ray.loc[ch_r_d].values
        sig_t_ray_ch = sig_ray.loc[ch_t_d].values
        
        a_m_ch = a_m.loc[ch_r_d].values
    
        # Create the y axis (height/range)
        y_lbin_cal, y_ubin_cal, y_llim_cal, y_ulim_cal, y_vals_cal, y_label_cal = \
            make_axis.polarization_calibration_y(
                heights = data.Height_levels_Calibration.loc[ch_r_d].values, 
                ranges = data.Range_levels_Calibration.loc[ch_r_d].values,
                y_lims = args['y_lims_calibration'], 
                use_dis = args['use_distance'])
    
        # Create the y axis (height/range)
        y_lbin_ray, y_ubin_ray, y_llim_ray, y_ulim_ray, y_vals_ray, y_label_ray = \
            make_axis.polarization_calibration_y(
                heights = data.Height_levels_Rayleigh.loc[ch_r_d].values, 
                ranges = data.Range_levels_Rayleigh.loc[ch_r_d].values,
                y_lims = args['y_lims_rayleigh'], 
                use_dis = args['use_distance'])
    
        # Smoothing
        if args['smooth']:
            # Smoothing averaged sectors
            x_r_m45, _ = \
                sliding_average_1D(y_vals = sig_r_m45_ch, 
                                   x_vals = y_vals_cal,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])
    
            x_t_m45, _ = \
                sliding_average_1D(y_vals = sig_t_m45_ch, 
                                   x_vals = y_vals_cal,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])
                
            x_r_p45, _ = \
                sliding_average_1D(y_vals = sig_r_p45_ch, 
                                   x_vals = y_vals_cal,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])
    
            x_t_p45, _ = \
                sliding_average_1D(y_vals = sig_t_p45_ch, 
                                   x_vals = y_vals_cal,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])
    
            x_r_ray, _ = \
                sliding_average_1D(y_vals = sig_r_ray_ch, 
                                   x_vals = y_vals_ray,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])    
                
            x_t_ray, _ = \
                sliding_average_1D(y_vals = sig_t_ray_ch, 
                                   x_vals = y_vals_ray,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])  
        
        else:
            x_r_m45 = sig_r_m45_ch
            x_t_m45 = sig_t_m45_ch
            x_r_p45 = sig_r_p45_ch
            x_t_p45 = sig_t_p45_ch
            x_r_ray = sig_r_ray_ch
            x_t_ray = sig_t_ray_ch
            
        
        avg_r_m45 = average.region(sig = x_r_m45, 
                                   x_vals = y_vals_cal, 
                                   calibr = args['calibration_height'], 
                                   hwin = args['half_calibration_window'], 
                                   axis = 0,
                                   squeeze = True)
        
        avg_t_m45 = average.region(sig = x_t_m45, 
                                   x_vals = y_vals_cal, 
                                   calibr = args['calibration_height'], 
                                   hwin = args['half_calibration_window'], 
                                   axis = 0,
                                   squeeze = True)
        
        avg_r_p45 = average.region(sig = x_r_p45, 
                                   x_vals = y_vals_cal, 
                                   calibr = args['calibration_height'], 
                                   hwin = args['half_calibration_window'], 
                                   axis = 0,
                                   squeeze = True)
        
        avg_t_p45 = average.region(sig = x_t_p45, 
                                   x_vals = y_vals_cal, 
                                   calibr = args['calibration_height'], 
                                   hwin = args['half_calibration_window'], 
                                   axis = 0,
                                   squeeze = True)
        
        avg_r_ray = average.region(sig = x_r_ray, 
                                   x_vals = y_vals_ray, 
                                   calibr = args['calibration_height'], 
                                   hwin = args['half_calibration_window'], 
                                   axis = 0,
                                   squeeze = True)
        
        avg_t_ray = average.region(sig = x_t_ray, 
                                   x_vals = y_vals_ray, 
                                   calibr = args['calibration_height'], 
                                   hwin = args['half_calibration_window'], 
                                   axis = 0,
                                   squeeze = True)
    
        eta_m45_prf = (x_r_m45 / x_t_m45)
    
        eta_p45_prf = (x_r_p45 / x_t_p45)
        
        eta_prf = np.sqrt(x_r_p45 * x_r_m45) / K_ch
        
        eta_m45 = (avg_r_m45 / avg_t_m45)
    
        eta_p45 = (avg_r_p45 / avg_t_p45)
        
        eta = np.sqrt(eta_p45 * eta_m45) / K_ch
        
        delta_s_prf = (x_r_ray / x_t_ray) #/ eta
        
        delta_s = (avg_r_ray / avg_t_ray) #/ eta
        
        psi = (eta_p45 - eta_m45) / (eta_p45 + eta_m45)
        
        kappa = 1.
        
        epsilon = np.rad2deg(0.5 * np.arcsin(np.tan(0.5 * np.arcsin(psi) / kappa)))
        
        # kappa = np.tan(0.5 * np.arcsin(psi)) / np.sin(2. * np.deg2rad(epsilon)) 
    
            
        # Create the x axis (calibration)
        x_llim_cal, x_ulim_cal = \
            make_axis.polarization_calibration_cal_x(
                ratio_m = eta_m45_prf[slice(y_lbin_cal,y_ubin_cal+1)], 
                ratio_p = eta_p45_prf[slice(y_lbin_cal,y_ubin_cal+1)],
                x_lims_cal = args['x_lims_calibration'])
            
        # Create the x axis (rayleigh)
        x_llim_ray, x_ulim_ray = \
            make_axis.polarization_calibration_ray_x(
                ratio = delta_s_prf[slice(y_lbin_ray,y_ubin_ray+1)],
                x_lims_ray = args['x_lims_rayleigh'])
        
                
        # Make title
        title_cal = make_title.polarization_calibration(
            start_date = data.RawData_Start_Date_Calibration,
            start_time = data.RawData_Start_Time_UT_Calibration, 
            end_time = data.RawData_Stop_Time_UT_Calibration, 
            lidar = data.Lidar_Name_Calibration, 
            channel_r = ch_r, 
            channel_t = ch_t, 
            zan = data.Laser_Pointing_Angle_Calibration,
            lat = data.Latitude_degrees_north_Calibration, 
            lon = data.Longitude_degrees_east_Calibration, 
            elv = data.Altitude_meter_asl_Calibration)
        
        title_ray = make_title.polarization_calibration(
            start_date = data.RawData_Start_Date_Rayleigh,
            start_time = data.RawData_Start_Time_UT_Rayleigh, 
            end_time = data.RawData_Stop_Time_UT_Rayleigh, 
            lidar = data.Lidar_Name_Rayleigh, 
            channel_r = ch_r, 
            channel_t = ch_t, 
            zan = data.Laser_Pointing_Angle_Rayleigh,
            lat = data.Latitude_degrees_north_Rayleigh, 
            lon = data.Longitude_degrees_east_Rayleigh, 
            elv = data.Altitude_meter_asl_Rayleigh)
       
    
        # Make filename
        fname = f'pcl_{data.Measurement_ID_Calibration}_{ch_r}_to_{ch_t}.png'
    
        raise Exception("Halted!")
    
        # Make the plot# Read the quicklook file
    data = xr.open_dataset(args['input_file'])
    
    # Delete all existing png files within
    pngs = glob.glob(os.path.join(args['output_folder'],'*.png'))
    for file in pngs:
        os.remove(file)
    
    # Extract signal
    sig = data.Range_Corrected_Signals
    sig = sig.copy().where(sig != nc.default_fillvals['f8'])
    
    # Extract Attenuated Backscatter
    atb = data.Attenuated_Backscatter
    
    # Extract signal time, channels, and bins
    time = sig.time.values
    bins = sig.bins.values
    
    # Check if the parsed channels exist
    channels = check_channels(sel_channels = args['channels'], 
                              all_channels = data.channel.values)
    
    # iterate over the channels
    for ch in channels:
        
        ch_d = dict(channel = ch)
        sig_ch = sig.loc[ch_d].values
        atb_ch = atb.loc[ch_d].values
    
        # Create the y axis (height/range)
        x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
            make_axis.rayleigh_x(heights = data.Height_levels.loc[ch_d].values, 
                                 ranges = data.Range_levels.loc[ch_d].values,
                                 x_lims = args['x_lims'], 
                                 use_dis = args['use_distance'])
    
        # Smoothing
        if args['smooth']:
            y_vals_sm, y_errs = \
                sliding_average_1D(y_vals = sig_ch.copy(), 
                                   x_vals = x_vals,
                                   x_sm_lims = args['smoothing_range'],
                                   x_sm_hwin = args['half_window'],
                                   expo = args['smooth_exponential'])
    
        else:
            y_vals_sm = sig_ch.copy()
            y_errs = np.nan * y_vals_sm
    
        # Normalization of the signals
        coef = normalize.to_a_point(sig = y_vals_sm, 
                                    sig_b = atb_ch.copy(), 
                                    x_vals = x_vals,
                                    norm = args['reference_height'],
                                    hwin = args['half_reference_window'],
                                    axis = 0)
            
        # Create the y axis (signal)
        y_llim, y_ulim, y_label = \
            make_axis.rayleigh_y(sig = coef * y_vals_sm[slice(x_lbin,x_ubin+1)], 
                                 atb = atb_ch.copy(), 
                                 y_lims = args['y_lims'] , 
                                 use_lin = args['use_lin_scale'])
        
                
        # Make title
        title = make_title.rayleigh(start_date = data.RawData_Start_Date,
                                    start_time = data.RawData_Start_Time_UT, 
                                    end_time = data.RawData_Stop_Time_UT, 
                                    lidar = data.Lidar_Name, 
                                    channel = ch, 
                                    zan = data.Laser_Pointing_Angle,
                                    lat = data.Latitude_degrees_north, 
                                    lon = data.Longitude_degrees_east, 
                                    elv = data.Altitude_meter_asl)
    
        # Make filename
        fname = f'ray_{data.Measurement_ID}_{ch}.png'
    
        # Make the plot
        fpath = make_plot.rayleigh(dir_out = args['output_folder'], 
                                   fname = fname, title = title,
                                   dpi_val = args['dpi'],
                                   use_lin = args['use_lin_scale'],
                                   x_refr = args['reference_height'],
                                   refr_hwin = args['half_reference_window'],
                                   x_vals = x_vals, 
                                   y1_vals = y_vals_sm,
                                   y2_vals = atb_ch.copy(),
                                   y1_errs = y_errs,
                                   coef = coef,
                                   x_lbin = x_lbin, x_ubin = x_ubin,
                                   x_llim = x_llim, x_ulim = x_ulim, 
                                   y_llim = y_llim, y_ulim = y_ulim, 
                                   x_label = x_label, y_label = y_label,
                                   x_tick = args['x_tick']) 
        
        fpath = make_plot.polarization_calibration(
            dir_out = args['output_folder'], 
            fname = fname, title_1 = title_cal, title_2 = title_ray, 
            dpi_val = args['dpi'],
            x_refr = args['calibration_height'],
            refr_hwin = args['half_calibration_window'],
            y21_vals = delta_s, y22_vals = mldr,   
            y11_vals = eta_m45_prf, y12_vals = eta_p45_prf,
            y1_vals = y_vals_cal, y2_vals = y_vals_ray, 
            y_tick = args['y_tick'])  
    
    print('-----------------------------------------')
    print(' ')
    
    return()

if __name__ == '__main__':
    # Get the command line argument information
    args = call_parser()
    
    # Check the command line argument information
    args = check_parser(args)
    
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
