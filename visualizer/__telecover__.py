#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, os, sys, glob
import xarray as xr
import numpy as np
import netCDF4 as nc
from .readers.parse_tlc_args import call_parser, check_parser
from .readers.check import check_channels
from .plotting import make_axis, make_title, make_plot
from .tools.smoothing import sliding_average_2D, sliding_average_1D
from .tools import normalize

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print('Initializing the Telecover Test...')
    print('-----------------------------------------')
    
    # Read the quicklook file
    data = xr.open_dataset(args['input_file'])
    ranges = data.Range_levels
    
    # Extract signal
    if 'Range_Corrected_Signals_North_Sector' in data.keys():
        sig_n = data['Range_Corrected_Signals_North_Sector']
        sig_e = data['Range_Corrected_Signals_East_Sector']
        sig_s = data['Range_Corrected_Signals_South_Sector']
        sig_w = data['Range_Corrected_Signals_West_Sector']
        
        sig_n = sig_n.copy().where(sig_n != nc.default_fillvals['f8'])
        sig_e = sig_e.copy().where(sig_e != nc.default_fillvals['f8'])
        sig_s = sig_s.copy().where(sig_s != nc.default_fillvals['f8'])
        sig_w = sig_w.copy().where(sig_w != nc.default_fillvals['f8'])
        
        if args['use_non_rangecor'] == True:
            sig_n = sig_n.copy() / (ranges.copy() ** 2.)
            sig_e = sig_e.copy() / (ranges.copy() ** 2.)
            sig_s = sig_s.copy() / (ranges.copy() ** 2.)
            sig_w = sig_w.copy() / (ranges.copy() ** 2.)
    
    else:
        sig_n = []
        sig_e = []
        sig_s = []
        sig_w = []
        
    if 'Range_Corrected_Signals_Outer_Ring' in data.keys():
        sig_o = data['Range_Corrected_Signals_Outer_Ring']  
        sig_i = data['Range_Corrected_Signals_Inner_Ring']
        
        sig_o = sig_o.copy().where(sig_o != nc.default_fillvals['f8'])
        sig_i = sig_i.copy().where(sig_i != nc.default_fillvals['f8'])
        
        if args['use_non_rangecor'] == True:
            sig_o = sig_o.copy() / (ranges.copy() ** 2.)
            sig_i = sig_i.copy() / (ranges.copy() ** 2.)    
    
    else:
        sig_o = []
        sig_i = []
    
    # Check if the parsed channels exist
    channels = \
        check_channels(sel_channels = args['channels'], 
                       all_channels = data.channel.values,
                       exclude_field_type = args['exclude_field_type'], 
                       exclude_scattering_type = args['exclude_scattering_type'], 
                       exclude_detection_mode = args['exclude_detection_mode'], 
                       exclude_channel_subtype = args['exclude_channel_subtype'])
    
    # iterate over the channels
    for ch in channels:
        print(f"-- channel: {ch}")
        
        ch_d = dict(channel = ch)
        
        # Create the y axis (height/range)
        x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
            make_axis.rayleigh_x(heights = data.Height_levels.loc[ch_d].values, 
                                 ranges = data.Range_levels.loc[ch_d].values,
                                 x_lims = args['x_lims'], 
                                 use_dis = args['use_distance'])
    
        if isinstance(sig_n,list) == False:
            iters = np.min([sig_n.time_n.size,
                            sig_e.time_e.size,
                            sig_s.time_s.size,
                            sig_w.time_w.size])
    
            y_n = sig_n.loc[ch_d].values
            y_e = sig_e.loc[ch_d].values
            y_s = sig_s.loc[ch_d].values
            y_w = sig_w.loc[ch_d].values
    
            # Store the extra iteration if it exists 
            if y_n.shape[0] == iters + 1:
                y_n2 = y_n[iters,:]
            else:
                y_n2 = []
    
            if y_e.shape[0] == iters + 1:
                y_e2 = y_e[iters,:]
            else:
                y_e2 = []
            
            if y_s.shape[0] == iters + 1:
                y_s2 = y_s[iters,:]
            else:
                y_s2 = []
                
            if y_w.shape[0] == iters + 1:
                y_w2 = y_w[iters,:]
            else:
                y_w2 = []
                
            # Averaging sector signals
            y_m_n = np.nanmean(y_n[:iters,:], axis = 0) 
            y_m_e = np.nanmean(y_e[:iters,:], axis = 0) 
            y_m_s = np.nanmean(y_s[:iters,:], axis = 0) 
            y_m_w = np.nanmean(y_w[:iters,:], axis = 0)     
    
            # Smoothing
            if args['smooth']:
                # Smoothing averaged sectors
                y_m_sm_n, _ = \
                    sliding_average_1D(y_vals = y_m_n, 
                                       x_vals = x_vals,
                                       x_sm_lims = args['smoothing_range'],
                                       x_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
                y_m_sm_e, _ = \
                    sliding_average_1D(y_vals = y_m_e, 
                                       x_vals = x_vals,
                                       x_sm_lims = args['smoothing_range'],
                                       x_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
                    
                y_m_sm_s, _ = \
                    sliding_average_1D(y_vals = y_m_s, 
                                       x_vals = x_vals,
                                       x_sm_lims = args['smoothing_range'],
                                       x_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
                y_m_sm_w, _ = \
                    sliding_average_1D(y_vals = y_m_w, 
                                       x_vals = x_vals,
                                       x_sm_lims = args['smoothing_range'],
                                       x_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
                # Smoothing the extra sector
                if len(y_n2) > 0:
                    y_sm_n2, _ = \
                        sliding_average_1D(y_vals = y_n2, 
                                           x_vals = x_vals,
                                           x_sm_lims = args['smoothing_range'],
                                           x_sm_hwin = args['half_window'],
                                           expo = args['smooth_exponential'])
                else:
                    y_sm_n2 = []
                    
                if len(y_e2) > 0:
                    y_sm_e2, _ = \
                        sliding_average_1D(y_vals = y_e2, 
                                           x_vals = x_vals,
                                           x_sm_lims = args['smoothing_range'],
                                           x_sm_hwin = args['half_window'],
                                           expo = args['smooth_exponential'])
                else:
                    y_sm_e2 = []
                    
                if len(y_s2) > 0:
                    y_sm_s2, _ = \
                        sliding_average_1D(y_vals = y_s2, 
                                           x_vals = x_vals,
                                           x_sm_lims = args['smoothing_range'],
                                           x_sm_hwin = args['half_window'],
                                           expo = args['smooth_exponential'])
                else:
                    y_sm_s2 = []
                    
                if len(y_w2) > 0:
                    y_sm_w2, _ = \
                        sliding_average_1D(y_vals = y_w2, 
                                           x_vals = x_vals,
                                           x_sm_lims = args['smoothing_range'],
                                           x_sm_hwin = args['half_window'],
                                           expo = args['smooth_exponential'])
                else:
                    y_sm_w2 = []
                    
                # Smoothing unaveraged sectors
                y_sm_n = \
                    sliding_average_2D(z_vals = y_n.T, 
                                       y_vals = x_vals,
                                       y_sm_lims = args['smoothing_range'],
                                       y_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
                y_sm_e = \
                    sliding_average_2D(z_vals = y_e.T, 
                                       y_vals = x_vals,
                                       y_sm_lims = args['smoothing_range'],
                                       y_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
                    
                y_sm_s = \
                    sliding_average_2D(z_vals = y_s.T, 
                                       y_vals = x_vals,
                                       y_sm_lims = args['smoothing_range'],
                                       y_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
                y_sm_w = \
                    sliding_average_2D(z_vals = y_w.T, 
                                       y_vals = x_vals,
                                       y_sm_lims = args['smoothing_range'],
                                       y_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
            
            else:
                y_m_sm_n = y_m_n
                y_m_sm_e = y_m_e
                y_m_sm_s = y_m_s
                y_m_sm_w = y_m_w
                
                if len(y_n2) > 0: y_sm_n2 = y_n2
                else: y_sm_n2 = []
               
                if len(y_e2) > 0: y_sm_e2 = y_e2
                else: y_sm_e2 = []
    
                if len(y_s2) > 0: y_sm_s2 = y_s2
                else: y_sm_s2 = []
                
                if len(y_w2) > 0: y_sm_w2 = y_w2
                else: y_sm_w2 = []
                
                y_sm_n = y_n.T
                y_sm_e = y_e.T
                y_sm_s = y_s.T
                y_sm_w = y_w.T
                
    
            y_l_sm_n = np.nanmin(y_sm_n[:,:iters], axis = 1)
            y_l_sm_e = np.nanmin(y_sm_e[:,:iters], axis = 1)
            y_l_sm_s = np.nanmin(y_sm_s[:,:iters], axis = 1)
            y_l_sm_w = np.nanmin(y_sm_w[:,:iters], axis = 1)
    
            y_u_sm_n = np.nanmax(y_sm_n[:,:iters], axis = 1)
            y_u_sm_e = np.nanmax(y_sm_e[:,:iters], axis = 1)
            y_u_sm_s = np.nanmax(y_sm_s[:,:iters], axis = 1)
            y_u_sm_w = np.nanmax(y_sm_w[:,:iters], axis = 1)
            
            # Normalization for smoothed signals
            n_coef, n_bin = normalize.to_a_point(sig = y_m_sm_n, 
                                                 sig_b = np.ones(y_m_sm_n.shape), 
                                                 x_vals = x_vals,
                                                 norm = args['normalization_height'],
                                                 hwin = args['half_normalization_window'],
                                                 axis = 0)
        
            e_coef, _ = normalize.to_a_point(sig = y_m_sm_e, 
                                             sig_b = np.ones(y_m_sm_e.shape), 
                                             x_vals = x_vals,
                                             norm = args['normalization_height'],
                                             hwin = args['half_normalization_window'],
                                             axis = 0)
        
            s_coef, _ = normalize.to_a_point(sig = y_m_sm_s, 
                                             sig_b = np.ones(y_m_sm_s.shape), 
                                             x_vals = x_vals,
                                             norm = args['normalization_height'],
                                             hwin = args['half_normalization_window'],
                                             axis = 0)
        
            w_coef, _ = normalize.to_a_point(sig = y_m_sm_w, 
                                             sig_b = np.ones(y_m_sm_w.shape), 
                                             x_vals = x_vals,
                                             norm = args['normalization_height'],
                                             hwin = args['half_normalization_window'],
                                             axis = 0)
        
            # Create the y axis (signal)
            y_llim, y_ulim, y_llim_nr, y_ulim_nr = \
                make_axis.telecover_y(sig = [y_m_sm_n[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_e[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_s[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_w[slice(x_lbin,x_ubin+1)]],
                                      sig_nr = [n_coef * y_m_sm_n[slice(x_lbin,x_ubin+1)],
                                                e_coef * y_m_sm_e[slice(x_lbin,x_ubin+1)],
                                                s_coef * y_m_sm_s[slice(x_lbin,x_ubin+1)],
                                                w_coef * y_m_sm_w[slice(x_lbin,x_ubin+1)]],
                                      y_lims = args['y_lims'])
            
                    
            # Make title
            title = make_title.telecover(start_date = data.RawData_Start_Date,
                                        start_time = data.RawData_Start_Time_UT, 
                                        end_time = data.RawData_Stop_Time_UT, 
                                        lidar = data.Lidar_Name, 
                                        channel = ch, 
                                        zan = data.Laser_Pointing_Angle,
                                        loc = data.Lidar_Location,
                                        sm_lims = args['smoothing_range'],
                                        sm_hwin = args['half_window'],
                                        sm_expo = args['smooth_exponential'])
        
        
            # Make filename
            fname = f'{data.Measurement_ID}_{data.Lidar_Name}_tlc_sectors_{ch}_ATLAS_{__version__}.png'
        
            # Make the plot
            fpath = \
                make_plot.telecover_sec(dir_out = args['output_folder'], 
                                        fname = fname, title = title,
                                        dpi_val = args['dpi'],
                                        x_refr = args['normalization_height'],
                                        refr_hwin = args['half_normalization_window'],
                                        x_vals = x_vals, 
                                        y1_vals = y_m_sm_n, 
                                        y2_vals = y_m_sm_e, 
                                        y3_vals = y_m_sm_s, 
                                        y4_vals = y_m_sm_w,
                                        y1_extr = y_sm_n2, 
                                        y2_extr = y_sm_e2, 
                                        y3_extr = y_sm_s2, 
                                        y4_extr = y_sm_w2, 
                                        y1_lvar = y_l_sm_n, 
                                        y2_lvar = y_l_sm_e, 
                                        y3_lvar = y_l_sm_s, 
                                        y4_lvar = y_l_sm_w,
                                        y1_uvar = y_u_sm_n, 
                                        y2_uvar = y_u_sm_e, 
                                        y3_uvar = y_u_sm_s, 
                                        y4_uvar = y_u_sm_w,
                                        coef_1 = n_coef, 
                                        coef_2 = e_coef, 
                                        coef_3 = s_coef, 
                                        coef_4 = w_coef,
                                        x_lbin = x_lbin, x_ubin = x_ubin,
                                        x_llim = x_llim, x_ulim = x_ulim, 
                                        y_llim = y_llim, y_ulim = y_ulim, 
                                        y_llim_nr = y_llim_nr, 
                                        y_ulim_nr = y_ulim_nr, 
                                        x_label = x_label,
                                        x_tick = args['x_tick'])  
    
        if isinstance(sig_o,list) == False:
            iters = np.min([sig_o.time_o.size,
                            sig_i.time_i.size])
    
            y_o = sig_o.loc[ch_d].values
            y_i = sig_i.loc[ch_d].values
                
            # Averaging sector signals
            y_m_o = np.nanmean(y_o[:iters,:], axis = 0) 
            y_m_i = np.nanmean(y_i[:iters,:], axis = 0) 
    
            # Smoothing
            if args['smooth']:
                # Smoothing averaged sectors
                y_m_sm_o, _ = \
                    sliding_average_1D(y_vals = y_m_o, 
                                       x_vals = x_vals,
                                       x_sm_lims = args['smoothing_range'],
                                       x_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
                y_m_sm_i, _ = \
                    sliding_average_1D(y_vals = y_m_i, 
                                       x_vals = x_vals,
                                       x_sm_lims = args['smoothing_range'],
                                       x_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
                    
                # Smoothing unaveraged sectors
                y_sm_o = \
                    sliding_average_2D(z_vals = y_o.T, 
                                       y_vals = x_vals,
                                       y_sm_lims = args['smoothing_range'],
                                       y_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
                y_sm_i = \
                    sliding_average_2D(z_vals = y_i.T, 
                                       y_vals = x_vals,
                                       y_sm_lims = args['smoothing_range'],
                                       y_sm_hwin = args['half_window'],
                                       expo = args['smooth_exponential'])
    
            else:
                y_m_sm_o = y_m_o
                y_m_sm_i = y_m_i
                
                y_sm_o = y_o.T
                y_sm_i = y_i.T
                    
            # Calculate the atmospheric variability limits
            y_l_sm_o = np.nanmin(y_sm_o[:,:iters], axis = 1)
            y_l_sm_i = np.nanmin(y_sm_i[:,:iters], axis = 1)
    
            y_u_sm_o = np.nanmax(y_sm_o[:,:iters], axis = 1)
            y_u_sm_i = np.nanmax(y_sm_i[:,:iters], axis = 1)
            
            # Normalization for smoothed signals
            o_coef, n_bin = normalize.to_a_point(sig = y_m_sm_o, 
                                                 sig_b = np.ones(y_m_sm_o.shape), 
                                                 x_vals = x_vals,
                                                 norm = args['normalization_height'],
                                                 hwin = args['half_normalization_window'],
                                                 axis = 0)
        
            i_coef, _ = normalize.to_a_point(sig = y_m_sm_i, 
                                             sig_b = np.ones(y_m_sm_i.shape), 
                                             x_vals = x_vals,
                                             norm = args['normalization_height'],
                                             hwin = args['half_normalization_window'],
                                             axis = 0)
        
            # Create the y axis (signal)
            y_llim, y_ulim, y_llim_nr, y_ulim_nr = \
                make_axis.telecover_y(sig = [y_m_sm_o[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_i[slice(x_lbin,x_ubin+1)]],
                                      sig_nr = [o_coef * y_m_sm_o[slice(x_lbin,x_ubin+1)],
                                                i_coef * y_m_sm_i[slice(x_lbin,x_ubin+1)]],
                                      y_lims = args['y_lims'])
            
                    
            # Make title
            title = make_title.telecover(start_date = data.RawData_Start_Date,
                                         start_time = data.RawData_Start_Time_UT, 
                                         end_time = data.RawData_Stop_Time_UT, 
                                         lidar = data.Lidar_Name, 
                                         channel = ch, 
                                         zan = data.Laser_Pointing_Angle,
                                         loc = data.Lidar_Location,
                                         sm_lims = args['smoothing_range'],
                                         sm_hwin = args['half_window'],
                                         sm_expo = args['smooth_exponential'])
        
            # Make filename
            fname = f'{data.Measurement_ID}_{data.Lidar_Name}_tlc_rings_{ch}_ATLAS_{__version__}.png'
        
            # Make the plot
            fpath = \
                make_plot.telecover_sec(dir_out = args['output_folder'], 
                                        fname = fname, title = title,
                                        dpi_val = args['dpi'],
                                        x_refr = args['normalization_height'],
                                        refr_hwin = args['half_normalization_window'],
                                        x_vals = x_vals, 
                                        y1_vals = y_m_sm_o, 
                                        y2_vals = y_m_sm_i, 
                                        y1_lvar = y_l_sm_o, 
                                        y2_lvar = y_l_sm_i, 
                                        y1_uvar = y_u_sm_o, 
                                        y2_uvar = y_u_sm_i, 
                                        coef_1 = o_coef, 
                                        coef_2 = i_coef, 
                                        x_lbin = x_lbin, x_ubin = x_ubin,
                                        x_llim = x_llim, x_ulim = x_ulim, 
                                        y_llim = y_llim, y_ulim = y_ulim, 
                                        y_llim_nr = y_llim_nr, 
                                        y_ulim_nr = y_ulim_nr, 
                                        x_label = x_label,
                                        x_tick = args['x_tick']) 
    
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
