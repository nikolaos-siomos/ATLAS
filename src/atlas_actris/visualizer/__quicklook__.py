#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import pandas as pd
import warnings, os
from visualizer.readers.check import check_channels
from helper_functions.printouts import print_header
from visualizer.plotting import make_axis
from visualizer.make_text import GenerateText, Libraries
from visualizer.plotting.plot_utils import perform_color_reduction
from visualizer.plotting.plot_quicklook import generate_plot

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def clean_quicklook_plots(plot_dir,pattern):
    if not os.path.isdir(plot_dir):
        return

    for filename in os.listdir(plot_dir):
        if "pattern" not in filename:
            continue

        file_path = os.path.join(plot_dir, filename)

        if os.path.isfile(file_path):
            os.remove(file_path)
            
def generate_quicklooks(data_pack, caller_info, settings):
    
    process_qck = caller_info['process_qck']
    
    for key in process_qck:
        if key in data_pack['time_info'].keys():
    
            print_header(f"Start generating {key} quicklooks...")
            
            plot_dir = os.path.join(caller_info["output_folder"], "plots")
            os.makedirs(plot_dir, exist_ok=True)
            
            clean_quicklook_plots(plot_dir, pattern = f"_qck_{key}")  

            profiles = data_pack['profile'][key]

            ranges = data_pack['range'][key]
            height_asl = data_pack['height_asl'][key]

            system_info = data_pack['system_info'][key]
            channel_info = data_pack['channel_info'][key]
            time_info = data_pack['time_info'][key]
                
            # # Extract signal time, channels, and bins
            time = profiles.time.values
            start_time = pd.to_datetime(time_info.sel(parameters="start_time").values)
            
            end_time = pd.to_datetime(time_info.sel(parameters="end_time").values)

            delta_t = (end_time - start_time).total_seconds()
            
            # Check if the parsed channels exist
            channels = \
                check_channels(sel_channels = settings['channels'], 
                               all_channels = profiles.channel.values,
                               exclude_telescope_type = settings['exclude_telescope_type'], 
                               exclude_channel_type = settings['exclude_channel_type'], 
                               exclude_acquisition_mode = settings['exclude_acquisition_mode'], 
                               exclude_channel_subtype = settings['exclude_channel_subtype'])
        
            # Create the x axis (time)
            x_lbin, x_ubin, x_tick, t_vals, t_tick, = \
                make_axis.quicklook_x(t_lims = settings['t_lims'],
                                      t_tick = settings['t_tick'], 
                                      time = time)
        
            lib = Libraries(
                system_info=system_info,
                channel_info=channel_info,
                time_info=time_info,
                settings=settings,
                )
            
            # iterate over the channels
            for ch in channels:
                
                print(f"-- channel: {ch}")
        
                ch_d = dict(channel = ch)
                
                sig_ch = profiles.loc[ch_d].values
                
                ranges_ch = ranges.copy().loc[ch_d].values
                heights_ch = height_asl.copy().loc[ch_d].values
                    
                # Create the y axis (height/range)
                y_lbin, y_ubin, y_llim, y_ulim, y_vals, y_label = \
                    make_axis.quicklook_y(
                        heights = heights_ch, 
                        ranges = ranges_ch,  
                        y_lims = settings['y_lims'], 
                        use_dis = caller_info['use_range']
                        )
        
                # Smoothing
                if settings['smooth']:
                    if not isinstance(settings['smoothing_window'],list):
                        from visualizer.tools.smoothing import sliding_average_2D_fast as smooth_2D
                    else:
                        from visualizer.tools.smoothing import sliding_average_2D as smooth_2D
        
                    z_vals, _ = smooth_2D(
                        z_vals = sig_ch, 
                        y_vals = y_vals,
                        y_sm_lims = settings['smoothing_range'],
                        y_sm_win = settings['smoothing_window'],
                        expo = settings['smooth_exponential']
                        )
                else:
                    z_vals = sig_ch
        
                # Create the z axis (signal)
                z_llim, z_ulim, z_vals = \
                    make_axis.quicklook_z(z_vals = z_vals, 
                                          y_vals = y_vals,
                                          z_lims = settings['z_lims'] , 
                                          use_log = False,
                                          z_max_zone = settings['z_max_zone'],
                                          z_min_zone = settings['z_min_zone'])

                text_generator = GenerateText(lib = lib, atlas_channel_id = ch)
                
                # Make title
                title = text_generator.make_quicklook_title()
          
                # Make filename
                filename = text_generator.make_filename(qa_test=f'qck_{key}')              

                # Make the plot
                plot_path = generate_plot(
                    dir_out = os.path.join(
                    caller_info['output_folder'],'plots'), 
                    fname = f"{filename}.png",
                    title = title,
                    dpi_val = caller_info['dpi'],
                    use_log = False,
                    delta_t = delta_t,
                    t_vals = t_vals, y_vals = y_vals, 
                    z_vals = z_vals, 
                    x_lbin = x_lbin, x_ubin = x_ubin, 
                    y_lbin = y_lbin, y_ubin = y_ubin, 
                    y_llim = y_llim, y_ulim = y_ulim, 
                    z_llim = z_llim, z_ulim = z_ulim,
                    y_label = y_label, 
                    t_tick = t_tick, x_tick = x_tick,
                    y_tick = settings['y_tick']
                    )  
            
                # Perform color reduction        
                perform_color_reduction(
                    color_reduction = True, 
                    plot_path = plot_path
                    )

             
                # # Add metadata to the quicklook plot
                # plot_metadata = make_plot.get_plot_metadata(metadata = metadata, 
                #                                             args = args, 
                #                                             channel = ch,
                #                                             meas_type = f'qck_{key}', 
                #                                             version = __version__)
                
                # make_plot.add_plot_metadata(plot_path = plot_path, 
                #                             plot_metadata = plot_metadata)

            print('-----------------------------------------')
            print(' ')

    