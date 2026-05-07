#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, os, sys
import numpy as np
from .readers.parse_tlc_args import call_parser, check_parser
from .readers.check import check_channels
from .readers.read_prepro import unpack
from .plotting import make_axis, make_title, make_plot
from .writters import make_header, export_ascii 
from .tools import sector, curve_fit
              
# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print('Initializing the Ring Telecover Test...')
    print('-----------------------------------------')

    profiles, metadata = unpack(args['input_file'])
    
    if 'sig_i' in profiles.keys():
    
        # Check if the parsed channels exist
        channels = \
            check_channels(sel_channels = args['channels'], 
                           all_channels = metadata['atlas_channel_id'],
                           exclude_telescope_type = args['exclude_telescope_type'], 
                           exclude_channel_type = args['exclude_channel_type'], 
                           exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                           exclude_channel_subtype = args['exclude_channel_subtype'])
        
        # iterate over the channels
        for ch in channels:
            print(f"-- channel: {ch}")
            
            ch_d = dict(channel = ch)
            
            ranges_ch = profiles['ranges'].copy().loc[ch_d].values
            heights_ch = profiles['heights'].copy().loc[ch_d].values
            
            sig_i = profiles['sig_i']
            sig_o = profiles['sig_o']
            
            # Create the x axis (height/range)
            x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_tick, x_label = \
                make_axis.telecover_x(heights = heights_ch, 
                                      ranges = ranges_ch,
                                      x_lims = args['x_lims'], 
                                      x_tick = args['x_tick'],
                                      use_dis = args['use_range'],
                                      telescope_type = ch[4])
          
            args['smoothing_range'] = [x_llim, x_ulim]

            iters = np.min([sig_o.time_o.size,
                            sig_i.time_i.size])
        
            extra_sec = {'O' : False,
                         'I' : False}
        
            i_o_ratio = np.nanmean(sig_i.loc[ch_d].values.copy(), axis = 0) /\
                np.nanmean(sig_o.loc[ch_d].values.copy(), axis = 0)
                
            y_m_i = np.nanmean(sig_i.loc[ch_d].values.copy(), axis = 0)
            y_m_o = np.nanmean(sig_o.loc[ch_d].values.copy(), axis = 0)
            sector_dev = np.nanstd(np.vstack((y_m_o,y_m_i)), axis = 0)
                
            llim = 1.
            ulim = 3.
            min_win = 0.5
            max_win = 2.
            
            if args['normalization_region'][0] < llim:
                llim = args['normalization_region'][0]
            if args['normalization_region'][0] > ulim:
                ulim = args['normalization_region'][1]
            if args['normalization_region'][1] - args['normalization_region'][0] < min_win:
                min_win = args['normalization_region'][1] - args['normalization_region'][0]
            if args['normalization_region'][1] - args['normalization_region'][0] > max_win:
                max_win = args['normalization_region'][1] - args['normalization_region'][0]
            
            # # Check for a fit range for the Δ90 calibration
            # stats, masks = \
            #     curve_fit.statistics(y1 = sector_dev,
            #                          y2 = np.ones(i_o_ratio.shape),
            #                          x  = x_vals,
            #                          cancel_stats = ["sph"],
            #                          cross_check_type = 'apply_everywhere',
            #                          cross_check_crit = 'between_both_thresholds',
            #                          cross_check_all_points = False)
                
                    
            # norm_region, idx, fit = \
            #     curve_fit.scan(mfit = mfit,
            #                    dflt_region = args['normalization_region'],
            #                    prefered_range = "near")
                
            coef_o, y_m_o, y_sm_o, y_m_sm_o, y_l_sm_o, y_u_sm_o, \
            coef_extra_o, y_extra_o, y_extra_sm_o, extra_sec['O'] = \
                sector.process(x = x_vals, 
                               y = sig_o.loc[ch_d].values, 
                               iters = iters, 
                               smooth = args['smooth'], 
                               x_sm_lims = [x_llim, x_ulim],
                               x_sm_win = args['smoothing_window'],
                               expo = args['smooth_exponential'],
                               region = args['normalization_region'])
                
            coef_i, y_m_i, y_sm_i, y_m_sm_i, y_l_sm_i, y_u_sm_i, \
            coef_extra_i, y_extra_i, y_extra_sm_i, extra_sec['I'] = \
                sector.process(x = x_vals, 
                               y = sig_i.loc[ch_d].values, 
                               iters = iters, 
                               smooth = args['smooth'], 
                               x_sm_lims = [x_llim, x_ulim],
                               x_sm_win = args['smoothing_window'],
                               expo = args['smooth_exponential'],
                               region = args['normalization_region'])
    
    
            y_llim, y_ulim, y_llim_nr, y_ulim_nr = \
                make_axis.telecover_y(sig = [y_m_sm_o[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_i[slice(x_lbin,x_ubin+1)]],
                                      sig_nr = [coef_o * y_m_sm_o[slice(x_lbin,x_ubin+1)],
                                                coef_i * y_m_sm_i[slice(x_lbin,x_ubin+1)]],
                                      y_lims = args['y_lims'])
            
                    
                    
            # Make title
            title = make_title.telecover(channel = ch, 
                                         metadata = metadata, 
                                         args = args, 
                                         iters = iters,
                                         is_ring = True) 
        
            # Make filename
            # Make filename
            fname = make_plot.make_filename(metadata = metadata, 
                                            channel = ch, 
                                            meas_type = 'tlc_rin', 
                                            version = __version__)
            
            # Make the plot
            plot_path = \
                make_plot.telecover_rin(dir_out = os.path.join(args['output_folder'],'plots'), 
                                        fname = f"{fname}.png", title = title,
                                        dpi_val = args['dpi'],
                                        color_reduction = args['color_reduction'],
                                        norm_region = args['normalization_region'],
                                        use_nonrc = args['use_non_rangecor'],
                                        x_vals = x_vals, 
                                        y1_raw = y_m_o, 
                                        y2_raw = y_m_i, 
                                        y1_vals = y_m_sm_o, 
                                        y2_vals = y_m_sm_i, 
                                        y1_extr = y_extra_sm_o, 
                                        y2_extr = y_extra_sm_i, 
                                        y1_extr_raw = y_extra_o, 
                                        y2_extr_raw = y_extra_i, 
                                        y1_lvar = y_l_sm_o, 
                                        y2_lvar = y_l_sm_i, 
                                        y1_uvar = y_u_sm_o, 
                                        y2_uvar = y_u_sm_i, 
                                        coef_1 = coef_o, 
                                        coef_2 = coef_i, 
                                        coef_extra_1 = coef_extra_o, 
                                        coef_extra_2 = coef_extra_i, 
                                        extra_sec = extra_sec,
                                        ranges = ranges_ch,
                                        x_lbin = x_lbin, x_ubin = x_ubin,
                                        x_llim = x_llim, x_ulim = x_ulim, 
                                        y_llim = y_llim, y_ulim = y_ulim, 
                                        y_llim_nr = y_llim_nr, 
                                        y_ulim_nr = y_ulim_nr, 
                                        x_label = x_label,
                                        x_tick = x_tick,
                                        use_last = args['use_last'],
                                        iters = iters)
            
            sectors = {'O' : y_m_o,
                       'I' : y_m_i}
            
            sectors_e = {'O' : y_extra_o,
                         'I' : y_extra_i}
            
            # Make ascii file header
            header = \
                make_header.telecover(channel = ch, 
                                      metadata = metadata, 
                                      iters = iters, 
                                      extra_sec = extra_sec)
            
            # Export to ascii (Volker's format)        
            export_ascii.telecover(dir_out = args['output_folder'], 
                                   fname = f"{fname}.txt", 
                                   header = header,
                                   iters = 1, 
                                   alt = x_vals, 
                                   sectors = sectors, 
                                   sectors_e = sectors_e)
    
            # Add metadata to the plots
            plot_metadata = make_plot.get_plot_metadata(metadata = metadata, 
                                                        args = args, 
                                                        channel = ch,
                                                        meas_type = 'tlc_rin', 
                                                        version = __version__)
            
            make_plot.add_plot_metadata(plot_path = plot_path, 
                                        plot_metadata = plot_metadata)
    
    print('-----------------------------------------')
    print(' ')
    
    return()

if __name__ == '__main__':
    
    sys.path.append('../')
    
    from .version import __version__
    
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args, __version__)
    
