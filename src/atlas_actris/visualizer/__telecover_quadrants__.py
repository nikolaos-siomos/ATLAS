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
    print('Initializing the Quadr. Telecover Test...')
    print('-----------------------------------------')

    profiles, metadata = unpack(args['input_file'])
        
    if 'sig_n' in profiles.keys():
    
        # Check if the parsed channels exist
        channels = \
            check_channels(sel_channels = args['channels'], 
                           all_channels = metadata['atlas_channel_id'],
                           exclude_telescope_type = args['exclude_telescope_type'], 
                           exclude_channel_type = args['exclude_channel_type'], 
                           exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                           exclude_channel_subtype = args['exclude_channel_subtype'])
        
        # Iterate over the channels
        for ch in channels:
            print(f"-- channel: {ch}")
            
            ch_d = dict(channel = ch)
            
            ranges_ch = profiles['ranges'].copy().loc[ch_d].values
            heights_ch = profiles['heights'].copy().loc[ch_d].values
            
            sig_n = profiles['sig_n']
            sig_e = profiles['sig_e']
            sig_s = profiles['sig_s']
            sig_w = profiles['sig_w']
            
            # Create the x axis (height/range)
            x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_tick, x_label = \
                make_axis.telecover_x(heights = heights_ch, 
                                      ranges = ranges_ch,
                                      x_lims = args['x_lims'], 
                                      x_tick = args['x_tick'],
                                      use_dis = args['use_range'],
                                      telescope_type = ch[4])
            
            args['smoothing_range'] = [x_llim, x_ulim]
        
            iters = np.min([sig_n.time_n.size,
                            sig_e.time_e.size,
                            sig_s.time_s.size,
                            sig_w.time_w.size])
    
            extra_sec = {'N' : False,
                         'E' : False,
                         'S' : False,
                         'W' : False}
            
            n_s_ratio = np.nanmean(sig_n.loc[ch_d].values.copy(), axis = 0) /\
                np.nanmean(sig_s.loc[ch_d].values.copy(), axis = 0)
            
            y_m_n = np.nanmean(sig_n.loc[ch_d].values.copy(), axis = 0)
            y_m_e = np.nanmean(sig_e.loc[ch_d].values.copy(), axis = 0)
            y_m_s = np.nanmean(sig_s.loc[ch_d].values.copy(), axis = 0)
            y_m_w = np.nanmean(sig_w.loc[ch_d].values.copy(), axis = 0)
            sector_dev = np.nanstd(np.vstack((y_m_n,y_m_e,y_m_s,y_m_w)), axis = 0)
    
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
            # rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef = \
            #     curve_fit.statistics(y1 = sector_dev,
            #                          y2 = np.ones(n_s_ratio.shape), 
            #                          x  = x_vals,
            #                          min_win = min_win,
            #                          max_win = max_win,
            #                          step = 0.1,
            #                          llim = llim,
            #                          ulim = ulim,
            #                          rsem_lim = 0.05,
            #                          cancel_shp = True,
            #                          cross_check_type = 'both',
            #                          cross_check_crit = 'both',
            #                          cross_check_all_points = False)
            
                    
            # norm_region, idx, fit = \
            #     curve_fit.scan(mfit = mfit,
            #                    dflt_region = args['normalization_region'],
            #                    prefered_range = "near")
    
            coef_n, y_m_n, y_sm_n, y_m_sm_n, y_l_sm_n, y_u_sm_n, \
            coef_extra_n, y_extra_n, y_extra_sm_n, extra_sec['N'] = \
                sector.process(x = x_vals, 
                               y = sig_n.loc[ch_d].values.copy(), 
                               iters = iters, 
                               smooth = args['smooth'], 
                               x_sm_lims = [x_llim, x_ulim],
                               x_sm_win = args['smoothing_window'],
                               expo = args['smooth_exponential'],
                               region = args['normalization_region'])
                
            coef_e, y_m_e, y_sm_e, y_m_sm_e, y_l_sm_e, y_u_sm_e, \
            coef_extra_e, y_extra_e, y_extra_sm_e, extra_sec['E'] = \
                sector.process(x = x_vals, 
                               y = sig_e.loc[ch_d].values.copy(), 
                               iters = iters, 
                               smooth = args['smooth'], 
                               x_sm_lims = [x_llim, x_ulim],
                               x_sm_win = args['smoothing_window'],
                               expo = args['smooth_exponential'],
                               region = args['normalization_region'])
                
            coef_s, y_m_s, y_sm_s, y_m_sm_s, y_l_sm_s, y_u_sm_s, \
            coef_extra_s, y_extra_s, y_extra_sm_s, extra_sec['S'] = \
                sector.process(x = x_vals, 
                               y = sig_s.loc[ch_d].values.copy(), 
                               iters = iters, 
                               smooth = args['smooth'], 
                               x_sm_lims = [x_llim, x_ulim],
                               x_sm_win = args['smoothing_window'],
                               expo = args['smooth_exponential'],
                               region = args['normalization_region'])
                
            coef_w, y_m_w, y_sm_w, y_m_sm_w, y_l_sm_w, y_u_sm_w, \
            coef_extra_w, y_extra_w, y_extra_sm_w, extra_sec['W'] = \
                sector.process(x = x_vals, 
                               y = sig_w.loc[ch_d].values.copy(), 
                               iters = iters, 
                               smooth = args['smooth'], 
                               x_sm_lims = [x_llim, x_ulim],
                               x_sm_win = args['smoothing_window'],
                               expo = args['smooth_exponential'],
                               region = args['normalization_region'])
    
            # Create the y axis (signal)
            y_llim, y_ulim, y_llim_nr, y_ulim_nr = \
                make_axis.telecover_y(sig = [y_m_sm_n[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_e[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_s[slice(x_lbin,x_ubin+1)],
                                             y_m_sm_w[slice(x_lbin,x_ubin+1)]],
                                      sig_nr = [coef_n * y_m_sm_n[slice(x_lbin,x_ubin+1)],
                                                coef_e * y_m_sm_e[slice(x_lbin,x_ubin+1)],
                                                coef_s * y_m_sm_s[slice(x_lbin,x_ubin+1)],
                                                coef_w * y_m_sm_w[slice(x_lbin,x_ubin+1)]],
                                      y_lims = args['y_lims'])
            
                    
            # Make title
            title = make_title.telecover(channel = ch, 
                                         metadata = metadata, 
                                         args = args, 
                                         iters = iters)        
        
            # Make filename
            fname = make_plot.make_filename(metadata = metadata, 
                                            channel = ch, 
                                            meas_type = 'tlc', 
                                            version = __version__)
            
            # Make the plot
            plot_path = \
                make_plot.telecover_sec(dir_out = os.path.join(args['output_folder'],'plots'), 
                                        fname = f"{fname}.png", 
                                        title = title,
                                        dpi_val = args['dpi'],
                                        color_reduction = args['color_reduction'],
                                        norm_region = args['normalization_region'],
                                        use_nonrc = args['use_non_rangecor'],
                                        x_vals = x_vals, 
                                        y1_raw = y_m_n, 
                                        y2_raw = y_m_e, 
                                        y3_raw = y_m_s, 
                                        y4_raw = y_m_w,   
                                        y1_vals = y_m_sm_n, 
                                        y2_vals = y_m_sm_e, 
                                        y3_vals = y_m_sm_s, 
                                        y4_vals = y_m_sm_w,
                                        y1_extr = y_extra_sm_n, 
                                        y2_extr = y_extra_sm_e, 
                                        y3_extr = y_extra_sm_s, 
                                        y4_extr = y_extra_sm_w,
                                        y1_extr_raw = y_extra_n, 
                                        y2_extr_raw = y_extra_e, 
                                        y3_extr_raw = y_extra_s, 
                                        y4_extr_raw = y_extra_w,
                                        y1_lvar = y_l_sm_n, 
                                        y2_lvar = y_l_sm_e, 
                                        y3_lvar = y_l_sm_s, 
                                        y4_lvar = y_l_sm_w,
                                        y1_uvar = y_u_sm_n, 
                                        y2_uvar = y_u_sm_e, 
                                        y3_uvar = y_u_sm_s, 
                                        y4_uvar = y_u_sm_w,
                                        coef_1 = coef_n, 
                                        coef_2 = coef_e, 
                                        coef_3 = coef_s, 
                                        coef_4 = coef_w,
                                        coef_extra_1 = coef_extra_n, 
                                        coef_extra_2 = coef_extra_e, 
                                        coef_extra_3 = coef_extra_s, 
                                        coef_extra_4 = coef_extra_w,
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
            
            sectors = {'N' : y_m_n,
                       'E' : y_m_e,
                       'S' : y_m_s,
                       'W' : y_m_w}
            
            sectors_e = {'N' : y_extra_n,
                         'E' : y_extra_e,
                         'S' : y_extra_s,
                         'W' : y_extra_w}
            
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
                                                        meas_type = 'tlc', 
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
    
