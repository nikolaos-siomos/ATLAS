#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings, sys
import numpy as np
from .readers.parse_pcb_args import call_parser, check_parser
from .readers.check import check_channels_no_exclude as check_channels
from .readers.check import find_rt_channels
from .readers.read_prepro import unpack
from .plotting import make_axis, make_title, make_plot
from .tools import average, curve_fit
from .writters import make_header, export_ascii 
import os

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
    # Check the command line argument information
    args = check_parser(args)
    
    print('-----------------------------------------')
    print('Initializing the pcb. Calibration...')
    print('-----------------------------------------')
    
    # Read the pcb cal file
    profiles, metadata = unpack(args['input_file'])

    channels_r, channels_t= find_rt_channels(ch_r = args['ch_r'], 
                                             ch_t = args['ch_t'], 
                                             channels = metadata['atlas_channel_id'])
    
    # Check if the parsed channels exist
    channels_r = check_channels(sel_channels = args['ch_r'], 
                                all_channels = channels_r) 
    channels_t = check_channels(sel_channels = args['ch_t'], 
                                all_channels = channels_t)
    
    G_R_def = len(channels_r) * [1.]
    G_T_def = len(channels_r) * [1.]
    H_R_def = len(channels_r) * [1.]
    H_T_def = len(channels_r) * [1.]

    for i in range(len(channels_r)):
        if channels_r[i][5] == 'c' and channels_t[i][5] == 'p':
            H_R_def[i] = -1.
        if channels_r[i][5] == 'p' and channels_t[i][5] == 'c':
            H_T_def[i] = -1.
        if channels_r[i][5] == 't' and channels_t[i][5] == 'p':
            H_R_def[i] =  0.
        if channels_r[i][5] == 'p' and channels_t[i][5] == 't':
            H_T_def[i] =  0.
        if channels_r[i][5] == 't' and channels_t[i][5] == 'c':
            H_R_def[i] =  0.
            H_T_def[i] = -1.
        if channels_r[i][5] == 'c' and channels_t[i][5] == 't':
            H_R_def[i] = -1.
            H_T_def[i] =  0.

    # Extract pair values
    if args['K'] == None:
        K = len(channels_r) * [1.]
        args['K'] = K
    else:
        K = args['K']

    if args['G_R'] == None:
        G_R = G_R_def
        args['G_R'] = G_R
    else:
        G_R = args['G_R']

    if args['G_T'] == None:
        G_T = G_T_def
        args['G_T'] = G_T
    else:
        G_T = args['G_T']

    if args['H_R'] == None:
        H_R = H_R_def
        args['H_R'] = H_R
    else:
        H_R = args['H_R']

    if args['H_T'] == None:
        H_T = H_T_def
        args['H_T'] = H_T
    else:
        H_T = args['H_T']

    if args['R_to_T_transmission_ratio'] == None:
        TR_to_TT = len(channels_r) * [1.]
        args['R_to_T_transmission_ratio'] = TR_to_TT
    else:
        TR_to_TT = args['R_to_T_transmission_ratio']

    # Iterate over the channels
    for i in range(len(channels_r)):
                
        ch_r = channels_r[i]
        ch_t = channels_t[i]
        K_ch = float(K[i])
        G_R_ch = float(G_R[i])
        G_T_ch = float(G_T[i])
        H_R_ch = float(H_R[i])
        H_T_ch = float(H_T[i])
        G_R_def_ch = G_R_def[i]
        G_T_def_ch = G_T_def[i]
        H_R_def_ch = H_R_def[i]
        H_T_def_ch = H_T_def[i]
        TR_to_TT_ch = TR_to_TT[i]
        
        print(f"-- channels: {ch_r} & {ch_t}")

        ch_r_d = dict(channel = ch_r)
        ch_t_d = dict(channel = ch_t)
        
        sig_r_p45_ch = profiles['sig_p45'].copy().loc[ch_r_d].values
        sig_t_p45_ch = profiles['sig_p45'].copy().loc[ch_t_d].values
        sig_r_m45_ch = profiles['sig_m45'].copy().loc[ch_r_d].values
        sig_t_m45_ch = profiles['sig_m45'].copy().loc[ch_t_d].values
        sig_r_ray_ch = profiles['sig_ray'].copy().loc[ch_r_d].values
        sig_t_ray_ch = profiles['sig_ray'].copy().loc[ch_t_d].values
        
        delta_m_prf = profiles['mldr'].loc[ch_r_d].values
        
        ranges_ray_ch = profiles['ranges_ray'].copy().loc[ch_r_d].values
        heights_ray_ch = profiles['heights_ray'].copy().loc[ch_r_d].values
        
        ranges_cal_ch = profiles['ranges_cal'].copy().loc[ch_r_d].values
        heights_cal_ch = profiles['heights_cal'].copy().loc[ch_r_d].values
        

        # Create the y axis (height/range)
        x_lbin_cal, x_ubin_cal, x_llim_cal, x_ulim_cal, x_vals_cal, x_label_cal = \
            make_axis.polarization_calibration_x(
                heights = heights_cal_ch, 
                ranges = ranges_cal_ch,
                x_lims = args['x_lims_calibration'], 
                use_dis = args['use_range'])
    
        # Create the y axis (height/range)
        x_lbin_ray, x_ubin_ray, x_llim_ray, x_ulim_ray, x_vals_ray, x_label_ray = \
            make_axis.polarization_calibration_x(
                heights = heights_ray_ch, 
                ranges = ranges_ray_ch,
                x_lims = args['x_lims_rayleigh'], 
                use_dis = args['use_range'])
    
        # Smoothing
        if args['smooth']== True:
            if not isinstance(args['smoothing_window'],list):
                from .tools.smoothing import sliding_average_1D_fast as smooth_1D
            else:
                from .tools.smoothing import sliding_average_1D as smooth_1D

            y_r_m45_sm, _ = \
                smooth_1D(y_vals = sig_r_m45_ch, 
                          x_vals = x_vals_cal,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'])
    
            y_t_m45_sm, _ = \
                smooth_1D(y_vals = sig_t_m45_ch, 
                          x_vals = x_vals_cal,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'])
                
            y_r_p45_sm, _ = \
                smooth_1D(y_vals = sig_r_p45_ch, 
                          x_vals = x_vals_cal,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'])
    
            y_t_p45_sm, _ = \
                smooth_1D(y_vals = sig_t_p45_ch, 
                          x_vals = x_vals_cal,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'])
    
            y_r_rax_sm, _ = \
                smooth_1D(y_vals = sig_r_ray_ch, 
                          x_vals = x_vals_ray,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'])    
                
            y_t_rax_sm, _ = \
                smooth_1D(y_vals = sig_t_ray_ch, 
                          x_vals = x_vals_ray,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'])      
                
        else:
            y_r_m45_sm = sig_r_m45_ch
            y_t_m45_sm = sig_t_m45_ch
            y_r_p45_sm = sig_r_p45_ch
            y_t_p45_sm = sig_t_p45_ch
            y_r_rax_sm = sig_r_ray_ch
            y_t_rax_sm = sig_t_ray_ch
        
        eta_m45_prf = (y_r_m45_sm / y_t_m45_sm) 
    
        eta_p45_prf = (y_r_p45_sm / y_t_p45_sm)
        
        eta_prf = np.sqrt(eta_m45_prf * eta_p45_prf)
        
        delta_r_prf = (y_r_rax_sm / y_t_rax_sm)
        
        llim = 0.5
        ulim = 11.
        min_win = 1.
        max_win = 4.
        
        if args['calibration_region'][0] < llim:
            llim = args['calibration_region'][0]
        if args['calibration_region'][0] > ulim:
            ulim = args['calibration_region'][1]
        if args['calibration_region'][1] - args['calibration_region'][0] < min_win:
            min_win = args['calibration_region'][1] - args['calibration_region'][0]
        if args['calibration_region'][1] - args['calibration_region'][0] > max_win:
            max_win = args['calibration_region'][1] - args['calibration_region'][0]
        
        # # Check for a fit range for the Δ90 calibration
        # rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef = \
        #     curve_fit.statistics(y1 = eta_prf.copy(),
        #                          y2 = np.ones(eta_prf.shape), 
        #                          x  = x_vals_cal,
        #                          min_win = min_win,
        #                          max_win = max_win,
        #                          step = 0.1,
        #                          llim = llim,
        #                          ulim = ulim,
        #                          rsem_lim = 0.05,
        #                          cross_check_type = 'both',
        #                          cross_check_all_points = False,
        #                          cross_check_crit = 'both',
        #                          der_fac = 2.,
        #                          cancel_shp = True)
        

        # norm_region_cal, idx_cal, fit_cal = \
        #     curve_fit.scan(mfit = mfit,
        #                    dflt_region = args['calibration_region'],
        #                    prefered_range = "near")

        llim = 1.
        ulim = 11.
        min_win = 0.5
        max_win = 4.
        
        if args['rayleigh_region'][0] < llim:
            llim = args['rayleigh_region'][0]
        if args['rayleigh_region'][0] > ulim:
            ulim = args['rayleigh_region'][1]
        if args['rayleigh_region'][1] - args['rayleigh_region'][0] < min_win:
            min_win = args['rayleigh_region'][1] - args['rayleigh_region'][0]
        if args['rayleigh_region'][1] - args['rayleigh_region'][0] > max_win:
            max_win = args['rayleigh_region'][1] - args['rayleigh_region'][0]
        
        # # Check for a fit range for the Rayleigh calibration
        # rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef = \
        #     curve_fit.statistics(y1 = delta_r_prf.copy(),
        #                          y2 = delta_m_prf.copy(), 
        #                          x  = x_vals_cal,
        #                          min_win = min_win,
        #                          max_win = max_win,
        #                          step = 0.1,
        #                          llim = llim,
        #                          ulim = ulim,
        #                          rsem_lim = 0.05,
        #                          cross_check_type = 'both',
        #                          cross_check_all_points = False,
        #                          der_fac = 2.,
        #                          cancel_shp = True)    

        # norm_region_ray, idx_ray, fit_ray = \
        #     curve_fit.scan(mfit = mfit,
        #                    dflt_region = args['rayleigh_region'],
        #                    prefered_range = "far")
                        
        avg_r_m45, _, sem_r_m45 = \
            average.region(sig = sig_r_m45_ch, 
                           x_vals = x_vals_cal, 
                           region = args['calibration_region'], 
                           axis = 0,
                           squeeze = True)
        
        avg_t_m45, _, sem_t_m45 = \
            average.region(sig = sig_t_m45_ch, 
                           x_vals = x_vals_cal, 
                           region = args['calibration_region'], 
                           axis = 0,
                           squeeze = True)
        
        avg_r_p45, _, sem_r_p45 = \
            average.region(sig = sig_r_p45_ch, 
                           x_vals = x_vals_cal, 
                           region = args['calibration_region'], 
                           axis = 0,
                           squeeze = True)
        
        avg_t_p45, _, sem_t_p45 = \
            average.region(sig = sig_t_p45_ch, 
                           x_vals = x_vals_cal, 
                           region = args['calibration_region'], 
                           axis = 0,
                           squeeze = True)
                
        avg_r_ray, _, sem_r_ray = \
            average.region(sig = sig_r_ray_ch, 
                           x_vals = x_vals_ray, 
                           region = args['rayleigh_region'], 
                           axis = 0,
                           squeeze = True)
            
        avg_t_ray, _, sem_t_ray = \
            average.region(sig = sig_t_ray_ch, 
                           x_vals = x_vals_ray, 
                           region = args['rayleigh_region'], 
                           axis = 0,
                           squeeze = True)

        delta_m, _, _ = \
            average.region(sig = delta_m_prf, 
                           x_vals = x_vals_ray, 
                           region = args['rayleigh_region'], 
                           axis = 0,
                           squeeze = True) 
            
        avg_r_m45_i = np.random.normal(loc = avg_r_m45, scale = sem_r_m45, size = 200)
        avg_t_m45_i = np.random.normal(loc = avg_t_m45, scale = sem_t_m45, size = 200)
        avg_r_p45_i = np.random.normal(loc = avg_r_p45, scale = sem_r_p45, size = 200)
        avg_t_p45_i = np.random.normal(loc = avg_t_p45, scale = sem_t_p45, size = 200)
        avg_r_ray_i = np.random.normal(loc = avg_r_ray, scale = sem_r_ray, size = 200)
        avg_t_ray_i = np.random.normal(loc = avg_t_ray, scale = sem_t_ray, size = 200)

        eta_f_s_m45 = (avg_r_m45_i / avg_t_m45_i)
        
        eta_f_s_m45[0] = (avg_r_m45 / avg_t_m45)
        
        eta_f_s_p45 = (avg_r_p45_i / avg_t_p45_i)

        eta_f_s_p45[0] = (avg_r_p45 / avg_t_p45)
        
        eta_f_s = np.sqrt(eta_f_s_p45 * eta_f_s_m45)
        
        eta_s = eta_f_s / TR_to_TT_ch

        eta = eta_s / K_ch

        delta_s_prf = (y_r_rax_sm / y_t_rax_sm) / eta[0]

        delta_s = (avg_r_ray_i / avg_t_ray_i) / eta
            
        delta_s[0] = (avg_r_ray / avg_t_ray) / eta[0]

        delta_c_def_prf = (delta_s_prf * (G_T_def_ch + H_T_def_ch) - (G_R_def_ch + H_R_def_ch)) /\
            ((G_R_def_ch - H_R_def_ch) - delta_s_prf * (G_T_def_ch - H_T_def_ch))
        
        delta_c_def = (delta_s * (G_T_def_ch+ H_T_def_ch) - (G_R_def_ch + H_R_def_ch)) /\
            ((G_R_def_ch - H_R_def_ch) - delta_s * (G_T_def_ch - H_T_def_ch))
            
        delta_c_prf = (delta_s_prf * (G_T_ch + H_T_ch) - (G_R_ch + H_R_ch)) /\
            ((G_R_ch - H_R_ch) - delta_s_prf * (G_T_ch - H_T_ch))

        delta_c = (delta_s * (G_T_ch + H_T_ch) - (G_R_ch + H_R_ch)) /\
            ((G_R_ch - H_R_ch) - delta_s * (G_T_ch - H_T_ch))
                    
        psi = (eta_f_s_p45 - eta_f_s_m45) / (eta_f_s_p45 + eta_f_s_m45)
        
        kappa = 1.
        
        epsilon = np.rad2deg(0.5 * np.arcsin(np.tan(0.5 * np.arcsin(psi) / kappa)))
        # kappa = np.tan(0.5 * np.arcsin(psi)) / np.sin(2. * np.deg2rad(epsilon)) 
        
        delta_l = (delta_c - delta_m) / (1. - delta_c * delta_m)
        
        # delta_l_err = delta_c_err * (1. - delta_m) * (1. + delta_c) / \
        #     (1. - delta_m * delta_c)**2
                      
        # base_delta_v = np.ceil(1E3 * delta_m) / 1E3
        # delta_v = np.hstack((np.arange(base_delta_v, 0.021, 0.001),
        #                      np.arange(0.02, 0.31, 0.01)))
        
        # err_v = delta_l[0]
        err_p = 0.025
        delta_p_err, delta_p, R, sr_lim = pldr_error(delta_m = delta_m, 
                                                     delta_v_err = delta_l[0], 
                                                     delta_p_err_ulim = err_p)
        
        # alpha = (1. + delta_m)**2 * (err_v - err_p)
        # beta = (1. + delta_m) * (2. * err_p * (1. + delta_v + err_v / 2.) - err_v * (1. + delta_m))
        # gamma = - err_p * (1. + delta_v) * (1. + delta_v + err_v)
        
        # sr_lim = (-beta - np.sqrt(beta**2 - 4. * alpha * gamma)) / (2. * alpha)
 
        # Create the y axis (calibration)
        y_llim_cal, y_ulim_cal, y_label_cal = \
            make_axis.polarization_calibration_cal_y(
                ratio_m = eta_f_s_m45[0], ratio_p = eta_f_s_p45[0],
                y_lims_cal = args['y_lims_calibration'])
            
        # Create the y axis (rayleigh)
        y_llim_ray, y_ulim_ray, y_label_ray = \
            make_axis.polarization_calibration_ray_y(
                ratio = delta_c_def[0], y_lims_ray = args['y_lims_rayleigh'])
        
                
        # Make title
        title = make_title.polarization_calibration(channel_r = ch_r, 
                                                    channel_t = ch_t, 
                                                    metadata = metadata, 
                                                    args = args)
        
        # Make plot filename
        fname = make_plot.make_filename(metadata = metadata, 
                                        channel = ch_r,
                                        extra_channel = ch_t,
                                        meas_type = 'pcb', 
                                        version = __version__)
                
        # Make filename
        plot_path = \
            make_plot.polarization_calibration(dir_out = os.path.join(args['output_folder'],'plots'), 
                                               fname = f"{fname}.png", title = title,
                                               dpi_val = args['dpi'],
                                               color_reduction = args['color_reduction'],
                                               cal_region = args['calibration_region'],
                                               vdr_region = args['rayleigh_region'],
                                               x_vals_cal = x_vals_cal, 
                                               x_vals_vdr = x_vals_ray, 
                                               y1_vals = eta_prf, 
                                               y2_vals = eta_p45_prf, 
                                               y3_vals = eta_m45_prf, 
                                               y4_vals = delta_c_def_prf,
                                               y5_vals = delta_c_prf,
                                               y6_vals = delta_m_prf,
                                               eta = eta[0], 
                                               eta_f_s = eta_f_s[0], 
                                               eta_s = eta_s[0], 
                                               delta_m = delta_m,
                                               delta_c_def = delta_c_def[0],
                                               delta_c = delta_c[0],
                                               delta_l = delta_l[0],
                                               epsilon = epsilon[0],
                                               sr_lim = sr_lim,
                                               err_p = err_p,
                                               eta_err = np.std(eta[1:]), 
                                               eta_f_s_err = np.std(eta_f_s[1:]), 
                                               eta_s_err = np.std(eta_s[1:]), 
                                               delta_c_def_err = np.std(delta_c_def[1:]),
                                               delta_c_err = np.std(delta_c[1:]),
                                               delta_l_err = np.std(delta_l[1:]),
                                               epsilon_err = np.std(epsilon[1:]),
                                               x_lbin_cal = x_lbin_cal,
                                               x_ubin_cal = x_ubin_cal, 
                                               x_llim_cal = x_llim_cal,
                                               x_ulim_cal = x_ulim_cal, 
                                               y_llim_cal = y_llim_cal, 
                                               y_ulim_cal = y_ulim_cal, 
                                               x_lbin_vdr = x_lbin_ray, 
                                               x_ubin_vdr = x_ubin_ray, 
                                               x_llim_vdr = x_llim_ray, 
                                               x_ulim_vdr = x_ulim_ray, 
                                               y_llim_vdr = y_llim_ray, 
                                               y_ulim_vdr = y_ulim_ray, 
                                               y_label_cal = y_label_cal, 
                                               x_label_cal = x_label_cal, 
                                               K = K_ch,
                                               G_R = G_R_ch,
                                               G_T = G_T_ch,
                                               H_R = H_R_ch,
                                               H_T = H_T_ch,
                                               x_tick_cal = args['x_tick_calibration'],
                                               y_label_vdr = y_label_ray, 
                                               x_label_vdr = x_label_ray, 
                                               x_tick_vdr = args['x_tick_rayleigh'])  
    
        # Make ascii file header
        header = \
            make_header.polarisation_calibration(channel_r = ch_r,
                                                 channel_t = ch_t,
                                                 metadata = metadata,
                                                 K = K_ch,
                                                 G_R = G_R_ch,
                                                 G_T = G_T_ch,
                                                 H_R = H_R_ch,
                                                 H_T = H_T_ch)
        
        # Export to ascii (Volker's format)        
        export_ascii.polarisation_calibration(dir_out = args['output_folder'], 
                                              fname = f"{fname}.txt", 
                                              alt_cal = x_vals_cal,                                              
                                              alt_ray = x_vals_ray,
                                              r_p45 = sig_r_p45_ch,
                                              t_p45 = sig_t_p45_ch,
                                              r_m45 = sig_r_m45_ch,
                                              t_m45 = sig_t_m45_ch,   
                                              ray_r = sig_r_ray_ch,
                                              ray_t = sig_t_ray_ch,
                                              header = header)
        
        # Add metadata to the plot
        plot_metadata_r = make_plot.get_plot_metadata(metadata = metadata, 
                                                      args = args, 
                                                      channel = ch_r,
                                                      meas_type = 'pcb', 
                                                      version = __version__,
                                                      data_source_id = 'r')
       
        plot_metadata_t = make_plot.get_plot_metadata(metadata = metadata, 
                                                      args = args, 
                                                      channel = ch_t,
                                                      meas_type = 'pcb', 
                                                      version = __version__,
                                                      data_source_id = 't')
        
        make_plot.add_plot_metadata(plot_path = plot_path, 
                                    plot_metadata = plot_metadata_r,
                                    plot_metadata_extra = plot_metadata_t)
        
    print('-----------------------------------------')
    print(' ')
    
    return()

def pldr_error(delta_m, delta_v_err, delta_p_ulim = 0.3, delta_p_err_ulim = 0.025):
    
    R = np.arange(1.01, 3., 0.001)
    delta_p = np.arange(0., delta_p_ulim + 0.001, 0.001)
    
    # delta_p_err = np.zeros((len(delta_p), len(R)))
    
    # for i in range(len(delta_p)):
    sq_term_nom = (delta_v_err + delta_p[:,np.newaxis]) * (1. + delta_m)**2 * np.power(R[np.newaxis,:], 2)
    sq_term_denom = (1. + delta_m)**2 * np.power(R[np.newaxis,:], 2)
    
    lin_term_nom = (1. + delta_m)*(delta_v_err * (delta_p[:,np.newaxis] - 2. * delta_m) - delta_p[:,np.newaxis] * (1. + delta_m)) * R[np.newaxis,:]
    lin_term_denom = -(1. + delta_m)*(delta_v_err + 1. + delta_m) * R[np.newaxis,:]
    
    const_term_nom = -delta_m * (delta_p[:,np.newaxis] - delta_m) * delta_v_err
    const_term_denom = (delta_p[:,np.newaxis] - delta_m) * delta_v_err
        
    delta_p_err = (sq_term_nom + lin_term_nom + const_term_nom) /\
        (sq_term_denom + lin_term_denom + const_term_denom) - delta_p[:,np.newaxis]
        
    delta_p_err[np.abs(delta_p_err) > delta_p_err_ulim] = np.nan
    
    if not np.isnan(delta_p_err[-1,:]).all():
        if delta_v_err > 0.0001:
            min_bsc_ratio = R[np.nanargmax(delta_p_err[-1,:])]
        elif delta_v_err < -0.0001:
            min_bsc_ratio = R[np.nanargmin(delta_p_err[-1,:])]
        else:
            min_bsc_ratio = 1.01
    else:
        min_bsc_ratio = np.nan
        
    return(delta_p_err, delta_p, R, min_bsc_ratio)

if __name__ == '__main__':
    
    sys.path.append('../')
    
    from .version import __version__
    
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args, __version__)
    

# scat_lim_1 = (-beta + np.sqrt(beta**2 - 4. * alpha * gamma)) / (2. * alpha)
# R = np.arange(1.005, 10.005, 0.005)

# def d_p(delta_m, delta_v, R):
    
#     delta_p = np.nan * np.zeros((delta_v.size, R.size))
    
#     for i in range(delta_v.size):
#         crit_1 = (1. + delta_v[i]) * delta_m / ((1. + delta_m) * delta_v[i])
#         crit_2 = (1. + delta_v[i]) / (1. + delta_m)
        
#         for j in range(R.size):
#             if (R[j] >= crit_1 and R[j] >= crit_2) or (R[j] <= -crit_1 and R[j] <= -crit_2):
#                 delta_p[i,j] = \
#                     ((1. + delta_m) * delta_v[i] * R[j] - (1. + delta_v[i]) * delta_m) /\
#                     ((1. + delta_m) * R[j] - (1. + delta_v[i]))
#     return(delta_p)

# delta_p_cor = d_p(delta_m = delta_m, delta_v = delta_v, R = R)
# delta_p_off = d_p(delta_m = delta_m, delta_v = delta_v + err_v, R = R)

# from matplotlib import pyplot as plt
# [X, Y] = np.meshgrid(delta_v, R)
# Z = (delta_p_off - delta_p_cor)
# plt.pcolormesh(X,Y,Z.T[:-1,:-1], vmin =0, vmax=0.05)
# plt.colorbar(label= 'PLDR error')
# plt.title(f'VLDR Offset: {np.round(err_v, decimals = 4)} ')
# # plt.plot(delta_v, scat_lim_1)
# plt.plot(delta_v, sr_lim, c = 'tab:orange')
# plt.xlabel('VLDR')
# plt.ylabel('Scattering Ratio')
# plt.show()
