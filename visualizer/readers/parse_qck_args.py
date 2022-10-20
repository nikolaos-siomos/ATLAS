#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:19:18 2022

@author: nick
"""

import argparse
import os

def call_parser():
        
    """Collects the information included as commandline arguments. 
    """

    print('Parsing Quicklook arguments...')
    
    parser = argparse.ArgumentParser(
    	description='arguments ')
    
    parser.add_argument('-i', '--input_file', metavar = 'input_file', 
                        type = str, nargs = '?', 
                        help = 'The path to the input file ')

    parser.add_argument('-o', '--output_folder', metavar = 'output_folder', 
                        type = str, nargs = '?', 
                        help = 'The path to the output folder where the results and plots subfolders will be placed ')

    parser.add_argument('-d', '--delete', metavar = 'delete',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then the input file will be DELETED after processing in order to save space. Use with care! ')

    parser.add_argument('--dpi', metavar = 'dpi',
                        type = int, nargs = 1, default = 300, 
                        help = 'The dots per inch (dpi) resolution of the exported figures. Defaults to 100 ')

    parser.add_argument('--use_log_scale', metavar = 'use_log_scale',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, a logarithmic scale will be used for the visualization of the signal levels ')

    parser.add_argument('--use_distance', metavar = 'use_distance',
                        type = bool, default = True, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, the y axis of the quicklook will correspond to the distance between the laser pulse and the telescope (vertical range) ')

    parser.add_argument('-c', '--channels', metavar = 'channels',
                        type = str, nargs = '+', default = None, 
                        help = 'Type one or more channel names (e.g. xpar0355) here in order to open the figures in interactive mode ')

    parser.add_argument('--exclude_field_type', metavar = 'exclude_field_type',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel field types that you want to EXCLUDE (None: None, x: unspecified, n: near field, f: far field ). Nothing is excluded by default in ATLAS preprocessor ')

    parser.add_argument('--exclude_detection_mode', metavar = 'exclude_detection_mode',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel detection mode types that you want to EXCLUDE (None: None, a: analogue, p: photon). Nothing is excluded by default in ATLAS preprocessor ')

    parser.add_argument('--exclude_scattering_type', metavar = 'exclude_scattering_type',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel scattering types that you want to EXCLUDE (None: None, p: co-polar linear analyzer, c: cross-polar linear analyzer, t: total (no depol), o: co-polar circular analyzer, x: cross-polar circular analyzer, v: vibrational Raman, r: rotational Raman, a: Cabannes, f: fluorescence). Nothing is excluded by default in ATLAS preprocessor ')

    parser.add_argument('--exclude_channel_subtype', metavar = 'exclude_channel_subtype',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel scattering types that you want to EXCLUDE (None: None, r: Signal Reflected from a PBS, t: Signal Transmitted through a PBS, n: N2 Ramal line, o: O2 Ramal line, w: H2O Ramal line, c: CH4 Ramal line, h: High Rotational Raman, l: Low Rotational Raman, a: Mie (aerosol) HSRL signal, m: Molecular HSRL signal, b: Broadband Fluorescence, s: Spectral Fluorescence, x: No specific subtype). Nothing is excluded by default in ATLAS preprocessor ')
    
    parser.add_argument('--x_lims', metavar = 'x_lims',
                        type = int, nargs = 2, default = [None, None], 
                        help = 'The x axis limits (lower and upper). Use two integers corresponding to the first and last timeframe (not date!) that will be plotted. Use 1 to start from the first timeframe. If values below 1 or above the total number of timeframes are used, they will be ignored')

    parser.add_argument('--x_tick', metavar = 'x_tick',
                        type = int, nargs = 1, default = None, 
                        help = 'The x axis finest tick in number of timeframes. ')

    parser.add_argument('--t_lims', metavar = 't_lims',
                        type = str, nargs = 2, default = [None, None], 
                        help = 'The x axis limits in time units. Use the following format hh:mm eg. --t_lims 17:38 21:56')
    
    parser.add_argument('--t_tick', metavar = 't_tick',
                        type = int, nargs = 1, default = None, 
                        help = 'The t axis finest tick in minutes. ')

    parser.add_argument('--y_lims', metavar = 'y_lims',
                        type = float, nargs = 2, default = [0., 14.], 
                        help = 'The y axis limits in km (lower and upper). If use_distance is called, the limits correspond to distance. Defaults to the altitude/distance of the lowermost (lower) and the highermost (upper) bin. If values below 0 or above the maximum altitude/distance are used, they will be ignored')

    parser.add_argument('--y_tick', metavar = 'y_tick',
                        type = int, nargs = 1, default = 1, 
                        help = 'The y axis finest tick in km. Defaults to 1km ')
    
    parser.add_argument('--z_lims', metavar = 'z_lims',
                        type = float, nargs = 2, default = [None, None], 
                        help = 'The colorscale limits of the normalized RC signal. Defaults to 0 (lower) 1 (upper) when use_log_scale is False. If use_log_scale is true then the lower limit becomes 1E-5 ')

    parser.add_argument('--z_max_zone', metavar = 'z_max_zone',
                        type = float, nargs = 2, default = [0.1, 2], 
                        help = 'Provide the zone (min and max height/distance) in km that will be used for the calculation of the max signal value. The signals are normalized to 1 with this value, that is the uppermost limit of the colorscale. Particularly useful in order to avoid scaling the colors with a cloud ')

    parser.add_argument('--smooth', metavar = 'smooth',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, a sliding average smoothing will be applied on the signals for better visualization ')

    parser.add_argument('--smooth_exponential', metavar = 'smooth',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, an exponentially increasing wind (base 2) will be applied. Defaut to a linearly increasing window ')

    parser.add_argument('--smoothing_range', metavar = 'smoothing_range',
                        type = float, nargs = 2, default = [1., 14.], 
                        help = 'First and last altitude/distance boundaries (depending on the selection of use_distance) where smoothing will be applied, in km. Defines the signal region to be smoothed. If they exceed the current signal boundaries the existing boundaries will be used instead ')

    parser.add_argument('--half_window', metavar = 'half_window',
                        type = float, nargs = 2, default = [100., 100.], 
                        help = 'Half smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the same value twice to apply a constant window ')

    args = vars(parser.parse_args())

    return(args)

def check_parser(args):
    

    mandatory_args = ['input_file']

    mandatory_args_abr = ['-i']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            raise Exception(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')            

    if not os.path.exists(args['input_file']):
        raise Exception(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}!")  
    
    if '_qck_' not in os.path.basename(args['input_file']):
        raise Exception('---- Error: Measurement filename not understood! The filename should contain the _qck_ field (quicklook)')
    
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..','..', 'plots')
        args['output_folder'] = out_path
        os.makedirs(args['output_folder'], exist_ok = True)
    elif not os.path.exists(args['output_folder'] ):
        raise Exception(f"The provided output folder {args['output_folder']} does not exist! Please use an existing folder or don't provide one and let the the parser create the default output folder ") 
        
    return(args)

def view_parser(args):
    
    print(" ")
    print("-- Quicklook arguments!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")
    
    return()