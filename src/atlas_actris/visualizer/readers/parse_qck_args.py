#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:19:18 2022

@author: nick
"""

import argparse
import os
import shutil

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

    parser.add_argument('--color_reduction', metavar = 'color_reduction',
                        type = bool, default = True, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If set to True, and provided that Imagemagick is installed, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True. A warning will be raised if Imagemagick is not installed. ')

    parser.add_argument('--use_log_scale', metavar = 'use_log_scale',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If set to True, a logarithmic scale will be used for the visualization of the z_axis (signal levels). It affects the selection of the default z_lims. Defaults to False.')

    parser.add_argument('--use_range', metavar = 'use_range',
                        type = bool, default = True, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, the y axis of the quicklook will correspond to the distance between the laser pulse and the telescope (vertical range) ')

    parser.add_argument('-c', '--channels', metavar = 'channels',
                        type = str, nargs = '+', default = None, 
                        help = 'Type one or more channel names (e.g. xpar0355) here in order to open the figures in interactive mode ')

    parser.add_argument('--exclude_telescope_type', metavar = 'exclude_telescope_type',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel field types that you want to EXCLUDE (None: None, x: unspecified, n: near field, f: far field ). Nothing is excluded by default in ATLAS preprocessor ')

    parser.add_argument('--exclude_acquisition_mode', metavar = 'exclude_acquisition_mode',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel detection mode types that you want to EXCLUDE (None: None, a: analogue, p: photon). Nothing is excluded by default in ATLAS preprocessor ')

    parser.add_argument('--exclude_channel_type', metavar = 'exclude_channel_type',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel scattering types that you want to EXCLUDE (None: None, p: co-polar linear analyzer, c: cross-polar linear analyzer, t: total (no depol), o: co-polar circular analyzer, x: cross-polar circular analyzer, v: vibrational Raman, r: rotational Raman, a: Cabannes, f: fluorescence). Nothing is excluded by default in ATLAS preprocessor ')

    parser.add_argument('--exclude_channel_subtype', metavar = 'exclude_channel_subtype',
                        type = str, nargs = '+', default = [], 
                        help = 'Provide all the channel scattering types that you want to EXCLUDE (None: None, r: Signal Reflected from a PBS, t: Signal Transmitted through a PBS, n: N2 Ramal line, o: O2 Ramal line, w: H2O Ramal line, c: CH4 Ramal line, h: High Rotational Raman, l: Low Rotational Raman, a: Mie (aerosol) HSRL signal, m: Molecular HSRL signal, b: Broadband Fluorescence, s: Spectral Fluorescence, x: No specific subtype). Nothing is excluded by default in ATLAS preprocessor ')
    
    parser.add_argument('--t_lims', metavar = 't_lims',
                        type = int, nargs = 2, default = [None, None], 
                        help = 'The x axis limits in time units. Use the following format: hhmm for both limits (not meant to be used for daylong quicklooks). Defaults to: automatic selection.')
    
    parser.add_argument('--t_tick', metavar = 't_tick',
                        type = int, nargs = 1, default = None, 
                        help = 'The t axis finest tick in minutes. ')

    parser.add_argument('--y_lims', metavar = 'y_lims',
                        type = float, nargs = 2, default = [0., 14.], 
                        help = 'The y axis limits in km (lower and upper). If use_range is called, the limits correspond to distance. Defaults to the altitude/distance of the lowermost (lower) and the highermost (upper) bin. If values below 0 or above the maximum altitude/distance are used, they will be ignored')

    parser.add_argument('--y_tick', metavar = 'y_tick',
                        type = int, nargs = 1, default = 1, 
                        help = 'The y axis finest tick in km. Defaults to 1km ')
    
    parser.add_argument('--z_lims', metavar = 'z_lims',
                        type = float, nargs = 2, default = [None, 1.], 
                        help = 'The colorscale limits of the normalized RC signal. Defaults to 0 (lower) 1 (upper) when use_log_scale is False. If use_log_scale is true then the lower limit becomes 1E-5 ')

    parser.add_argument('--z_max_zone', metavar = 'z_max_zone',
                        type = float, nargs = 2, default = [0.1, 2], 
                        help = 'Provide the zone (min and max height/distance) in km that will be used for the calculation of the max signal value. The signals are normalized to 1 with this value, that is the uppermost limit of the colorscale. Particularly useful in order to avoid scaling the colors with a cloud ')

    parser.add_argument('--z_min_zone', metavar = 'z_min_zone',
                        type = float, nargs = 2, default = [2., 10.], 
                        help = 'Provide the zone (min and max height/distance) in km that will be used for the calculation of the min signal value. This is used in case the z_lims are not explicitly provided when the colorscale is logarithmic')

    parser.add_argument('--smooth', metavar = 'smooth',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If set to True, a sliding average smoothing filter will be applied on the signals across y axis for better visualization. The smoothing_exponential, smoothing_exponential, and smoothing_window will be ignored if smooth is set to False. Defaults to: False')

    parser.add_argument('--smoothing_range', metavar = 'smoothing_range',
                        type = float, nargs = 2, default = [0., 14.], 
                        help = 'Set the first and last altitude/distance boundaries in km where smoothing should be applied. If they exceed the actual signal boundaries the actual boundaries will be used instead. Defaults to: 0.05, 14.')

    parser.add_argument('--smoothing_window', metavar = 'smoothing_window',
                        type = float, nargs = "?", default = 100., 
                        help = 'The full smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the only one value twice to apply a constant window. Defaults to: smoothing_window = 100.')

    parser.add_argument('--smooth_exponential', metavar = 'smooth',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'This variable is ignored if the upper and lower values of smoothing_window are the same. Choose one of: True: a smoothing window that exponentially increases with altitude/distance will be applied, False:  a smoothing window that exponentially increases with altitude/distance will be applied. Defaults to: True.')

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

    if args['color_reduction'] == True:
        if shutil.which('convert') == None:
            raise Warning("---- Color_reduction was set tot True for the Rayleigh fit test but Imagemagick was not found. No color reduction will be performed. ")
            args['color_reduction'] = False
        
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..','..')
        args['output_folder'] = out_path
        os.makedirs(os.path.join(args['output_folder'], 'plots'), exist_ok = True)
    elif not os.path.exists(args['output_folder'] ):
        raise Exception(f"The provided output folder {args['output_folder']} does not exist! Please use an existing folder or don't provide one and let the the parser create the default output folder ") 
    
    if args['t_lims'][0] != None:
        if len(str(args['t_lims'][0])) != 4: raise Exception("Wrong formating of the lower of the t_lims at the quicklook section of the settings file")

    if args['t_lims'][-1] != None:
        if len(str(args['t_lims'][-1])) != 4: raise Exception("Wrong formating of the upper of the t_lims at the quicklook section of the settings file")

    if args['t_lims'][0] != None and args['t_lims'][-1] != None:
        if len(str(args['t_lims'][0])) > len(str(args['t_lims'][-1])): raise Exception("The upper value of t_lims is lower than the lower value at the quicklook section of the settings file")
                                            
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