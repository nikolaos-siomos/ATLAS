#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:24:00 2022

@author: nick
"""

import argparse
import os

def call_parser():
        
    """Collects the information included as commandline arguments. 
    """

    print('Parsing Pol. Calibration arguments...')
    
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
                        type = int, nargs = '?', default = 300, 
                        help = 'The dots per inch (dpi) resolution of the exported figures. Defaults to 100 ')

    parser.add_argument('--use_distance', metavar = 'use_distance',
                        type = bool, default = True, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, the x axis of the quicklook will correspond to the distance between the laser pulse and the telescope (vertical range) ')

    parser.add_argument('--y_lims_calibration', metavar = 'y_lims_clibration',
                        type = float, nargs = 2, default = [None, None], 
                        help = 'The y axis limits (lower and upper) of the gain ratios at +-45. ')

    parser.add_argument('--y_lims_rayleigh', metavar = 'y_lims_rayleigh',
                        type = float, nargs = 2, default = [None, None], 
                        help = 'The y axis limits (lower and upper) of the gain ratios at +-45. ')

    parser.add_argument('--x_lims_calibration', metavar = 'x_lims_calibration',
                        type = float, nargs = 2, default = [0., 5.], 
                        help = 'The x axis limits in km (lower and upper) for the pcb. calibration plot. If use_distance is called, the limits correspond to distance. Defaults to 0 km (lower) and 14 km (upper) If values below 0 or above the maximum signal altitude/distance are used, they will be ignored')

    parser.add_argument('--x_tick_calibration', metavar = 'x_tick_calibration',
                        type = int, nargs = '?', default = 0.5, 
                        help = 'The x axis finest tick in km for the pcb. calibration plot. Defaults to 0.5km ')

    parser.add_argument('--x_lims_rayleigh', metavar = 'x_lims_rayleigh',
                        type = float, nargs = 2, default = [0., 10.], 
                        help = 'The x axis limits in km (lower and upper) for the Rayleigh VDR plot. If use_distance is called, the limits correspond to distance. Defaults to 0 km (lower) and 14 km (upper) If values below 0 or above the maximum signal altitude/distance are used, they will be ignored')

    parser.add_argument('--x_tick_rayleigh', metavar = 'x_tick_rayleigh',
                        type = int, nargs = '?', default = 1, 
                        help = 'The x axis finest tick in km for the Rayleigh VDR plot. Defaults to 1km ')

    parser.add_argument('--ch_r', metavar = 'ch_r',
                        type = str, nargs = '+',
                        help = 'Type one or more channel names (e.g. xpar0355) here that correspond to the reflected channels used for the polarization calibration calculation. The number of reflected channels must be the same as the number of the respective transmitted channels ')

    parser.add_argument('--ch_t', metavar = 'ch_t',
                        type = str, nargs = '+', 
                        help = 'Type one or more channel names (e.g. xpat0355) here that correspond to the transmitted channels used for the polarization calibration calculation. The number of transmitted channels must be the same as the number of the respective reflected channels ')

    parser.add_argument('--calibration_height', metavar = 'calibration_height',
                        type = float, nargs = '?', default = 3., 
                        help = 'The calibration height/distance where the signals will be normalized for the Rayleigh fit. If use_distance is called, the limits correspond to distance. Defaults to 5 km ')

    parser.add_argument('--half_calibration_window', metavar = 'half_calibration_window',
                        type = float, nargs = '?', default = 500., 
                        help = 'The half window in meters used for the calibration. Defaults to 500 m ')

    parser.add_argument('--rayleigh_height', metavar = 'rayleigh_height',
                        type = float, nargs = '?', default = 9., 
                        help = 'The height/distance where the atmosphere contains only insignificant amounts of aerosols. If use_distance is called, the limits correspond to distance. Defaults to 9 km ')

    parser.add_argument('--half_rayleigh_window', metavar = 'half_rayleigh_window',
                        type = float, nargs = '?', default = 500., 
                        help = 'The half window in meters used for the molecular caluclations. Defaults to 100 m ')

    parser.add_argument('--smooth', metavar = 'smooth',
                        type = bool, default = True, 
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

    parser.add_argument("--K", metavar = 'K',
                        type = float, nargs = '+', default = None, 
                        help = 'The K value for each channel pair. Defaults to 1 for all channels ')

    parser.add_argument("--G_R", metavar = 'G_R',
                        type = float, nargs = '+', default = None, 
                        help = 'The G value for the reflected channel of the pair. Defaults to 1 for all channels (no receiver optics + emitted pcb. state correction) ')

    parser.add_argument("--G_T", metavar = 'G_T',
                        type = float, nargs = '+', default = None, 
                        help = 'The G value for the transmitted channel of the pair. Defaults to 1 for all channels (no receiver optics + emitted pcb. state correction) ')

    parser.add_argument("--H_R", metavar = 'H_R',
                        type = float, nargs = '+', default = None, 
                        help = 'The H value for the reflected channel of the pair. Defaults to 1 or -1 for all co-polar (p) and cross-polar (c) reflected channels, respectively (no receiver optics + emitted pcb. state correction) ')

    parser.add_argument("--H_T", metavar = 'H_T',
                        type = float, nargs = '+', default = None, 
                        help = 'The H value for the transmitted channel of the pair. Defaults to 1 or -1 for all co-polar (p) and cross-polar (c) transmitted channels, respectively (no receiver optics + emitted pcb. state correction) ')

    parser.add_argument("--R_to_T_transmission_ratio", metavar = 'R_to_T_transmission_ratio',
                        type = float, nargs = '+', default = None, 
                        help = 'The transmission ratio between the R to T channels per pair (T_R/T_T). The transmission in path R or T is 1. if no filter was applied during calibration. Defaults to 1 for all pairs')

    args = vars(parser.parse_args())

    return(args)

def check_parser(args):
    
    mandatory_args = ['input_file']

    mandatory_args_abr = ['-i']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            raise Exception(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')            

    if not os.path.exists(args['input_file']):
        raise Exception(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}")  
    
    if '_pcb_' not in os.path.basename(args['input_file']):
        raise Exception('---- Error: Measurement filename not understood! The filename should contain the _pcb_ field (polarization_calibration)')
    
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..','..')
        args['output_folder'] = out_path
        os.makedirs(os.path.join(args['output_folder'], 'plots'), exist_ok = True)
        os.makedirs(os.path.join(args['output_folder'], 'ascii'), exist_ok = True)
    elif not os.path.exists(args['output_folder'] ):
        raise Exception(f"The provided output folder {args['output_folder']} does not exist! Please use an existing folder or don't provide one and let the the parser create the default output folder ") 

    if args['ch_r'] != None and args['ch_t'] != None:
        if len(args['ch_r']) != len(args['ch_t']):
            raise Exception('---- Error: The number of reflected channels is different from the number of transmitted channels! Please provide pairs of trasmitted and reflected channels')
        
    return(args)

def view_parser(args):
    
    print(" ")
    print("-- Pol. Calibration arguments!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")
    
    return()