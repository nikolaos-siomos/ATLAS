#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:24:00 2022

@author: nick
"""

import argparse
import os
import shutil

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

    parser.add_argument('--color_reduction', metavar = 'color_reduction',
                        type = bool, default = True, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If set to True, and provided that Imagemagick is installed, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True. A warning will be raised if Imagemagick is not installed. ')

    parser.add_argument('--use_range', metavar = 'use_range',
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
                        type = float, nargs = 2, default = [0., 15.], 
                        help = 'The x axis limits in km (lower and upper) for the pcb. calibration plot. If use_range is called, the limits correspond to distance. Defaults to 0 km (lower) and 15 km (upper) If values below 0 or above the maximum signal altitude/distance are used, they will be ignored')

    parser.add_argument('--x_tick_calibration', metavar = 'x_tick_calibration',
                        type = int, nargs = '?', default = 2., 
                        help = 'The x axis finest tick in km for the pcb. calibration plot. Defaults to 2 km ')

    parser.add_argument('--x_lims_rayleigh', metavar = 'x_lims_rayleigh',
                        type = float, nargs = 2, default = [0., 30.], 
                        help = 'The x axis limits in km (lower and upper) for the Rayleigh VDR plot. If use_range is called, the limits correspond to distance. Defaults to 0 km (lower) and 15 km (upper) If values below 0 or above the maximum signal altitude/distance are used, they will be ignored')

    parser.add_argument('--x_tick_rayleigh', metavar = 'x_tick_rayleigh',
                        type = int, nargs = '?', default = 2., 
                        help = 'The x axis finest tick in km for the Rayleigh VDR plot. Defaults to 2 km ')

    parser.add_argument('--ch_r', metavar = 'ch_r',
                        type = str, nargs = '+',
                        help = 'Type one or more channel names (e.g. xpar0355) here that correspond to the reflected channels used for the polarization calibration calculation. The number of reflected channels must be the same as the number of the respective transmitted channels ')

    parser.add_argument('--ch_t', metavar = 'ch_t',
                        type = str, nargs = '+', 
                        help = 'Type one or more channel names (e.g. xpat0355) here that correspond to the transmitted channels used for the polarization calibration calculation. The number of transmitted channels must be the same as the number of the respective reflected channels ')

    parser.add_argument('--calibration_region', metavar = 'calibration_region',
                        type = float, nargs = 2, default = [4., 6.], 
                        help = 'The lower and upper limits of the region used for Δ90 calibration. If use_range is called, the limits correspond to distance. Defaults to: 2., 4. km ')

    parser.add_argument('--rayleigh_region', metavar = 'rayleigh_height',
                        type = float, nargs = 2, default = [4, 10.], 
                        help = 'The lower and upper limits of the region used for the comparison with the Rayleigh atmosphere. Defaults to: 8.5, 9.5 ')

    parser.add_argument('--smooth', metavar = 'smooth',
                        type = bool, default = True, 
                        action = argparse.BooleanOptionalAction,
                        help = 'Refer to the smooth option in the quicklook section. Defaults to: True')

    parser.add_argument('--smoothing_range', metavar = 'smoothing_range',
                        type = float, nargs = 2, default = [0., 31.], 
                        help = 'Refer to the smooth option in the quicklook section Defaults to: 0.05, 31.')

    parser.add_argument('--smoothing_window', metavar = 'smoothing_window',
                        type = float, nargs = "?", default = 500., 
                        help = 'The full smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the only one value twice to apply a constant window. Defaults to: smoothing_window = 100.')

    parser.add_argument('--smooth_exponential', metavar = 'smooth',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'Refer to the smooth option in the quicklook section. Defaults to: False.')

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
    
    if args['color_reduction'] == True:
        if shutil.which('convert') == None:
            raise Warning("---- Color_reduction was set tot True for the Rayleigh fit test but Imagemagick was not found. No color reduction will be performed. ")
            args['color_reduction'] = False
            
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
