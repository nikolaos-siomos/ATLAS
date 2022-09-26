#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:25:41 2022

@author: nick
"""

import argparse
import os

def call_parser():
        
    """Collects the information included as commandline arguments. 
    """

    print('Parsing Intercomparison arguments...')
    
    parser = argparse.ArgumentParser(
    	description='arguments ')
    
    parser.add_argument('-i', '--input_files', metavar = 'input_files', 
                        type = str, nargs = 2, 
                        help = 'The path to the two input files (place the reference system second) ')

    parser.add_argument('-o', '--output_folder', metavar = 'output_folder', 
                        type = str, nargs = '?', 
                        help = 'The path to the output folder where the results and plots subfolders will be placed ')

    parser.add_argument('-d', '--delete', metavar = 'delete',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then the input file will be DELETED after processing in order to save space. Use with care! ')

    parser.add_argument('--dpi', metavar = 'dpi',
                        type = int, nargs = '?', default = 100, 
                        help = 'The dots per inch (dpi) resolution of the exported figures. Defaults to 100 ')

    parser.add_argument('--use_lin_scale', metavar = 'use_lin_scale',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, a linear scale will be used for the z axis (signal) ')

    parser.add_argument('--use_distance', metavar = 'use_distance',
                        type = bool, default = True, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, the y axis of the quicklook will correspond to the distance between the laser pulse and the telescope (vertical range) ')

    parser.add_argument('--y_lims', metavar = 'y_lims',
                        type = float, nargs = 2, default = [None, None], 
                        help = 'The y axis limits (lower and upper) of the normalized RC signal. Defaults to 0 (lower) 1.2 (upper) when use_lin_scale is True. If use_lin_scale is true then the lower limit becomes 1E-5 ')

    parser.add_argument('--x_lims', metavar = 'x_lims',
                        type = float, nargs = 2, default = [0., 14.], 
                        help = 'The x axis limits in km (lower and upper). If use_distance is called, the limits correspond to distance. Defaults to 0 km (lower) and 14 km (upper) If values below 0 or above the maximum signal altitude/distance are used, they will be ignored')

    parser.add_argument('--x_tick', metavar = 'x_tick',
                        type = int, nargs = '?', default = 1, 
                        help = 'The x axis finest tick in km. Defaults to 1km ')

    parser.add_argument('--channels_1', metavar = 'channels_1',
                        type = str, nargs = '+', default = None, 
                        help = 'Type one or more channel names from the ones included in the first file (e.g. xpar0355) here in order to define a channel to be compared. Each channel here corresponds to one channel provided in channel_2 ')

    parser.add_argument('--channels_2', metavar = 'channels_2',
                        type = str, nargs = '+', default = None, 
                        help = 'Type one or more channel names from the ones included in the second file (e.g. xpar0355) here in order to define a channel to be compared. Each channel here corresponds to one channel provided in channel_1 ')

    parser.add_argument('--normalization_height', metavar = 'normalization_height',
                        type = float, nargs = '?', default = 9., 
                        help = 'The reference height/distance where the signals will be normalized for the Rayleigh fit. If use_distance is called, the limits correspond to distance. Defaults to 9 km ')

    parser.add_argument('--half_normalization_window', metavar = 'half_normalization_window',
                        type = float, nargs = '?', default = 500., 
                        help = 'The half window in meters used for the normalization. Defaults to 100 m ')

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
    
    mandatory_args = ['input_files','output_folder','channels_1', 'channels_2']

    mandatory_args_abr = ['-i','-o','--channels_1','--channels_2']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            raise Exception(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')            
        
    print(" ")
    print("-- Intercomparison arguments!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_files'][0]):
        raise Exception(f"-- Error: The path to the first input file does not exists. Please provide a valid input file path. Current Path: {args['input_file'][0]}")  

    if not os.path.exists(args['input_files'][1]):
        raise Exception(f"-- Error: The path to the second input file does not exists. Please provide a valid input file path. Current Path: {args['input_file'][1]}")  
        
    if os.path.basename(args['input_files'][0])[:3] != 'ray' or \
        os.path.basename(args['input_files'][1])[:3] != 'ray':
        raise Exception('---- Error: Measurement filename not understood! Please start the filename with ray (rayleigh fit)')
    
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..','atlas_visualizer', 'cmp')
        args['output_folder'] = out_path
        os.makedirs(args['output_folder'], exist_ok = True)
    elif not os.path.exists(args['output_folder'] ):
        raise Exception(f"The provided output folder {args['output_folder']} does not exist! Please use an existing folder or don't provide one and let the the parser create the default output folder ") 

    if len(args['channels_1']) != len(args['channels_2']):
        raise Exception("-- Error: The number of channels provided in channels_1 is diffirent than the number of channels provided in channels_2. Please provide pair of channels that are going to be intercompared")          
    
    if len(args['channels_1']) != len(args['channels_2']):
        raise Exception("-- Error: The number of channels provided in channels_1 is diffirent than the number of channels provided in channels_2. Please provide pair of channels that are going to be intercompared")          

    return(args)