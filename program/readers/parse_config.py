"""
@authors: Nikolaos Siomos (nikolaos.siomos@lmu.de)

================
General:
    Parses the command line arguments provided when calling __main_ 
    The information is stored in a python dictionary
    
Returns:
    
    args:
        A dictionary with all the information provided as command line arguments
        
"""

import argparse
import os
import sys
import numpy as np

def main_parser():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""

    print('-----------------------------------------')
    print('Start parsing command line options...')
    print('-----------------------------------------')
    
    parser = argparse.ArgumentParser(
    	description='arguments ')
    
    parser.add_argument('-i', '--input_file', metavar = 'input_file', 
                        type = str, nargs = '?', 
                        help = 'The path to the input file ')

    parser.add_argument('-o', '--output_folder', metavar = 'output_folder', 
                        type = str, nargs = '?', 
                        help = 'The path to the output folder where the results and plots subfolders will be placed. Defaults to the ./out/atlas folder inside the parent folder')

    parser.add_argument('-d', '--debug', metavar = 'debug',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then debugging files will be generated in ./results/debug folder. These included the metadata gathered from the configuration file, the licel header, and the combination of the two. Default to False ')

    parser.add_argument('-q', '--quicklook', metavar = 'quicklook', 
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then quiclook preprocessed files will be generated (original temporal resolution). ')
                                  
    # parser.add_argument('-e', '--error_simulation', metavar = 'error_simulation', 
    #                     type = bool,  default = False, 
    #                     action = argparse.BooleanOptionalAction,
    #                     help = 'If called then a Monte Carlo Error simulation will be perform (warning! the processing time will increase proportionally to the number of iterations (100 by default). ')

    parser.add_argument('--signal_smoothing', metavar = 'signal_smoothing', 
                        type = bool,  default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then a simple sliding average smoothing will be performed. Not relevant to QA file preprocessing')

    parser.add_argument('--skip_dead_time_correction', metavar = 'skip_dead_time_correction', 
                        type = bool,  default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then the dead time correction will not be performed. ')

    parser.add_argument('--skip_dark_subtraction', metavar = 'skip_dark_subtraction', 
                        type = bool,  default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then the dark background will not be subtracted. ')

    parser.add_argument('--vertical_trimming', metavar = 'vertical_trimming', 
                        type = bool,  default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then bins above a certain altitude (20km by default) will be removed. ')
    
    # parser.add_argument('--iterations', metavar = 'iterations', 
    #                     type = int, nargs = '?', default = 100,
    #                     help = 'The number of iterations for the Monte Carlo error simulation. Applicable only if -e is set to True. Default to 100 ')

    parser.add_argument('--timescale', metavar = 'timescale', 
                        type = int, nargs = '?', default = -1,
                        help = "The temporal window applied for the signal averaging in seconds. Not relevant for QA file preprocessing. If provided then the profiles will be averaged in as many timeframes fit among the measurement start time and end time. Select -1 (default) to average all the profiles. ")

    parser.add_argument('--vertical_limit', metavar = 'vertical_limit', 
                        type = float, nargs = '?', default = 20000.,
                        help = "The maximum altitude in m above which no calculations will be performed. Solar background calculations are performed prior to vertical signal trimming to enable background calculations up to the maximum signal altitude. Defaults to 20km ")                                                                                                        

    parser.add_argument('--smoothing_window', metavar = 'smoothing_window', 
                        type = int, nargs = '?', default = 25,
                        help = 'Smoothing window in bins, common for all channels. Defaults to 25 bins ')

    parser.add_argument('--smoothing_sbin', metavar = 'smoothing_sbin', 
                        type = int, nargs = '?', default = -1,
                        help = 'The starting bin for smoothing. No smoothing will be applied before it. Defaults to -1 that corresponds to the first signal bin ')

    parser.add_argument('--smoothing_ebin', metavar = 'smoothing_ebin', 
                        type = int, nargs = '?', default = -1,
                        help = 'The ending bin for smoothing. No smoothing will be applied after it. Defaults to -1 that corresponds to the last signal bin ')

    args = vars(parser.parse_args())

    mandatory_args = ['input_file']

    mandatory_args_abr = ['-f']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            print(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')
            sys.exit('-- Program stopped') 
            
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..','atlas')
        args['output_folder'] = out_path
        os.makedirs(out_path, exist_ok = True)
        
    print("-- The following values have been used!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_file']):
        sys.exit(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}!")  
    
    if os.path.basename(args['input_file'])[:3] not in ['ray','tlc','pcl','drk']:
        sys.exit('---- Error: Measurement filename not understood! Please start the filename with ray (rayleigh), tlc (telecover), or pcl(polarization calibration) ')
    
    if not os.path.exists(args['output_folder']):
        sys.exit(f"-- Error: The path to the output folder does not exists. Please provide a output folder file path. Current Path: {args['output_folder']}!")  
    

    return(args)

    # parser.add_argument('--smoothing_range_nodes', metavar = 'smoothing_range_nodes', 
    #                     type = float, nargs = '+', default = None,
    #                     help = 'Define the range nodes where the smoothing window changes in m. For example --smoothing_range_nodes 0 1000 4000 6000 Defines 3 different smoothing regions, one between 0 and 1000m, a second betweeen 1000 and 4000m, and a third between 4000 and 60000m ')

    # parser.add_argument('--smoothing_window_nodes', metavar = 'smoothing_window_nodes', 
    #                     type = float, nargs = '+', default = None,
    #                     help = 'Smoothing type. Choose among (0: No smoothing, 1: Linear sliding average with single window, 2: Linear sliding average with range dependent windows, 3: Removal of Outliers. Default to 0 ')

    # parser.add_argument('--smoothing_polynomial_order', metavar = 'smoothing_polynomial_order', 
    #                     type = int, nargs = '+', default = 1,
    #                     help = 'The order of polynomial smoothing. Will be applied only if the signal_smoothing is set to 2. Default to 1 (least squares fit) ')

    # parser.add_argument('--reference_height', metavar = 'reference_height', 
    #                     type = float, nargs = '+', default = np.nan,
    #                     help = "The height of a molecular region in m, used for the rayleigh fit test and scattering ratio calculations. If not provided an automated identification will be atempted. ")

def quicklook_parser():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""

    print('-----------------------------------------')
    print('Start parsing command line options...')
    print('-----------------------------------------')
    
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
                        type = int, nargs = 1, default = 100, 
                        help = 'The dots per inch (dpi) resolution of the exported figures. Defaults to 100 ')

    parser.add_argument('--use_log_scale', metavar = 'use_log_scale',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, a logarithmic scale will be used for the visualization of the signal levels ')

    parser.add_argument('--use_distance', metavar = 'use_distance',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, the y axis of the quicklook will correspond to the distance between the laser pulse and the telescope (vertical range) ')

    parser.add_argument('--x_lims', metavar = 'x_lims',
                        type = int, nargs = 2, default = [None, None], 
                        help = 'The x axis limits (lower and upper). Use two integers corresponding to the first and last timeframe (not date!) that will be plotted. Use 1 to start from the first timeframe. If values below 1 or above the total number of timeframes are used, they will be ignored')

    parser.add_argument('--x_tick', metavar = 'x_tick',
                        type = int, nargs = 1, default = None, 
                        help = 'The x axis finest tick in number of timeframes. ')

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
                        type = float, nargs = 2, default = [None, None], 
                        help = 'Provide the zone (min and max height/distance) in km that will be used for the calculation of the max signal value. The signals are normalized to 1 with this value, that is the uppermost limit of the colorscale. Particularly useful in order to avoid scaling the colors with a cloud ')

    parser.add_argument('-c', '--channels', metavar = 'channels',
                        type = str, nargs = '+', default = None, 
                        help = 'Type one or more channel names (e.g. xpar0355) here in order to open the figures in interactive mode ')

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

    mandatory_args = ['input_file']

    mandatory_args_abr = ['-i']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            print(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')
            sys.exit('-- Program stopped')    

    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..',
                                'plots','qck')
        args['output_folder'] = out_path
        os.makedirs(out_path, exist_ok = True)
        
    print("-- The following values have been used!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_file']):
        sys.exit(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}!")  
    
    if os.path.basename(args['input_file'])[:3] != 'qck':
        sys.exit('---- Error: Measurement filename not understood! Please start the filename with qck (quicklook)')
    
    if not os.path.exists(args['output_folder']):
        sys.exit(f"-- Error: The path to the output folder does not exists. Please provide a output folder file path. Current Path: {args['output_folder']}!")  
    

    return(args)

def rayleigh_parser():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""

    print('-----------------------------------------')
    print('Start parsing command line options...')
    print('-----------------------------------------')
    
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
                        type = int, nargs = '?', default = 100, 
                        help = 'The dots per inch (dpi) resolution of the exported figures. Defaults to 100 ')

    parser.add_argument('--use_lin_scale', metavar = 'use_lin_scale',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, a linear scale will be used for the z axis (signal) ')

    parser.add_argument('--use_distance', metavar = 'use_distance',
                        type = bool, default = False, 
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

    parser.add_argument('-c', '--channels', metavar = 'channels',
                        type = str, nargs = '+', default = None, 
                        help = 'Type one or more channel names (e.g. xpar0355) here in order to open the figures in interactive mode ')

    parser.add_argument('--reference_height', metavar = 'reference_height',
                        type = float, nargs = '?', default = 9., 
                        help = 'The reference height/distance where the signals will be normalized for the Rayleigh fit. If use_distance is called, the limits correspond to distance. Defaults to 9 km ')

    parser.add_argument('--half_reference_window', metavar = 'half_reference_window',
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
    
    mandatory_args = ['input_file']

    mandatory_args_abr = ['-i']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            print(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')
            sys.exit('-- Program stopped')            
    
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..',
                                'plots','ray')
        args['output_folder'] = out_path
        os.makedirs(out_path, exist_ok = True)

    print("-- The following values have been used!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_file']):
        sys.exit(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}")  
    
    if os.path.basename(args['input_file'])[:3] != 'ray':
        sys.exit('---- Error: Measurement filename not understood! Please start the filename with ray (rayleigh fit)')
    
    if not os.path.exists(args['output_folder']):
        sys.exit(f"-- Error: The path to the output folder does not exists. Please provide a output folder file path. Current Path: {args['output_folder']}")  
    

    return(args)

def telecover_parser():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""

    print('-----------------------------------------')
    print('Start parsing command line options...')
    print('-----------------------------------------')
    
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
                        type = int, nargs = '?', default = 100, 
                        help = 'The dots per inch (dpi) resolution of the exported figures. Defaults to 100 ')

    parser.add_argument('--use_distance', metavar = 'use_distance',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, the y axis of the quicklook will correspond to the distance between the laser pulse and the telescope (vertical range) ')

    parser.add_argument('--y_lims', metavar = 'y_lims',
                        type = float, nargs = 2, default = [None, None], 
                        help = 'The y axis limits (lower and upper) of the normalized RC signal. Defaults to 0 (lower) 1.2 (upper) when use_lin_scale is True. If use_lin_scale is true then the lower limit becomes 1E-5 ')

    parser.add_argument('--x_lims', metavar = 'x_lims',
                        type = float, nargs = 2, default = [0., 3.], 
                        help = 'The x axis limits in km (lower and upper). If use_distance is called, the limits correspond to distance. Defaults to 0 km (lower) and 14 km (upper) If values below 0 or above the maximum signal altitude/distance are used, they will be ignored')

    parser.add_argument('--x_tick', metavar = 'x_tick',
                        type = int, nargs = '?', default = 0.25, 
                        help = 'The x axis finest tick in km. Defaults to 1km ')

    parser.add_argument('-c', '--channels', metavar = 'channels',
                        type = str, nargs = '+', default = None, 
                        help = 'Type one or more channel names (e.g. xpar0355) here in order to open the figures in interactive mode ')

    parser.add_argument('--normalization_height', metavar = 'normalization_height',
                        type = float, nargs = '?', default = 2., 
                        help = 'The reference height/distance where the signals will be normalized for the Rayleigh fit. If use_distance is called, the limits correspond to distance. Defaults to 9 km ')

    parser.add_argument('--half_normalization_window', metavar = 'half_normalization_window',
                        type = float, nargs = '?', default = 200., 
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
                        type = float, nargs = 2, default = [0.2, 14.], 
                        help = 'First and last altitude/distance boundaries (depending on the selection of use_distance) where smoothing will be applied, in km. Defines the signal region to be smoothed. If they exceed the current signal boundaries the existing boundaries will be used instead ')

    parser.add_argument('--half_window', metavar = 'half_window',
                        type = float, nargs = 2, default = [5., 200.], 
                        help = 'Half smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the same value twice to apply a constant window ')

    args = vars(parser.parse_args())
    
    mandatory_args = ['input_file']

    mandatory_args_abr = ['-i']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            print(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')
            sys.exit('-- Program stopped')            
    
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..',
                                'plots','tlc')
        args['output_folder'] = out_path
        os.makedirs(out_path, exist_ok = True)           
        
    print("-- The following values have been used!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_file']):
        sys.exit(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}")  
    
    if os.path.basename(args['input_file'])[:3] != 'tlc':
        sys.exit('---- Error: Measurement filename not understood! Please start the filename with tlc (telecover test)')
    
    if not os.path.exists(args['output_folder']):
        sys.exit(f"-- Error: The path to the output folder does not exists. Please provide a output folder file path. Current Path: {args['output_folder']}")  
    

    return(args)

def polarization_calibration_parser():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""

    print('-----------------------------------------')
    print('Start parsing command line options...')
    print('-----------------------------------------')
    
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
                        type = int, nargs = '?', default = 100, 
                        help = 'The dots per inch (dpi) resolution of the exported figures. Defaults to 100 ')

    parser.add_argument('--use_lin_scale', metavar = 'use_lin_scale',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called, a linear scale will be used for the z axis (signal) ')

    parser.add_argument('--use_distance', metavar = 'use_distance',
                        type = bool, default = False, 
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

    parser.add_argument('--ch_r', metavar = 'ch_r',
                        type = str, nargs = '+',
                        help = 'Type one or more channel names (e.g. xpar0355) here that correspond to the reflected channels used for the polarization calibration calculation. The number of reflected channels must be the same as the number of the respective transmitted channels ')

    parser.add_argument('--ch_t', metavar = 'ch_t',
                        type = str, nargs = '+', 
                        help = 'Type one or more channel names (e.g. xpat0355) here that correspond to the transmitted channels used for the polarization calibration calculation. The number of transmitted channels must be the same as the number of the respective reflected channels ')

    parser.add_argument('--calibration_height', metavar = 'calibration_height',
                        type = float, nargs = '?', default = 9., 
                        help = 'The calibration height/distance where the signals will be normalized for the Rayleigh fit. If use_distance is called, the limits correspond to distance. Defaults to 9 km ')

    parser.add_argument('--half_calibration_window', metavar = 'half_calibration_window',
                        type = float, nargs = '?', default = 500., 
                        help = 'The half window in meters used for the calibration. Defaults to 100 m ')

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

    parser.add_argument('-K', "--K", metavar = 'K',
                        type = float, nargs = '?', default = None, 
                        help = 'The K value for each channel pair. Defaults to 1 for all channels ')

    args = vars(parser.parse_args())
    
    mandatory_args = ['input_file','ch_r','ch_t']

    mandatory_args_abr = ['-i','-ch_r','ch_t']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            print(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')
            sys.exit('-- Program stopped')            
    
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..',
                                'plots','pcl')
        args['output_folder'] = out_path
        os.makedirs(out_path, exist_ok = True)
        
    print("-- The following values have been used!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_file']):
        sys.exit(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}")  
    
    if os.path.basename(args['input_file'])[:3] != 'pcl':
        sys.exit('---- Error: Measurement filename not understood! Please start the filename with ray (rayleigh fit)')
    
    if not os.path.exists(args['output_folder']):
        sys.exit(f"-- Error: The path to the output folder does not exists. Please provide a output folder file path. Current Path: {args['output_folder']}")  

    if len(args['ch_r']) != len(args['ch_t']):
        sys.exit('---- Error: The number of reflected channels is different from the number of transmitted channels! Please provide pairs of trasmitted and reflected channels')
        

    return(args)

def intercomparison_parser():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""

    print('-----------------------------------------')
    print('Start parsing command line options...')
    print('-----------------------------------------')
    
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
                        type = bool, default = False, 
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
    
    mandatory_args = ['input_files','output_folder']

    mandatory_args_abr = ['-i','-o']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            sys.exit(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')            

    print("-- The following values have been used!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_files'][0]):
        sys.exit(f"-- Error: The path to the first input file does not exists. Please provide a valid input file path. Current Path: {args['input_file'][0]}")  

    if not os.path.exists(args['input_files'][1]):
        sys.exit(f"-- Error: The path to the second input file does not exists. Please provide a valid input file path. Current Path: {args['input_file'][1]}")  
        
    if os.path.basename(args['input_files'][0])[:3] != 'ray' or \
        os.path.basename(args['input_files'][1])[:3] != 'ray':
        sys.exit('---- Error: Measurement filename not understood! Please start the filename with ray (rayleigh fit)')
    
    if not os.path.exists(args['output_folder']):
        sys.exit(f"-- Error: The path to the output folder does not exists. Please provide a output folder file path. Current Path: {args['output_folder']}")  

    if len(args['channels_1']) != len(args['channels_2']):
        sys.exit("-- Error: The number of channels provided in channels_1 is diffirent than the number of channels provided in channels_2. Please provide pair of channels that are going to be intercompared")          

    return(args)

def check_channels(sel_channels, all_channels):
    
    if not isinstance(sel_channels,type(None)):
        sel_channels = np.array(sel_channels)
        all_channels = np.array(all_channels)

        missing_ch = [ch not in all_channels for ch in sel_channels]
        
        if any(missing_ch):
            raise Exception("-- Error: The following provided channels do not exist: "+\
                            f"{sel_channels[missing_ch]} \n Please select one of:"+\
                            f"{all_channels}")
                
        channels = sel_channels
    
    else:
        channels = all_channels

    return(channels)

    