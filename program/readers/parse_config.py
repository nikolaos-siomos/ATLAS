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

def parse_config():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""

    print('-----------------------------------------')
    print('Start parsing command line options...')
    print('-----------------------------------------')
    
    parser = argparse.ArgumentParser(
    	description='arguments ')
    
    parser.add_argument('-f', '--input_file', metavar = 'input_file', 
                        type = str, nargs = '?', 
                        help = 'The path to the input file ')

    parser.add_argument('-c', '--config_folder', metavar = 'config_folder', 
                        type = str, nargs = '?', 
                        help = 'The path to the configuration file that contains the necessary metadata. This optional argument can be used if the settings folder must be placed out of the parent_folder (default)')            

    parser.add_argument('-o', '--output_folder', metavar = 'output_folder', 
                        type = str, nargs = '?', 
                        help = 'The path to the output folder where the results and plots subfolders will be placed ')

    parser.add_argument('-d', '--debug', metavar = 'debug',
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then debugging files will be generated in ./results/debug folder. These included the metadata gathered from the configuration file, the licel header, and the combination of the two. Default to False ')

    parser.add_argument('-q', '--quicklook', metavar = 'quicklook', 
                        type = bool, default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then quiclook preprocessed files will be generated (original temporal resolution). ')
                                  
    parser.add_argument('-e', '--error_simulation', metavar = 'error_simulation', 
                        type = bool,  default = False, 
                        action = argparse.BooleanOptionalAction,
                        help = 'If called then a Monte Carlo Error simulation will be perform (warning! the processing time will increase proportionally to the number of iterations (100 by default). ')

    parser.add_argument('-s', '--signal_smoothing', metavar = 'signal_smoothing', 
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
    
    parser.add_argument('-i','--iterations', metavar = 'iterations', 
                        type = int, nargs = '?', default = 100,
                        help = 'The number of iterations for the Monte Carlo error simulation. Applicable only if -e is set to True. Default to 100 ')

    parser.add_argument('-t', '--timescale', metavar = 'timescale', 
                        type = int, nargs = '?', default = -1,
                        help = "The temporal window applied for the signal averaging in seconds. Not relevant for QA file preprocessing. If provided then the profiles will be averaged in as many timeframes fit among the measurement start time and end time. Select -1 (default) to average all the profiles. ")

    parser.add_argument('-v', '--vertical_limit', metavar = 'vertical_limit', 
                        type = float, nargs = '?', default = 20000.,
                        help = "The maximum altitude in m above which no calculations will be performed. Solar background calculations are performed prior to vertical signal trimming to enable background calculations up to the maximum signal altitude. Defaults to 20km ")                                                                                                        

    parser.add_argument('-w','--smoothing_window', metavar = 'smoothing_window', 
                        type = int, nargs = '?', default = 25,
                        help = 'Smoothing window in bins, common for all channels. Defaults to 25 bins ')

    parser.add_argument('--smoothing_sbin', metavar = 'smoothing_sbin', 
                        type = int, nargs = '?', default = -1,
                        help = 'The starting bin for smoothing. No smoothing will be applied before it. Defaults to -1 that corresponds to the first signal bin ')

    parser.add_argument('--smoothing_ebin', metavar = 'smoothing_ebin', 
                        type = int, nargs = '?', default = -1,
                        help = 'The ending bin for smoothing. No smoothing will be applied after it. Defaults to -1 that corresponds to the last signal bin ')

    args = vars(parser.parse_args())

    mandatory_args = ['input_file','config_folder','output_folder']

    mandatory_args_abr = ['-f','-c','-o']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            print(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')
            sys.exit('-- Program stopped')            
        
    print("-- The following values have been used!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")

    if not os.path.exists(args['input_file']):
        sys.exit(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}!")  
    
    if os.path.basename(args['input_file'])[:3] not in ['ray','tlc','pcl','drk']:
        sys.exit('---- Error: Measurement filename not understood! Please start the filename with ray (rayleigh), tlc (telecover), or pcl(polarization calibration) ')
    
    if not os.path.exists(args['config_folder']):
        sys.exit(f"-- Error: The path to the configuration folder does not exist. Please provide a valid input file path. Path: {args['config_file']}!")  

    if not os.path.exists(args['output_folder']):
        os.makedirs(args['output_folder'])

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
