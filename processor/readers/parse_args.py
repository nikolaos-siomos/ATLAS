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

def call_parser():
        
    """Collects the information included as commandline arguments. 
    """

    print('Parsing ATLAS Preprocessor arguments...')
    
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

    parser.add_argument('--vertical_limit', metavar = 'vertical_limit', 
                        type = float, nargs = '?', default = 20.,
                        help = "The maximum altitude in km above which no calculations will be performed. Solar background calculations are performed prior to vertical signal trimming to enable background calculations up to the maximum signal altitude. Defaults to 20km ")                                                                                                        
    
    # parser.add_argument('--timescale', metavar = 'timescale', 
    #                     type = int, nargs = '?', default = -1,
    #                     help = "The temporal window applied for the signal averaging in seconds. Not relevant for QA file preprocessing. If provided then the profiles will be averaged in as many timeframes fit among the measurement start time and end time. Select -1 (default) to average all the profiles. ")
    
    # parser.add_argument('--signal_smoothing', metavar = 'signal_smoothing', 
    #                     type = bool,  default = False, 
    #                     action = argparse.BooleanOptionalAction,
    #                     help = 'If called then a simple sliding average smoothing will be performed. Not relevant to QA file preprocessing')

    # parser.add_argument('--smoothing_window', metavar = 'smoothing_window', 
    #                     type = int, nargs = '?', default = 25,
    #                     help = 'Smoothing window in bins, common for all channels. Defaults to 25 bins ')

    # parser.add_argument('--smoothing_sbin', metavar = 'smoothing_sbin', 
    #                     type = int, nargs = '?', default = -1,
    #                     help = 'The starting bin for smoothing. No smoothing will be applied before it. Defaults to -1 that corresponds to the first signal bin ')

    # parser.add_argument('--smoothing_ebin', metavar = 'smoothing_ebin', 
    #                     type = int, nargs = '?', default = -1,
    #                     help = 'The ending bin for smoothing. No smoothing will be applied after it. Defaults to -1 that corresponds to the last signal bin ')

    args = vars(parser.parse_args())

    return(args)

def check_parser(args):
    
    mandatory_args = ['input_file']

    mandatory_args_abr = ['-i']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            raise Exception(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>') 
        
    print(" ")
    print("-- ATLAS Preprocessor arguments!")
    print("-------------------------------------------------------------------")
    for key in args.keys():
        print(f"{key} = {args[key]}")
    print("-------------------------------------------------------------------")
    print("")

    if not os.path.exists(args['input_file']):
        raise Exception(f"-- Error: The path to the input file does not exists. Please provide a valid input file path. Current Path: {args['input_file']}!")  
    
    if os.path.basename(args['input_file'])[:3] not in ['ray','tlc','pcl','drk']:
        raise Exception('---- Error: Measurement filename not understood! Please start the filename with ray (rayleigh), tlc (telecover), or pcl(polarization calibration) ')
    
    if args['output_folder'] == None:
        out_path = os.path.join(os.path.dirname(args['input_file']),'..','atlas_preprocessor')
        args['output_folder'] = out_path
    os.makedirs(args['output_folder'], exist_ok = True)
    
    return(args)

    