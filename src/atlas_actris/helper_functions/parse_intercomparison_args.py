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
    
    print('Parsing master arguments...')    
    
    parser = argparse.ArgumentParser(
    	description='arguments ')
    

    parser.add_argument('-i', '--input_files', metavar='input_files', 
                        type=str, nargs=2,  default = None,
                        help='The paths to the files that will be compared')

    parser.add_argument('-o', '--output_folder', metavar='scc_converter', 
                        type=str, nargs='?', default = None,
                        help='The path to the results folder. This optional argument can be used if the results folder must be placed out of the parent_folder. Defaults to parent_folder/scc_converter ')

    parser.add_argument('-s', '--settings_file', metavar='settings_file', 
                        type=str, nargs='?', default = None,
                        help='The absolute path to the settings file that contains options to run ATLAS. ')            

    args = vars(parser.parse_args())
    
    return(args)

def check_parser(args):

    mandatory_args = ['input_files', 'output_folder', 'settings_file']

    mandatory_args_abr = ['-i','-o','-s']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            raise Exception(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')            
    
    for path in args['input_files']:
        if not os.path.exists(path):
            raise Exception(f"-- Error: The path to one of the input files does not exists. Please provide a valid input file path. Current Path: {args['input_file']}")  
        if '_ray_' not in os.path.basename(path):
            raise Exception('---- Error: Measurement filename not understood! The filename should contain the _ray_ field (rayleigh fit)')
    
        if not os.path.exists(args['output_folder']):
            raise Exception(f"-- Error: The path to output folder does not exists. Please provide a valid path. Current Path: {args['output_folder']}")  

        if not os.path.exists(args['settings_file']):
            raise Exception(f"-- Error: The path to settings_file does not exists. Please provide a valid path. Current Path: {args['settings_file']}")  

    return(args)


def substitute(org, rpl):
    
    for key in org.keys():
        if key in rpl.keys():
            if rpl[key] != None:
                org[key] = rpl[key]
                
    return(org)