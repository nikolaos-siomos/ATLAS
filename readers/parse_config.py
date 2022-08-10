"""
Created on Sun May 17 21:06:52 2020

@author: nick
"""

import argparse
import os
import sys

def parse_config():
        
    """Collects the information included as commandline arguments. 
    Current mandatory arguments: --parent_folder [-f]"""
    
    parser = argparse.ArgumentParser(
    	description='arguments ')
    
    parser.add_argument('-f', '--input_file', metavar='input_file', 
                        type=str, nargs='+', 
                        help='The path to the input file ')

    parser.add_argument('-c', '--config_folder', metavar='config_folder', 
                        type=str, nargs='+', 
                        help='The path to the configuration file that contains the necessary metadata. This optional argument can be used if the settings folder must be placed out of the parent_folder (default)')            

    parser.add_argument('-o', '--output_folder', metavar='output_folder', 
                        type=str, nargs='+', 
                        help='The path to the output folder where the results and plots subfolders will be placed ')

    parser.add_argument('-d', '--debug', metavar='debug', 
                        type=bool, nargs='+', 
                        help='If set to True then debugging files will be generated in ./results/debug folder. These included the metadata gathered from the configuration file, the licel header, and the combination of the two. Default to False ')
        
    args = vars(parser.parse_args())

    mandatory_args = ['input_file','config_folder','output_folder']

    scalar_args = ['input_file', 'config_folder', 'output_folder', 'debug']
    
    mandatory_args_abr = ['-f','-c','-o']
    
    for i in range(len(mandatory_args)):
        if not args[mandatory_args[i]]:
            print(f'-- Error: The mandatory argument {mandatory_args[i]} is not provided! Please provide it with: {mandatory_args_abr[i]} <path>')
            sys.exit('-- Program stopped')            

    for i in range(len(scalar_args)):
        if isinstance(args[scalar_args[i]],list):
            if len(args[scalar_args[i]]) > 1:
                print(f'-- Error: More than one objects provided for argument {scalar_args[i]}!')
                sys.exit('-- Program stopped')
            if len(args[scalar_args[i]]) == 1:
                args[scalar_args[i]] = args[scalar_args[i]][0]    

    optional_args = ['debug']
    
    default_values = [False]
    
    print("-------------------------------------------------------------------")
    print("-- Warning: The following default values have been used!")
    print("-------------------------------------------------------------------")
    for i in range(len(optional_args)):
        if not args[optional_args[i]]:
            args[optional_args[i]] = default_values[i]
            print(f"{optional_args[i]} = {args[optional_args[i]]}")
    print("-------------------------------------------------------------------")
            
    if not os.path.exists(args['config_folder']):
        sys.exit(f"-- Error: Path to the configuration file does not exists (defaults to <parent_folder>/config_file.ini). Path: {args['config_file']}!")  

    if not os.path.exists(args['output_folder']):
        os.makedirs(args['output_folder'])

    if not os.path.exists(os.path.join(args['output_folder'],'plots')):
        os.makedirs(os.path.join(args['output_folder'],'plots'))

    if not os.path.exists(os.path.join(args['output_folder'],'results')):
        os.makedirs(os.path.join(args['output_folder'],'results'))
        
    if args['debug'] not in [True, False]:
        sys.exit(f"-- Error: debug field should be boolean. Please use one of {[True, False]} with: -d <debug>")
        
    if not os.path.exists(os.path.join(args['output_folder'],'debug')) and args['debug']:
        os.makedirs(os.path.join(args['output_folder'],'debug'))

    return(args)