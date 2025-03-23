#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 21:44:25 2025

@author: nikos
"""

import argparse, configparser, os, re
import numpy as np

def call_parser():
        
    """Collects the information included as commandline arguments. 
    """
        
    parser = argparse.ArgumentParser(
    	description='arguments ')
    

    parser.add_argument('-i', '--ini_file', metavar = 'ini_file', 
                        type = str, nargs = '?',  default = None,
                        help = 'The path to the initialization file of the calling script')

    args = vars(parser.parse_args())

    if args['ini_file'] == None:
        raise Exception("-- Error: The ini_file argument was not provided")
    
    if not os.path.exists(args['ini_file']):
        raise Exception(f"-- Error: The provided initialization file does not exist:\n{args['ini_file']}\nPlease provide a valid filepath")

    return(args)

def parse_ini_file(path):
    parser = configparser.ConfigParser()
    parser.optionxform = str
    parser.read(path, encoding="utf-8")
    
    parser_args = dict()
    for section in parser.sections():
        parser_args[section] = read_section(parser[section], 
                                            dtype = object, 
                                            squeeze = True)
    return(parser_args)

def read_section(section, dtype=object, skip_vars=[], squeeze = False, exception = []):
    # Reads the whole or part of the section and returns a Pandas Series
    map_info = dict()
    
    exceptions = {'slice_rayleigh' : str}

    for key in section:
        if key not in skip_vars:
            arr = [i.strip() for i in re.split(',', section[key]) if i !='']
            if len(arr) > 0:
                if key not in exceptions.keys():
                    arr = assume_type(arr)
                else:
                    arr = assume_type(arr, exception = exceptions[key])
                if squeeze and len(arr) == 1 and key not in exception:
                    map_info[key] = arr[0]
                else:
                    map_info[key] = arr
    
    return(map_info)

def assume_type(arr, exception = None):
    
    for i in range(len(arr)):
        if exception == None:
            if arr[i] == 'True':
                arr[i] = True
            elif arr[i] == 'False':
                arr[i] = False
            elif '.' in arr[i] and arr[i].replace('.', '', 1).replace('-', '', 1).isdigit() == True:
                arr[i] = float(arr[i])
            elif arr[i].isdigit():
                arr[i] = int(arr[i])       
        else:
           arr[i] = exception(arr[i])
        
    return(arr)

def comma_split(var, dtype):
    
    if var != '':
        var = re.split(',', var)
    
        var = np.array([item.strip() for item in var], 
                       dtype = dtype) #trimming the spaces
    else:
        var=[]
    return(var)

def check_parser_args(parser_args):
    
    mandatory_sections = ['autodetect_paths', 'options']
    
    mandatory_args = {'autodetect_paths' : ['main_data_folder'],
                      'options' : ['file_format']}
    
    path_bundle = {'explicit_paths' : ['parent_folder', 
                                       'atlas_configuration_file',
                                       'atlas_settings_file',
                                       'radiosonde_folder']}
    
    path_args = {'autodetect_paths' : ['main_data_folder'],
                 'explicit_paths' : ['parent_folder', 
                                     'atlas_configuration_file',
                                     'atlas_settings_file',
                                     'radiosonde_folder']}
        
    for section in mandatory_sections:
        if section not in parser_args.keys():
            raise Exception(f"-- Error: The mandatory section {section} was not found in the initialization file of the calling script")
            
    for section in mandatory_args.keys():
        for mandatory_arg in mandatory_args[section]:
            if mandatory_arg not in parser_args[section].keys():
                raise Exception(f"-- Error: The mandatory field {mandatory_arg} was not provided in the {section} section of the initialization file of the calling script")

    for section in path_bundle.keys():
        if len(parser_args[section]) != 0:
            if not all(key in parser_args[section].keys() for key in path_bundle[section]):
                raise Exception("-- Error: The explicit paths have been partially provided. Please either provide all of them together or none of them")

    for section in path_args.keys():
        for path_arg in path_args[section]:
            if path_arg in parser_args[section].keys():
                if not os.path.exists(parser_args[section][path_arg]):
                    raise Exception(f"-- Error: The provided {path_arg} path does not exist:\n{parser_args[section][path_arg]}")
        
    return()