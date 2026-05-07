#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 21:43:45 2025

@author: nikos
"""

import os, glob
import numpy as np
#from ..visualizer.writters import export_html
#from ..scc_converter.readers.read_config import config

def default_station_folder(parser_args):
    
    station_folder = os.path.join(parser_args['main_data_folder'],  parser_args['scc_station_id'])
    
    if not os.path.exists(station_folder):
        raise Exception(f"-- Error: No folder named after the provided station_id ({parser_args['station_id']}) was detected in the main_data_folder")

    return station_folder 

def find_parent_folders(station_folder):
    
    parent_folders = [os.path.basename(fld) for fld in glob.glob(os.path.join(station_folder,'*'), recursive = False)]
    if len(parent_folders) == 0:
        raise Exception("-- Error: No folders detected in the main_data_folder")

    return parent_folders

def autodetect_paths(parser_args):

    explicit_folder_check = any([parser_args.get('parent_folder') is None,
                                 parser_args.get('atlas_configuration_file') is None,
                                 parser_args.get('atlas_settings_file') is None,
                                 parser_args.get('radiosonde_folder') is None])
                                 
    if parser_args['main_data_folder'] is not None and explicit_folder_check:
                
        # Check if station id is provided. If not, station folder dis the same as the main_data_folder
        if parser_args.get('scc_station_id') is None:
            station_folder = parser_args['main_data_folder']
        else:
            station_folder = default_station_folder(parser_args)
            
        if parser_args.get('scc_configuration_id') is None:

            parent_folders = find_parent_folders(station_folder)
            config_ids = [fld.split('_')[2] for fld in parent_folders if len(fld.split('_')) >= 4 and fld.split('_')[0].isdigit() and fld.split('_')[1].isdigit() and fld.split('_')[2].isdigit()]
            lidar_ids = [fld.split('_')[0] for fld in parent_folders if len(fld.split('_')) >= 4 and fld.split('_')[0].isdigit() and fld.split('_')[1].isdigit() and fld.split('_')[2].isdigit()]
            
            unique_config_ids, ind_unique = np.unique(config_ids, return_index=True)
            unique_lidar_ids = np.array(lidar_ids)[ind_unique]
            
            if len(unique_config_ids) == 0:
                raise Exception("-- Error: The scc_configuration_id was not provided and cannot be infered from the parent folders. Please make sure that the parent folders are named correctly: <scc_lidar_id>_<scc_version_id>_<scc_configuration_id>_<data_identifier>")
            elif len(unique_config_ids) == 1:
                parser_args['scc_configuration_id'] = unique_config_ids[0]
                parser_args['scc_lidar_id'] = unique_lidar_ids[0]
            else:
                unique_config_id = universal_menu(
                    "Please select the SCC configuration ID from the list below: ",
                    choices=np.sort(unique_config_ids))
                print(f"Selected SCC configuration: {unique_config_id}")
                parser_args['scc_configuration_id'] = unique_config_id
                ind_config = np.where(unique_config_ids == unique_config_id)[0][0]
                parser_args['scc_lidar_id'] = unique_lidar_ids[ind_config]
        
       
        if parser_args.get('data_identifier') is None:
            parent_folders = [os.path.basename(fld) for fld in glob.glob(os.path.join(station_folder,f"*_*_{parser_args['scc_configuration_id']}_*"), recursive = False)]
            if len(parent_folders) == 0:
                raise Exception(f"-- Error: No parent folders detected matching the provided SCC configuration ID {parser_args['scc_configuration_id']}.")
        
            data_identifiers = ['_'.join(fld.split('_')[3:]) for fld in parent_folders if len(fld.split('_')) >= 4 and '.' not in fld]
            
            if len(data_identifiers) == 0:
                raise Exception("-- Error: It seems that a data_identifier was not used in the name of any parent folders. Please make sure that the parent folders are named correctly: <scc_lidar_id>_<scc_version_id>_<scc_configuration_id>_<data_identifier>")
            elif len(data_identifiers) == 1:
                parser_args['data_identifier'] = data_identifiers[0]
            else:
                data_identifier = universal_menu(
                    "Please select the date identifier from the list below: ",
                    choices=np.sort(data_identifiers))
                print(f"Selected identifier: {data_identifier}")
                parser_args['data_identifier'] = data_identifier
              
        if parser_args.get('parent_folder') is None:
            parent_folders = glob.glob(os.path.join(station_folder,f"*_*_{parser_args['scc_configuration_id']}_{parser_args['data_identifier']}"), recursive = False)
            
            parent_folder_name = os.path.basename(parent_folders[0])
            scc_lidar_id = parent_folder_name.split('_')[0]
            parser_args['scc_lidar_id'] = scc_lidar_id
            
            if len(parent_folders) == 0:
                raise Exception(f"-- Error: The provided SCC configuration ID ({parser_args['scc_configuration_id']}) and data_identifier ({parser_args['data_identifier']}) do not point to any valid parent folder")
            parser_args['parent_folder'] = parent_folders[0]
            
        if parser_args.get('atlas_configuration_file') is None and parser_args['export_hoi_cfg'] == '0':
	    
            config_files, expexted_patterns = infer_expected_paths(
                station_folder = station_folder, 
                sub_folder = "configurations", 
                label = "config_file", 
                config_id = parser_args['scc_configuration_id'], 
                data_identifier = parser_args['data_identifier']
                )
            if len(config_files) == 0:
                file_list = '\n'.join([f'config_file_{pat}.ini' for pat in expexted_patterns])
                raise Exception(f"-- Error: The ATLAS config file was not found. Expected path:\n{file_list}")
            else:
                parser_args['atlas_configuration_file'] = config_files[0]
                print(f"--Selected config file:\n{parser_args['atlas_configuration_file']}\n")

        elif parser_args.get('atlas_configuration_file') == None and parser_args['export_hoi_cfg'] in ['1','2']:
            parser_args['atlas_configuration_file'] = os.path.join(station_folder, "configurations", f"config_file_{parser_args['scc_configuration_id']}_exported.ini")
                                
        if parser_args.get('atlas_settings_file') is None:
            settings_files, expexted_patterns = infer_expected_paths(
                station_folder = station_folder, 
                sub_folder = "settings", 
                label = "settings_file", 
                config_id = parser_args['scc_configuration_id'], 
                data_identifier = parser_args['data_identifier']
                )
            if len(settings_files) == 0:
                file_list = '\n'.join([f'settings_file_{pat}.ini' for pat in expexted_patterns])
                print(f"-- Warning: The ATLAS settings file was not found. Expected path:\n{file_list}")
            else:
                parser_args['atlas_settings_file'] = settings_files[0]
                print(f"--Selected settings file:\n{parser_args['atlas_settings_file']}\n")
        
        if parser_args.get('radiosonde_folder') is None:
        
            radiosonde_folder_option_1 = os.path.join(station_folder, 'radiosondes', parser_args['scc_lidar_id'])
            radiosonde_folder_option_2 = os.path.join(station_folder, 'radiosondes')
            
            if os.path.exists(radiosonde_folder_option_1):
                radiosonde_folder = radiosonde_folder_option_1
            elif os.path.exists(radiosonde_folder_option_2):
                radiosonde_folder = radiosonde_folder_option_2
            else:
                raise Exception(f"-- Error: The radiosonde folder was not found. Expected one of the following paths:\n{radiosonde_folder_option_1}\n{radiosonde_folder_option_2}")
            
            parser_args['radiosonde_folder'] = radiosonde_folder
            
    return(parser_args)


def export_report(parser_args):
    
    print('-----------------------------------------')
    print('Creating QA report...')
    print('-----------------------------------------')

    path_cfg = parser_args['atlas_configuration_file']
    
    file_format = parser_args['file_format']
        
    # Reading of the configuration file    
    cfg = config(path = path_cfg, file_format = file_format, 
                 operation_mode = 'labeling') 
        
    # File format
    file_format = parser_args['file_format']
        
    if 'scc_configuration_id' in parser_args.keys():
        scc_configuration_id = cfg.system['configuration_id']
    else:
        scc_configuration_id = None
        
    if 'scc_station_id' in parser_args.keys():
        scc_station_id = cfg.system['station_id']
    else:
        scc_station_id = None

    if 'expert_analyst' in parser_args.keys():
        expert_analyst = parser_args['expert_analyst']
    else:
        expert_analyst = None

    if 'export_all' in parser_args.keys():
        export_all = parser_args['export_all']
    else:
        export_all = False

    data_identifier = parser_args['data_identifier']
  
    html_filename = export_html.make_filename(data_identifier = data_identifier, 
                                              expert_analyst = expert_analyst,
                                              scc_station_id = scc_station_id,
                                              scc_configuration_id = scc_configuration_id)

    photon_only = file_format in ['polly_xt', 'polly_xt_first']

    reports_folder = os.path.join(parser_args["output_folder"], 'reports')
    plots_folder = os.path.join(parser_args["output_folder"], 'plots')
    
    os.makedirs(reports_folder, exist_ok = True)

    export_html.QA_report(plots_folder = plots_folder, 
                          html_filename = os.path.join(reports_folder, html_filename),
                          photon_only = photon_only, 
                          export_all = export_all)
    

    return()

def universal_menu(prompt, choices):
    """
    Displays a numbered menu and returns the selected choice.
    Works in IPython, Jupyter, and any terminal.
    """
    print(prompt)
    for i, choice in enumerate(choices, start=1):
        print(f"{i}. {choice}")

    while True:
        selection = input("Enter number: ").strip()
        if selection.isdigit():
            idx = int(selection)
            if 1 <= idx <= len(choices):
                return choices[idx - 1]
        print("Invalid selection, try again.")
    
    return()
def infer_expected_paths(station_folder, sub_folder, label, config_id, data_identifier):
    
    #Oredered on selection priority
    suffix_dict = {}
    suffix_dict[0]  = ["*", "*", config_id, data_identifier]
    suffix_dict[1]  = ["*", "*", config_id, data_identifier,'*']
    suffix_dict[2]  = [config_id, data_identifier]
    suffix_dict[3]  = [config_id, data_identifier, '*']
    suffix_dict[4]  = ["*", "*", config_id]
    suffix_dict[5]  = ["*", "*", config_id, 'exported','*']
    suffix_dict[6]  = [config_id]
    suffix_dict[7]  = [config_id, 'exported','*']
    
    expected_suffixes = ["_".join(filter(None, suffix_dict[key])) \
                         for key in suffix_dict.keys()]

    expected_paths = [os.path.join(station_folder, sub_folder, f"{label}_{suf}.ini") for suf in expected_suffixes]

    # Set least priority for a general config_file.ini or settings_file.ini
    expected_paths = expected_paths + \
        [os.path.join(station_folder, sub_folder, f"{label}.ini")] 

    files = []
    for path in expected_paths:
        files.extend(glob.glob(path, recursive = False)) 

    
    return(files, expected_suffixes)
