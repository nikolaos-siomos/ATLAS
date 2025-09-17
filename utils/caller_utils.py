#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 21:43:45 2025

@author: nikos
"""

import os, glob
import numpy as np
# from visualizer.writters import export_html

def autodetect_paths(parser_args):

    explicit_folder_check = any([parser_args.get('parent_folder') == None,
                                 parser_args.get('atlas_configuration_file') == None,
                                 parser_args.get('atlas_settings_file') == None,
                                 parser_args.get('radiosonde_folder') == None])
                                 
    if parser_args['main_data_folder'] != None and explicit_folder_check:
        
        if parser_args.get('scc_station_id') != None:
            station_folder = os.path.join(parser_args['main_data_folder'],  parser_args['scc_station_id'])
            
            if not os.path.exists(station_folder):
                raise Exception(f"-- Error: No folder named after the provided station_id ({parser_args['station_id']}) was detected in the main_data_folder")
        else:
            station_folder = parser_args['main_data_folder']
            
        if parser_args.get('scc_configuration_id') == None:
            parent_folders = [os.path.basename(fld) for fld in glob.glob(os.path.join(station_folder,'*'), recursive = False)]
            if len(parent_folders) == 0:
                raise Exception("-- Error: No folders detected in the main_data_folder")
        
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
        
       
        if parser_args.get('data_identifier') == None:
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
              
        if parser_args.get('parent_folder') == None:
            parent_folders = glob.glob(os.path.join(station_folder,f"*_*_{parser_args['scc_configuration_id']}_{parser_args['data_identifier']}"), recursive = False)
            
            parent_folder_name = os.path.basename(parent_folders[0])
            scc_lidar_id = parent_folder_name.split('_')[0]
            parser_args['scc_lidar_id'] = scc_lidar_id
            
            if len(parent_folders) == 0:
                raise Exception(f"-- Error: The provided SCC configuration ID ({parser_args['scc_configuration_id']}) and data_identifier ({parser_args['data_identifier']}) do not point to any valid parent folder")
            parser_args['parent_folder'] = parent_folders[0]
            
        if parser_args.get('atlas_configuration_file') == None and parser_args['export_hoi_cfg'] == '0':
            expected_path = os.path.join(station_folder, 'configurations', f"config_file_*{parser_args['scc_configuration_id']}*.ini")
            config_files = glob.glob(expected_path, recursive = False)
            if len(config_files) == 0:
                raise Exception(f"-- Error: The ATLAS configuration file was not found. Expected path:\n{expected_path}")
            elif len(config_files) == 1:
                parser_args['atlas_configuration_file'] = config_files[0]
            else:
                print("More than one potential configuration files detected:")
                config_filenames = [os.path.basename(path) for path in config_files]
                config_filename = universal_menu(
                    "More than one matching configuration files were detected. Please select one file from the list below: ",
                    choices=np.sort(config_filenames))
                print(f"Selected configuration file: {config_filename}")
                selected_index = config_filenames.index(config_filename)
                parser_args['atlas_configuration_file'] = config_files[selected_index]
        elif parser_args.get('atlas_configuration_file') == None and parser_args['export_hoi_cfg'] in ['1','2']:
            parser_args['atlas_configuration_file'] = os.path.join(station_folder, "configurations", f"config_file_{parser_args['scc_configuration_id']}_exported.ini")
                                
                
        if parser_args.get('atlas_settings_file') == None:
            expected_path = os.path.join(station_folder, 'settings', f"settings_file_*{parser_args['scc_configuration_id']}*.ini")
            settings_files = glob.glob(expected_path, recursive = False)
            if len(settings_files) == 0:
                print(f"-- Warning: The ATLAS settings file was not found. Expected path:\n{expected_path}")
            elif len(settings_files) == 1:
                parser_args['atlas_settings_file'] = settings_files[0]
            else:
                print("More than one potential settings files detected:")
                settings_filenames = [os.path.basename(path) for path in settings_files]
                settings_filename = universal_menu(
                    "More than one matching configuration files were detected. Please select one file from the list below: ",
                    choices=np.sort(settings_filenames))
                print(f"Selected settings file: {settings_filename}")
                selected_index = settings_filenames.index(settings_filename)
                parser_args['atlas_settings_file'] = settings_files[selected_index]

        if parser_args.get('radiosonde_folder') == None:
        
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

# def prepare_master_args(parser_args):
    
#     # Parent folder
#     parent_folder = parser_args['parent_folder']
    
#     mst_args = parse_mst(get_defaults = True)

#     mst_args['parent_folder'] = parser_args['parent_folder']
#     mst_args['config_file']   = parser_args['atlas_configuration_file']
#     mst_args['settings_file'] = parser_args['atlas_settings_file']
#     mst_args['radiosonde']    = parser_args['radiosonde_folder']
#     mst_args['file_format']   = parser_args['file_format']

#     if parser_args.get('nrm') != None:
#         nrm = parser_args['nrm']
#         mst_args['rayleigh_folder'] = os.path.join(parent_folder, nrm)

#     if parser_args.get('pcb') != None:
#         pcb = parser_args['pcb']
#         mst_args['pcb_cal_p45_folder'] = os.path.join(parent_folder, pcb, '+45')
#         mst_args['pcb_cal_m45_folder'] = os.path.join(parent_folder, pcb, '-45')

#     if parser_args.get('tlc') != None:
#         tlc = parser_args['tlc']
#         mst_args['telecover_sectors_folder']  = os.path.join(parent_folder, tlc)
        
#     if parser_args.get('tlc_rin') != None:
#         tlc_rin = parser_args['tlc_rin']
#         mst_args['telecover_rings_folder'] = os.path.join(parent_folder, tlc_rin)
        
#     if parser_args.get('drk') != None:
#         drk = parser_args['drk']
#         mst_args['dark_folder'] = os.path.join(parent_folder, drk)

#     if parser_args.get('process') != None:
#         mst_args['process'] = parser_args['process']
        
#     if parser_args.get('process_qck') != None:
#         mst_args['process_qck'] = parser_args['process_qck']

#     if parser_args.get('quick_run') != None:
#         mst_args['quick_run'] = parser_args['quick_run']

#     if parser_args.get('slice_rayleigh') != None:
#         mst_args['slice_rayleigh'] = parser_args['slice_rayleigh']
        
#     mst_args = check_mst(mst_args)

#     return(mst_args)

# def export_report(parser_args):
    
#     # Parent folder
#     parent_folder = parser_args['parent_folder']

#     # File format
#     file_format = parser_args['file_format']
    
#     if 'scc_configuration_id' in parser_args.keys():
#         scc_configuration_id = parser_args['scc_configuration_id']
#     else:
#         scc_configuration_id = ''
        
#     if 'scc_station_id' in parser_args.keys():
#         scc_station_id = parser_args['scc_station_id']
#     else:
#         scc_station_id = ''
        
#     if 'expert_analyst' in parser_args.keys():
#         expert_analyst = parser_args['expert_analyst']
#     else:
#         expert_analyst = ''
        
#     if 'export_all' in parser_args.keys():
#         export_all = parser_args['export_all']
#     else:
#         export_all = False
        
#     data_identifier = parser_args['data_identifier']

#     html_filename = export_html.make_filename(data_identifier = data_identifier, 
#                                               expert_analyst = expert_analyst,
#                                               scc_station_id = scc_station_id,
#                                               scc_configuration_id = scc_configuration_id)

#     photon_only = file_format in ['polly_xt', 'polly_xt_first']

#     reports_folder = os.path.join(parent_folder, 'reports')
#     os.makedirs(reports_folder, exist_ok = True)

#     export_html.QA_report(plots_folder = os.path.join(parent_folder, 'plots'), 
#                           html_filename = os.path.join(reports_folder, html_filename),
#                           photon_only = photon_only, 
#                           export_all = export_all)
        
#     return()

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