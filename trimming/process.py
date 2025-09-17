#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 23:30:39 2022

@author: nick
"""

from ..readers.read_config import config
import os
import numpy as np
from ..readers import read_files
from ..tools import modify, make
from ..tools.automate import check_telecover_sec, check_telecover_rin
from ..tools.automate import detect_overflows
import shutil

def rayleigh(args, version):
   
    path_ray = args['rayleigh_folder']
    path_drk = args['dark_folder']
    path_cfg = args['config_file']
    
    file_format = args['file_format']
    operation_mode = args['operation_mode']
    
    # Reading of the configuration file    
    cfg = config(path = path_cfg, file_format = file_format, 
                 operation_mode = operation_mode) 
    
    print('-----------------------------------------')
    print(f'Lidar: {cfg.system.lidar_name}')
    print('-----------------------------------------')
    print(' ')

    # Read the files in the dark folder
    sig_raw_d, shots_d, system_info_d, channel_info_d, time_info_d  = \
        read_files.dark(finput_drk = path_drk, file_format = file_format)
    
    # Read the files in the rayleigh folder
    sig_raw, shots, system_info, channel_info, time_info = \
        read_files.rayleigh(finput_ray = path_ray, file_format = file_format)
        
    sig_raw, shots, time_info = \
        modify.slice_in_time(sig_raw = sig_raw, 
                             shots = shots, 
                             time_info = time_info,
                             slice_reg = args['slice_rayleigh'])
    
    if not isinstance(sig_raw,list):
        # Remove channels that should be excluded according to the configuration file
        if not isinstance(sig_raw_d,list):
            sig_raw_d, shots_d, channel_info_d, cfg = \
                modify.trim_channels(cfg = cfg, sig = sig_raw_d, shots = shots_d, channel_info = channel_info_d, meas_type = 'drk')
        
        sig_raw, shots, channel_info, cfg = \
            modify.trim_channels(cfg = cfg, sig = sig_raw, shots = shots, channel_info = channel_info, meas_type = 'ray')

        # Add the information from the raw file headers to the configuration object
        if not isinstance(sig_raw_d,list):
            cfg = modify.merge_config(cfg = cfg, system_info = system_info_d, channel_info = channel_info_d)
        cfg = modify.merge_config(cfg = cfg, system_info = system_info, channel_info = channel_info)
    
        # Add default values to the configuration object when the respective variables are not provided in the configuration file
        cfg = modify.fill_defaults(cfg)
    
        # Screen profiles that have an iregularly low number of shots
        if not isinstance(sig_raw_d,list):
            sig_raw_d = modify.screen_low_shots(time_info_d, channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
        sig_raw = modify.screen_low_shots(time_info, channel_info, signal = sig_raw.copy(), shots = shots)
    
        # Convert analog channel units to mV (applicable mainly to licel)   
        if not isinstance(sig_raw_d,list):
            sig_raw_d = modify.unit_conv_bits_to_mV(channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
        
        sig_raw = modify.unit_conv_bits_to_mV(channel_info, signal = sig_raw.copy(), shots = shots)
    
    
        # Detect and Screen Overflows
        sig_raw, shots, time_info = \
            detect_overflows(sig = sig_raw.copy(), 
                             shots = shots.copy(),
                             channel_info = cfg.channels,
                             time_info = time_info,
                             method = args['trim_overflows'],
                             meas_type = 'ray')
    
        # Detect and Screen Overflows for the dark
        if not isinstance(sig_raw_d,list):
            sig_raw_d, shots_d, time_info_d = \
                detect_overflows(sig = sig_raw_d.copy(), 
                                 shots = shots_d.copy(),
                                 channel_info = cfg.channels,
                                 time_info = time_info_d,
                                 method = args['trim_overflows'],
                                 meas_type = 'drk')    
                
        # Creating the measurement ID
        meas_ID = make.meas_id(lr_id = cfg.system['station_id'], time = sig_raw.time)
              
        # Creating the paths and folders
        nc_path = make.path(output_folder = args['output_folder'], system_info = cfg.system, time = sig_raw.time, meas_type = 'ray', version = version)
        
        # Checking for radiosonde data
        if os.path.isdir(args['radiosonde']):
            nc_path_rs = radiosonde(args, time = sig_raw.time.values, version = version, lidar = cfg.system['lidar_name'], lr_id = cfg.system['station_id'])
        else: 
            nc_path_rs = os.path.join(args['output_folder'], os.path.basename(args['radiosonde']))
            shutil.copy(args['radiosonde'], nc_path_rs)
        args['radiosonde'] = os.path.basename(nc_path_rs)
            
        # Making the raw SCC file
        make.rayleigh_file(system_info = cfg.system.copy(), 
                           channel_info = cfg.channels.copy(), 
                           time_info = time_info, time_info_d = time_info_d,
                           nc_path = nc_path, meas_ID = meas_ID,  
                           P = args['ground_pressure'], 
                           T = args['ground_temperature'], 
                           radiosonde_file = args['radiosonde'], 
                           sig = sig_raw, sig_d = sig_raw_d,
                           shots = shots, shots_d = shots_d)
    
        # Creating debugging files from the configuration and licel input
        if args['debug']:
            make.debug_file(path = args['output_folder'], data = cfg.system, meas_type = 'ray', label = 'config_lidar', meas_ID = meas_ID, show_index = True, header = False)
            make.debug_file(path = args['output_folder'], data = cfg.channels, meas_type = 'ray', label = 'config_channels', meas_ID = meas_ID)
            make.debug_file(path = args['output_folder'], data = time_info, meas_type = 'ray', label = 'time_info', meas_ID = meas_ID)
            if not isinstance(sig_raw_d,list):
                make.debug_file(path = args['output_folder'], data = time_info_d, meas_type = 'ray', label = 'time_info_d', meas_ID = meas_ID)
    
        
        print('Succesfully generated a rayleigh QA file!')
        print('')
    
    else:
        nc_path = None
        nc_path_rs = None
    
    return([nc_path, nc_path_rs])

def telecover(args, version):
    
    path_sec = args['telecover_sectors_folder']
    path_rin = args['telecover_rings_folder']
    path_drk = args['dark_folder']
    path_cfg = args['config_file']

    file_format = args['file_format']
    operation_mode = args['operation_mode']
    
    files_per_sector = args['files_per_sector']
    files_per_ring = args['files_per_ring']

    # Checking the telecover folder    
    check_telecover_sec(path_sec, files_per_sector = files_per_sector)
    check_telecover_rin(path_rin, files_per_ring = files_per_ring)
    
    # Reading of the configuration file    
    cfg = config(path = path_cfg, file_format = file_format, 
                 operation_mode = operation_mode)  
    
    print('-----------------------------------------')
    print(f'Lidar: {cfg.system.lidar_name}')
    print('-----------------------------------------')
    print(' ')

    # Read the files in the dark folder
    sig_raw_d, shots_d, system_info_d, channel_info_d, time_info_d  = \
        read_files.dark(finput_drk = path_drk, file_format = file_format)

    # Read the files in the telecover folder
    sig_raw, shots, system_info, channel_info, time_info = \
        read_files.telecover(finput_sec = path_sec, finput_rin = path_rin, file_format = file_format, files_per_sector = files_per_sector, files_per_ring = files_per_ring)

    if not isinstance(sig_raw,list):
        # Remove channels that should be excluded according to the configuration file
        if not isinstance(sig_raw_d,list):
            sig_raw_d, shots_d, channel_info_d, cfg = \
                modify.trim_channels(cfg = cfg, sig = sig_raw_d, shots = shots_d, channel_info = channel_info_d, meas_type = 'drk')
        
        sig_raw, shots, channel_info, cfg = \
            modify.trim_channels(cfg = cfg, sig = sig_raw, shots = shots, channel_info = channel_info, meas_type = 'tlc')
      
        # Add the information from the raw file headers to the configuration object
        if not isinstance(sig_raw_d,list):
            cfg = modify.merge_config(cfg = cfg, system_info = system_info_d, channel_info = channel_info_d)
        cfg = modify.merge_config(cfg = cfg, system_info = system_info, channel_info = channel_info)
    
        # Add default values to the configuration object when the respective variables are not provided in the configuration file
        cfg = modify.fill_defaults(cfg)
    
        # Screen profiles that have an iregularly low number of shots
        if not isinstance(sig_raw_d,list):
            sig_raw_d = modify.screen_low_shots(time_info_d, channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
        sig_raw = modify.screen_low_shots(time_info, channel_info, signal = sig_raw.copy(), shots = shots)
    
        # Convert analog channel units to mV (applicable mainly to licel)   
        if not isinstance(sig_raw_d,list):
            sig_raw_d = modify.unit_conv_bits_to_mV(channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
        
        sig_raw = modify.unit_conv_bits_to_mV(channel_info, signal = sig_raw.copy(), shots = shots)
    
    
        # Detect and Screen Overflows
        sig_raw, shots, time_info = \
            detect_overflows(sig = sig_raw.copy(), 
                             shots = shots.copy(),
                             channel_info = cfg.channels,
                             time_info = time_info,
                             method = args['trim_overflows'],
                             meas_type = 'tlc')
            
        # Detect and Screen Overflows for the dark
        if not isinstance(sig_raw_d,list):
            sig_raw_d, shots_d, time_info_d = \
                detect_overflows(sig = sig_raw_d.copy(), 
                                 shots = shots_d.copy(),
                                 channel_info = cfg.channels,
                                 time_info = time_info_d,
                                 method = args['trim_overflows'],
                                 meas_type = 'drk')    
        # Creating the measurement ID
        meas_ID = make.meas_id(lr_id = cfg.system['station_id'], time = sig_raw.time)
              
        # Creating the paths and folders
        nc_path = make.path(output_folder = args['output_folder'], system_info = cfg.system, time = sig_raw.time, meas_type = 'tlc', version = version)
        
        # Making the raw SCC file
        make.telecover_file(system_info = cfg.system.copy(), 
                            channel_info = cfg.channels.copy(),
                            time_info = time_info, time_info_d = time_info_d,
                            nc_path = nc_path, meas_ID = meas_ID, 
                            sig = sig_raw, sig_d = sig_raw_d,
                            shots = shots, shots_d = shots_d)
    
        # Creating debugging files from the configuration and licel input
        if args['debug']:
            make.debug_file(path = args['output_folder'], data = cfg.system, meas_type = 'tlc', label = 'config_lidar', meas_ID = meas_ID, show_index = True, header = False)
            make.debug_file(path = args['output_folder'], data = cfg.channels, meas_type = 'tlc', label = 'config_channels', meas_ID = meas_ID)
            make.debug_file(path = args['output_folder'], data = time_info, meas_type = 'tlc', label = 'time_info', meas_ID = meas_ID)
            if not isinstance(sig_raw_d,list):
                make.debug_file(path = args['output_folder'], data = time_info_d, meas_type = 'tlc', label = 'time_info_d', meas_ID = meas_ID)
        
        
        print('Succesfully generated a telecover QA file!')
        print('')
        
    else:
        nc_path = None
        
    return([nc_path])

def polarization_calibration(args, version):
   
    path_ray = args['rayleigh_folder']
    path_p45 = args['pcb_cal_p45_folder']
    path_m45 = args['pcb_cal_m45_folder']
    path_stc = args['pcb_cal_stc_folder']
    path_drk = args['dark_folder']
    path_cfg = args['config_file']

    file_format = args['file_format']
    operation_mode = args['operation_mode']
    
    # Reading of the configuration file    
    cfg = config(path = path_cfg, file_format = file_format, 
                 operation_mode = operation_mode) 
    
    print('-----------------------------------------')
    print(f'Lidar: {cfg.system.lidar_name}')
    print('-----------------------------------------')
    print(' ')

    # Check if the rayleigh filename was provided
    if not args['rayleigh_filename']:
        raise Exception("-- Error: A polarization calibration measurement is being processed but the rayleigh filename was not provided in the arguments! Please prepare the rayleigh file fist and included it with: -l <rayleigh_filename>'")
        
    # Read the files in the dark folder
    sig_raw_d, shots_d, system_info_d, channel_info_d, time_info_d  = \
        read_files.dark(finput_drk = path_drk, file_format = file_format)

    # Read the files in the calibration folder
    sig_raw, shots, system_info, channel_info, time_info = \
        read_files.polarization_calibration(finput_ray = path_ray, finput_p45 = path_p45, finput_m45 = path_m45, finput_stc = path_stc, file_format = file_format)

    if not isinstance(sig_raw,list):
        # Remove channels that should be excluded according to the configuration file
        if not isinstance(sig_raw_d,list):
            sig_raw_d, shots_d, channel_info_d, cfg = \
                modify.trim_channels(cfg = cfg, sig = sig_raw_d, shots = shots_d, channel_info = channel_info_d, meas_type = 'pcb')
        
        sig_raw, shots, channel_info, cfg = \
            modify.trim_channels(cfg = cfg, sig = sig_raw, shots = shots, channel_info = channel_info, meas_type = 'pcb')

        # Add the information from the raw file headers to the configuration object
        if not isinstance(sig_raw_d,list):
            cfg = modify.merge_config(cfg = cfg, system_info = system_info_d, channel_info = channel_info_d)
        cfg = modify.merge_config(cfg = cfg, system_info = system_info, channel_info = channel_info)
    
        # Add default values to the configuration object when the respective variables are not provided in the configuration file
        cfg = modify.fill_defaults(cfg)
    
        # Screen profiles that have an iregularly low number of shots
        if not isinstance(sig_raw_d,list):
            sig_raw_d = modify.screen_low_shots(time_info_d, channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
        sig_raw = modify.screen_low_shots(time_info, channel_info, signal = sig_raw.copy(), shots = shots)
    
        # Convert analog channel units to mV (applicable mainly to licel)   
        if not isinstance(sig_raw_d,list):
            sig_raw_d = modify.unit_conv_bits_to_mV(channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
        
        sig_raw = modify.unit_conv_bits_to_mV(channel_info, signal = sig_raw.copy(), shots = shots)
    
        # Detect and Screen Overflows
        sig_raw, shots, time_info = \
            detect_overflows(sig = sig_raw.copy(), 
                             shots = shots.copy(),
                             channel_info = cfg.channels,
                             time_info = time_info,
                             method = args['trim_overflows'],
                             meas_type = 'pcb')
            
        # Detect and Screen Overflows for the dark
        if not isinstance(sig_raw_d,list):
            sig_raw_d, shots_d, time_info_d = \
                detect_overflows(sig = sig_raw_d.copy(), 
                                 shots = shots_d.copy(),
                                 channel_info = cfg.channels,
                                 time_info = time_info_d,
                                 method = args['trim_overflows'],
                                 meas_type = 'drk')
                
        # Creating the measurement ID
        meas_ID = make.meas_id(lr_id = cfg.system['station_id'], time = sig_raw.time)
              
        # Creating the paths and folders
        nc_path = make.path(output_folder = args['output_folder'], system_info = cfg.system, time = sig_raw.time, meas_type = 'pcb', version = version)
        
        # Making the raw SCC file
        make.polarization_calibration_file(system_info = cfg.system.copy(), 
                                           channel_info = cfg.channels.copy(), 
                                           time_info = time_info, time_info_d = time_info_d,
                                           nc_path = nc_path, meas_ID = meas_ID, 
                                           P = args['ground_pressure'], 
                                           T = args['ground_temperature'], 
                                           radiosonde_file = args['radiosonde'],
                                           rayleigh = args['rayleigh_filename'],  
                                           sig = sig_raw, sig_d = sig_raw_d,
                                           shots = shots, shots_d = shots_d)
    
        # Creating debugging files from the configuration and licel input
        if args['debug']:
            make.debug_file(path = args['output_folder'], data = cfg.system, meas_type = 'pcb', label = 'config_lidar', meas_ID = meas_ID, show_index = True, header = False)
            make.debug_file(path = args['output_folder'], data = cfg.channels, meas_type = 'pcb', label = 'config_channels', meas_ID = meas_ID)
            make.debug_file(path = args['output_folder'], data = time_info, meas_type = 'pcb', label = 'time_info', meas_ID = meas_ID)
            if not isinstance(sig_raw_d,list):
                make.debug_file(path = args['output_folder'], data = time_info_d, meas_type = 'pcb', label = 'time_info_d', meas_ID = meas_ID)
                
    
        print('Succesfully generated a calibration QA file!')
        print('')

    else:
        nc_path = None

    
    return([nc_path])

def dark(args, version):

    path_drk = args['dark_folder']
    path_cfg = args['config_file']

    file_format = args['file_format']
    operation_mode = args['operation_mode']

    cfg = config(path = path_cfg, file_format = file_format, 
                 operation_mode = operation_mode) 
    
    print('-----------------------------------------')
    print(f'Lidar: {cfg.system.lidar_name}')
    print('-----------------------------------------')
    print(' ')

    # Read the files in the dark folder
    sig_raw_d, shots_d, system_info_d, channel_info_d, time_info_d  = \
        read_files.dark(finput_drk = path_drk, file_format = file_format)

    if not isinstance(sig_raw_d,list):
        # Remove channels that should be excluded according to the configuration file
        sig_raw_d, shots_d, channel_info_d, cfg = \
            modify.trim_channels(cfg = cfg, sig = sig_raw_d, shots = shots_d, channel_info = channel_info_d, meas_type = 'drk')
        
    
        # Add the information from the raw file headers to the configuration object
        cfg = modify.merge_config(cfg = cfg, system_info = system_info_d, channel_info = channel_info_d)
    
        # Add default values to the configuration object when the respective variables are not provided in the configuration file
        cfg = modify.fill_defaults(cfg)
        
        # Screen profiles that have an iregularly low number of shots
        sig_raw_d = modify.screen_low_shots(time_info_d, channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
    
        # Convert analog channel units to mV (applicable mainly to licel)   
        sig_raw_d = modify.unit_conv_bits_to_mV(channel_info_d, signal = sig_raw_d.copy(), shots = shots_d)
        
        # Detect and Screen Overflows for the dark
        sig_raw_d, shots_d, time_info_d = \
            detect_overflows(sig = sig_raw_d.copy(), 
                             shots = shots_d.copy(),
                             channel_info = cfg.channels,
                             time_info = time_info_d,
                             method = args['trim_overflows'],
                             meas_type = 'drk')
        
        # Creating the measurement ID
        meas_ID = make.meas_id(lr_id = cfg.system['station_id'], time = sig_raw_d.time)
              
        # Creating the paths and folders
        nc_path = make.path(output_folder = args['output_folder'], system_info = cfg.system, time = sig_raw_d.time, meas_type = 'drk', version = version)
        
        # Making the raw SCC file
        make.dark_file(system_info = cfg.system.copy(), 
                       channel_info = cfg.channels.copy(),
                       time_info = time_info_d,
                       nc_path = nc_path, meas_ID = meas_ID, 
                       sig = sig_raw_d, shots = shots_d)
    
        # Creating debugging files from the configuration and licel input
        if args['debug']:
            make.debug_file(path = args['output_folder'], data = cfg.system, meas_type = 'drk', label = 'config_lidar', meas_ID = meas_ID, show_index = True, header = False)
            make.debug_file(path = args['output_folder'], data = cfg.channels, meas_type = 'drk', label = 'config_channels', meas_ID = meas_ID)
            make.debug_file(path = args['output_folder'], data = time_info_d, meas_type = 'drk', label = 'time_info_d', meas_ID = meas_ID)
                
        
        print('Succesfully generated a dark QA file!')
        print('')
    
    else:
        nc_path = None
        
    return([nc_path])

def radiosonde(args, time, version, lidar, lr_id):
   
    path_rs = args['radiosonde']
    
    delimiter = args['rsonde_delimiter']
    skip_header = args['rsonde_skip_header']
    skip_footer = args['rsonde_skip_footer']
    usecols = args['rsonde_column_index']
    units = args['rsonde_column_units']
    
    mtime = np.datetime64(time[0] + (time[-1] - time[0]) / 2., 'us').item()
    
    try:
        date, time, atmo = read_files.radiosonde(path_rs, delimiter = delimiter, 
                                                 skip_header = skip_header, 
                                                 skip_footer = skip_footer, 
                                                 usecols = usecols,
                                                 units = units,
                                                 mtime = mtime,
                                                 ground = args['rsonde_altitude'])
    except:
        raise Exception('Error detecting when reading the radiosonde file. Please revise the radiosonde parsing options in the settings_file (converter section)')
    
    if args['rsonde_altitude'] == None:
       
        print('-- Warning: The rsonde_altitude field was not provided in the settings file. Substituting the empty value with the first atlitude bin of the radiosonde file')
    
        args['rsonde_altitude'] = atmo.height.values[0]
        
    # # Reading radiosonde geodata with the lidar station values if a geodata argument is not provided  
    # if any([geodata_i == None for geodata_i in geodata]):
    #     raise Exception("-- Error: The rsonde_geodata field is mandatory when processing a radiosonde file (mode = A and the radiosonde folder exists). Please provide 3 floats that correspond to the radiosonde station latitude, longitude, and altitude eg: --rsonde_geodata 40.5 22.9 60.0")

    # Creating the radiosonde ID
    rsonde_ID = f"{date}{lr_id}{time[:4]}"
    
    # Creating the paths and folders
    nc_path = make.rs_path(output_folder = args['output_folder'], meas_ID = rsonde_ID)
    
    # Making the raw SCC file
    make.radiosonde_file(nc_path = nc_path, date = date, time = time, 
                         ground = args['rsonde_altitude'], 
                         lat = args['rsonde_latitude'], 
                         lon = args['rsonde_longitude'], 
                         st_name = args['rsonde_station_name'], 
                         wmo_id = args['rsonde_wmo_number'], 
                         wban_id = args['rsonde_wban_number'], 
                         atmo = atmo)
    
    print('Succesfully generated a radiosonde file!')
    print('')
    
    return(nc_path)