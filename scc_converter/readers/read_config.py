"""
@author: P. Paschou & N. Siomos
"""
import configparser
import re
import numpy as np
import pandas as pd
import sys

class config():
    
    def __init__(self, path, file_format):
        """Reads the config file at the given path"""
        
        parser = configparser.ConfigParser()
        parser.read(path, encoding="utf-8")

# Lidar
        if parser.has_section('Lidar'):

            self.meas = read_meas(parser['Lidar'], dtype = object)
            check_meas(meas_info = self.meas)

        else:
            
            raise Exception("-- Error: No lidar section is provided in the configuration files. Please include a section with at least the mandatory fields!")

# Channels
        if parser.has_section('Channels'):
            
            channel_section = read_channels(parser['Channels'], dtype = object)
            
            channel_section[channel_section == '_'] = np.nan

            check_channels(channel_info = channel_section, 
                           file_format = file_format)
            
            if 'recorder_channel_id' in channel_section.columns.values:
                if file_format == 'licel' or file_format == 'liel_matlab':
                    channels = [f'{channel_section.recorder_channel_id[i]}_L{str(int(channel_section.laser[i]))}'
                                for i in range(channel_section.index.size)]
                elif file_format == 'polly_xt':
                    channels = channel_section.recorder_channel_id.values

                channel_section.index = channels 
            
            self.channels = channel_section
            
        else:
            
            raise Exception("-- Error: No channel section is provided in the configuration files. Please include a section with at least the mandatory fields!")
            
# -------- END OF CLASS

def read_meas(section, dtype=object, skip_vars=[]):
    # Reads the whole or part of the section and returns a Pandas Series
    map_info = []
    
    val_list = []
    key_list = []
    
    for key in section:
        if key not in skip_vars:
            arr = [i.strip() for i in re.split(',', section[key]) if i !='']
            if len(arr) > 1:
                raise Exception(f'More than one values detected in the variable {key} of the configuration file. Please revise the {key} section of the configuration file!')
            elif len(arr) == 1:
                val_list.extend(arr)
                key_list.extend([key])
    
    if len(key_list) > 0:
        map_info = pd.Series(val_list, index = key_list)
    elif len(key_list) == 0:
        map_info = pd.Series()
    
    return(map_info)

def read_channels(section, dtype=object, skip_vars=[], squeeze = False):
    # Reads the whole or part of the section and returns a Pandas Series
    first = True
    map_info = []

    for key in section:
        if key not in skip_vars:
            arr = [i.strip() for i in re.split(',', section[key]) if i !='']
            if len(arr) > 0:
                temp = pd.DataFrame(arr, dtype = dtype, columns = [key])
                if first:
                    map_info = temp
                    first = False
                else:
                    check_len(map_info, temp, section, key)
                    map_info = pd.concat([map_info, temp], axis = 1)
    
    if len(map_info) > 0 and squeeze:
        map_info = map_info.squeeze()
    elif len(map_info) == 0:
        map_info = pd.DataFrame()
    
    return(map_info)

def comma_split(var, dtype):
    
    if var != '':
        var = re.split(',', var)
    
        var = np.array([item.strip() for item in var], 
                       dtype = dtype) #trimming the spaces
    else:
        var=[]
    return(var)

def check_channels(channel_info, file_format):

    if file_format == 'licel' or  file_format == 'licel_matlab':
            
        # recorder_channel_id check
        if 'recorder_channel_id' not in channel_info.columns.values:
            
            raise Exception("-- Error: The recorder_channel_id field is mandatory for licel systems! Please provide it in the configuration file. If this is not a licel system, please provide the correct file_format in the settings_file.")

        else:
            
            for recorder_channel_id in channel_info.recorder_channel_id:
                
                if isinstance(recorder_channel_id,str):
                    
                    if not (recorder_channel_id[:2] == 'BT' or recorder_channel_id[:2] == 'BC'):
                
                        raise Exception("-- Error: Provided recorder_channel_id not recognized. The first two letters must be either BT or BC. Please do not provide S2A ot S2P channels (currently not supported)!")
            
                else:
                    
                    raise Exception("-- Error: The recorder_channel_id provided in the configuration file must be a string. Please correct!")

        # Laser number check
        if 'laser' not in channel_info.columns.values:
            
            raise Exception("-- Error: The laser field is mandatory for licel systems! Please provide it in the configuration file")

        else:
            
            for laser in channel_info.laser:
                                
                if int(laser) not in [1, 2, 3, 4]:
            
                    raise Exception("-- Error: Provided laser number not recognized. Licel uses a laser number between 1 and 4!")

    if file_format == 'polly_xt':
            
        # recorder_channel_id check
        if 'recorder_channel_id' not in channel_info.columns.values:
            
            raise Exception("-- Error: Since version 0.4 the recorder_channel_id field is mandatory also for Polly XT systems! Please provide it in the configuration file. For Polly XTs use ascending integers starting from 1 and to the last channel following the order of channel from the raw netcdf files. ")

        else:
            
            for recorder_channel_id in channel_info.recorder_channel_id:
                
                if recorder_channel_id.isdigit() == False:
                    print(recorder_channel_id)
                    raise Exception("-- Error: The recorder_channel_id must be an intereger for PollyXT systems. Please revise the configuration file!")
            
                elif recorder_channel_id == 0:
                    
                    raise Exception("-- Error: The recorder_channel_id cnnot be zero. It must be an integer between 1 and up to the maximum number of available channels in the raw files. Please revise the configuration file!")


    mandatory = ['dead_time', 'trigger_delay_bins', 'telescope_type', 
                 'channel_type', 'channel_subtype']
    
    for mnd in mandatory:
        
        if mnd not in channel_info.columns.values:
        
            raise Exception(f"-- Error: The mandatory field {mnd} was not provided. Please include it in the configuration file!")

    # SCC ID check
    if 'scc_channel_id' not in channel_info.columns.values:
        channel_info.loc[:,'scc_channel_id'] = np.arange(0,channel_info.index.values.size,1)
        print('-- Warning: The scc_channel_id field is mandatory for the submission of the QA tests to CARS! Please make sure to include it if you plan to submit this set ') 
        
    # Dead time correction type check
    if 'dead_time_correction_type' not in channel_info.columns.values:
        
        channel_info.loc[:,'dead_time_correction_type'] = np.nan * np.zeros(channel_info.loc[:,'dead_time'].size, dtype = object)
        
        for i in range(channel_info.dead_time.size):
            if channel_info.dead_time[i] == channel_info.dead_time[i]:
                channel_info.loc[i,'dead_time_correction_type'] = '0'
            
    # Channel Bandwidth check
    if 'channel_bandwidth' not in channel_info.columns.values:
        
        channel_info.loc[:,'channel_bandwidth'] = '1.'
        
        print("-- Warning: The channel_bandwidth field was not provided. The default value (1 nm) has been used for all the channels. Please make sure that this corresponds to the actual system value")

    # Background Low
    if 'background_low_bin' not in channel_info.columns.values:
        print("-- Warning: The background_low_bin was not provided. For a pretriger range >= 400 (< 400) bins it will default to 100 bins (600 bins before the end of the profile). ")
        for i in range(channel_info.trigger_delay_bins.size):
            if int(channel_info.trigger_delay_bins[i]) <= -400:
                channel_info.loc[:,'background_low_bin'] = '100'

    # Background High
    if 'background_high_bin' not in channel_info.columns.values:
        print("-- Warning: The background_high_bin was not provided. For a pretriger range >= 400 (< 400) bins it will default to 100 bins before the end of the pretrigger range (end of the profile). ") 
        for i in range(channel_info.trigger_delay_bins.size):
            if int(channel_info.trigger_delay_bins[i]) <= -400:
                channel_info.loc[:,'background_high_bin'] = str(-int(channel_info.trigger_delay_bins[i]) - 100)
                
    return()
    
def check_meas(meas_info):

    mandatory = ['lidar_name', 'lidar_id']
    
    for mnd in mandatory:
        
        if mnd not in meas_info.index.values:
        
            raise Exception(f"-- Error: The mandatory field {mnd} was not provided. Please include it in the configuration file!")
            
    return()

def check_len(reference_var, testing_var, section, key):
    if len(reference_var) != len(testing_var):
        raise ValueError(f'Length inconsistencies detected in {section}, variable: {key}. All section variables must have the same length! Please revise the configuration file!')
    return()