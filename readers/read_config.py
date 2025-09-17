"""
@author: P. Paschou & N. Siomos
"""
import configparser
import re
import numpy as np
import pandas as pd
import sys

class config():
    
    def __init__(self, path, file_format, operation_mode):
        """Reads the config file at the given path"""
        
        parser = configparser.ConfigParser()
        parser.read(path, encoding="utf-8")

# Lidar
        if parser.has_section('System'):

            self.system = read_system(parser['System'], dtype = object)
            check_system(system_info = self.system, operation_mode = operation_mode)

        else:
            
            raise Exception("-- Error: No System section is provided in the configuration files. Please include a section with at least the mandatory fields!")

# Channels
        if parser.has_section('Channels'):
            
            channel_section = read_channels(parser['Channels'], dtype = object)
            
            check_channels(channel_info = channel_section, 
                           file_format = file_format,
                           operation_mode = operation_mode)
            
            if 'recorder_channel_id' in channel_section.columns.values:
                if file_format == 'licel' or file_format == 'licel_matlab' or file_format == 'licel_old2rack':
                    channels = [f'{channel_section.recorder_channel_id[i]}_L{str(int(channel_section.laser[i]))}'
                                for i in range(channel_section.index.size)]
                elif file_format == 'polly_xt' or file_format == 'polly_xt_first':
                    channels = channel_section.recorder_channel_id.values

                channel_section.index = channels 
            
            self.channels = channel_section
            
        else:
            
            raise Exception("-- Error: No channel section is provided in the configuration files. Please include a section with at least the mandatory fields!")
            
# -------- END OF CLASS

def read_system(section, dtype=object, skip_vars=[]):
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

def check_channels(channel_info, file_format, operation_mode):

    # Check the recorder_channel_id for licel systems
    if file_format == 'licel' or  file_format == 'licel_matlab':            
        if 'recorder_channel_id' not in channel_info.columns.values:
            raise Exception("-- Error: The recorder_channel_id field is mandatory for licel systems! Please provide it in the configuration file. If this is not a licel system, please provide the correct file_format in the settings_file.")

        else: 
            for recorder_channel_id in channel_info.recorder_channel_id:
                if isinstance(recorder_channel_id,str):
                    if not (recorder_channel_id[:2] == 'BT' or recorder_channel_id[:2] == 'BC'):
                        raise Exception("-- Error: Provided recorder_channel_id not recognized. The first two letters must be either BT or BC. Please do not provide S2A ot S2P channels (currently not supported)!")
                else: 
                    raise Exception("-- Error: The recorder_channel_id provided in the configuration file must be a string. Please correct!")

        # Check the laser number for licel systems
        if 'laser' not in channel_info.columns.values:
            raise Exception("-- Error: The laser field is mandatory for licel systems! Please provide it in the configuration file")
        else:
            for laser in channel_info.laser:
                if int(laser) not in [1, 2, 3]:
                    raise Exception("-- Error: Provided laser number not recognized. Licel uses a laser number between 1 and 3!")

    # Check the recorder_channel_id for PollyXTs
    if file_format == 'polly_xt' or file_format == 'polly_xt_first':
        if 'recorder_channel_id' not in channel_info.columns.values:
            raise Exception("-- Error: Since version 0.4 the recorder_channel_id field is mandatory also for Polly XT systems! Please provide it in the configuration file. For Polly XTs use ascending integers starting from 1 and to the last channel following the order of channel from the raw netcdf files. ")
        else:
            for recorder_channel_id in channel_info.recorder_channel_id:
                if recorder_channel_id.isdigit() == False:
                    raise Exception("-- Error: The recorder_channel_id must be an intereger for PollyXT systems. Please revise the configuration file!")
                elif recorder_channel_id == 0: 
                    raise Exception("-- Error: The recorder_channel_id cannot be zero. It must be an integer between 1 and up to the maximum number of available channels in the raw files. Please revise the configuration file!")

    # Check if all mandatory fields are provided
    if operation_mode == 'testing':
        mandatory = ['telescope_type', 'channel_type', 'channel_subtype']
    elif operation_mode == 'labeling':
        mandatory = ['scc_channel_id', 'telescope_type', 'channel_type', 'channel_subtype']
        
    
    for mnd in mandatory:
        if mnd not in channel_info.columns.values:
            raise Exception(f"-- Error: The mandatory field {mnd} was not provided. Please include it in the configuration file!")

    # Warn if any partially mandatory fields is not provided
    partially_optional = ['dead_time', 'daq_trigger_offset', 
                          'background_low_bin', 'background_high_bin']

    first = True
    for mnd in partially_optional:
        if mnd not in channel_info.columns.values:
            if first == True:
                print("-- Warning: The following partially optional fields were not provided in the configuration file. Default values will be used according to the manual. Please make sure that this corresponds to the actual system value:")
                first = False
            print(f"-- {mnd}")
    
    # Avoid high positive daq_trigger_offset
    if 'daq_trigger_offset' in channel_info.columns.values:
        if (channel_info.daq_trigger_offset.astype(float) > 50).any():
            raise Exception("-- Error: Unnaturally high positive daq_trigger_offset values (> 50 bins) were provided for one of the channels. Please note that for channels where the data aquisition starts BEFORE the Q-switch (pretriggering) the daq_trigger_offset must be negative! ")
    
    # Make sure the background_low_bin and background_high_bin are not set in the near range
    if ('background_low_bin' in channel_info.columns.values and 'background_high_bin' not in channel_info.columns.values) or \
        ('background_high_bin' in channel_info.columns.values and 'background_low_bin' not in channel_info.columns.values):
            raise Exception("-- Error: The background_low_bin and background_high_bin fields must be either provided together or not provided at all. Providing only one of them can lead to assigning unpredictable default values. ")

    if 'background_low_bin' in channel_info.columns.values:
        if (channel_info.background_high_bin.astype(float) - channel_info.background_low_bin.astype(float) < 100).any():
            raise Exception("-- Error: One of the background_high_bin values is less than 100 bins higher than the corresponding background_low_bin values. Please use a broader background correction range")
            
    if 'daq_trigger_offset' not in channel_info.columns.values and 'background_low_bin' in channel_info.columns.values:
        if ((channel_info.background_low_bin.astype(float) < 1000) & (channel_info.background_low_bin.astype(float) >= 0)).any():
            raise Exception("-- Error: One of the background_high_bin values is smaller than 1000 bins while the daq_trigger_offset is not provided. This will lead to a wrong background subraction ")
    elif 'daq_trigger_offset' in channel_info.columns.values  and 'background_low_bin' in channel_info.columns.values:
        if ((channel_info.daq_trigger_offset.astype(float) >= 0) & (channel_info.background_low_bin.astype(float) - channel_info.daq_trigger_offset.astype(float) < 1000) & (channel_info.background_low_bin.astype(float) >= 0)).any():
            raise Exception("-- Error: One of the background_low_bin values is less than 1000 bins above its corresponding absolute daq_trigger_offset. This will lead to a wrong background subraction ")
        if ((channel_info.daq_trigger_offset.astype(float) < 0) & (channel_info.background_low_bin.astype(float) + channel_info.daq_trigger_offset.astype(float) < 1000) & (channel_info.background_low_bin.astype(float) + channel_info.daq_trigger_offset.astype(float) >= 0) & (channel_info.background_low_bin.astype(float) >= 0)).any():             
            raise Exception("-- Error: One of the background_low_bin values is less than 1000 bins above its corresponding absolute daq_trigger_offset. This will lead to a wrong background subraction ")
        if ((channel_info.daq_trigger_offset.astype(float) < 0) & (channel_info.background_high_bin.astype(float) + channel_info.daq_trigger_offset.astype(float) < 1000) & (channel_info.background_high_bin.astype(float) + channel_info.daq_trigger_offset.astype(float) >= 0) & (channel_info.background_high_bin.astype(float) >= 0)).any():             
            raise Exception("-- Error: One of the background_high_bin values is less than 1000 bins above its corresponding absolute daq_trigger_offset. This will lead to a wrong background subraction ")
            
    # Check the field type
    allowed_telescope_types = ['n', 'f', 'x']
    for ch in channel_info.index:
        if channel_info.telescope_type.loc[ch] not in allowed_telescope_types:
            raise Exception(f"-- Error: The telescope_type '{channel_info.loc[ch,'telescope_type']}' of channel {channel_info.recorder_channel_id.loc[ch]} is not was not understood. Please revise the configuration file using values only among {allowed_telescope_types} for the telescope_type")
        
    # Check the channel type
    allowed_channel_types = ['p', 'c', 't', 'v', 'r', 'a', 'f']
    for ch in channel_info.index:
        if channel_info.channel_type.loc[ch] not in allowed_channel_types:
            raise Exception(f"-- Error: The channel_type '{channel_info.loc[ch,'channel_type']}' of channel {channel_info.recorder_channel_id.loc[ch]} was not understood. Please revise the configuration file using values only among {allowed_channel_types} for the channel_type")

    # Check the channel subtype
    allowed_channel_subtypes = ['r', 't', 'n', 'o', 'w', 'c', 'h', 'l', 'a', 'm', 'b', 's', 'x']
    for ch in channel_info.index:
        if channel_info.channel_subtype.loc[ch] not in allowed_channel_subtypes:
            raise Exception(f"-- Error: The channel_type '{channel_info.channel_subtype.loc[ch]}' of channel {channel_info.recorder_channel_id.loc[ch]} was not understood. Please revise the configuration file using values only among {allowed_channel_subtypes} for the channel_subtype")
          
    ### check for the subtype depending on the type
    allowed_channel_subtypes_per_type ={'p': ['r','t'],
                                        'c': ['r','t'],
                                        't': ['r','t','x'],
                                        'v': ['n', 'o', 'w', 'c'],
                                        'r': ['x','l','h'],
                                        'a': ['a', 'm'],
                                        'f': ['b', 's']}
    for ch in channel_info.index:
        if channel_info.loc[ch,'channel_subtype'] not in allowed_channel_subtypes_per_type[channel_info.loc[ch,'channel_type']]:
            raise Exception(f"-- Error: The channel_subtype '{channel_info.loc[ch,'channel_subtype']}' of channel {channel_info.recorder_channel_id.loc[ch]} is not applicable for the provided channel_type '{channel_info.loc[ch,'channel_type']}'. Please revise the configuration file using values only among {allowed_channel_subtypes_per_type[channel_info.loc[ch,'channel_type']]} for the channel_subtype")
   
    return()
    
def check_system(system_info, operation_mode):

    # Check if all mandatory fields are provided
    if operation_mode == 'testing':
        mandatory = ['lidar_name', 'station_id', 'station_name']
    else:
        mandatory = ['lidar_name', 'lidar_id', 'station_name', 'station_id', 
                     'version_name', 'version_id', 'configuration_name',
                     'configuration_id']
        
    
    for mnd in mandatory:
        if mnd not in system_info.index.values:
            raise Exception(f"-- Error: The mandatory field {mnd} was not provided. Please include it in the configuration file!")
            
    return()

def check_len(reference_var, testing_var, section, key):
    
    if len(reference_var) != len(testing_var):
        raise ValueError(f'Length inconsistencies detected in {section}, variable: {key}. All section variables must have the same length! Please revise the configuration file!')
    
    return()