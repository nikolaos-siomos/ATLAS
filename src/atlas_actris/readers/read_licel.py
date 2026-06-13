import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime as dt
from datetime import timedelta
import xarray as xr
from readers.check_file_format import is_licel_header
from utils.error_classes import FileReaderError
from utils.time_conversions import datetimes_to_iso
from utils.error_classes import CustomWarning

def read_body(buffer, num_channels, bins):
    
    """ Reads the information from the raw licel files below the header.
    Blocks are separated by #"""

    sig_arr = np.nan * np.zeros((num_channels, np.max(bins)))

    buffer.readline(2)
    
    for i in range(num_channels):
        if i > 0:
            buffer.read(2)
        sig_arr[i,:bins[i]] = np.fromfile(buffer, dtype="<i4", count = bins[i])
        
    return(sig_arr)

def read_header(buffer):

    """ Retrieves location and geometry relevant information from 
    the licel header [altitude, latitude, longitude, 
    zenith angle, azimuth angle] and laser relevant information from 
    the licel header [laser A repetion rate, laser B repetion rate if it exists
    laser C repetion rate if it exists]"""
    system_info = pd.Series()
     
    # Just skip the first line
    buffer.readline()
    
    # Get the information from the 2nd line
    secondline = str(buffer.readline(), encoding="utf-8").split()
    file_format = 'new'
    
    # Check if the file has the old 2 rack format
    if len(secondline) == 1:
        buffer.readline()
        secondline = str(buffer.readline(), encoding="utf-8").split()  
        file_format = 'old'
    
    start_date = secondline[1]
    start_time = secondline[2]

    end_date = secondline[3]
    end_time = secondline[4]
    
    stime = dt.strptime(start_date + ' ' + start_time, "%d/%m/%Y %H:%M:%S") # start meas
    etime = dt.strptime(end_date + ' ' + end_time, "%d/%m/%Y %H:%M:%S") # end meas

    system_info['station_altitude'] = float(secondline[5])    
    system_info['station_latitude'] = np.round(float(secondline[6]), 4)
    system_info['station_longitude'] = np.round(float(secondline[7]), 4)
    
    if len(secondline) > 8:
        system_info['zenith_angle'] = float(secondline[8])

    if len(secondline) > 9:
        system_info['azimuth_angle'] = float(secondline[9])
        
    # Get the information from the 3nd line
    thirdline = str(buffer.readline(), encoding="utf-8").split()
    
    num_channels = int(thirdline[4])
    
    system_info['laser_A_repetition_rate'] = float(thirdline[1])

    if len(thirdline) > 2:
        system_info['laser_B_repetition_rate'] = float(thirdline[3])
    else:
        system_info['laser_B_repetition_rate'] = np.nan
        
    if len(thirdline) > 5:
        system_info['laser_C_repetition_rate'] = float(thirdline[6])
    else:
        system_info['laser_C_repetition_rate'] = np.nan
        
    return(system_info, stime, etime, num_channels, file_format)

def read_channels(buffer, num_channels, file_format):
    
    """ Collects channel specific information from the licel header
    [analog/photon mode (0/1), laser number (A,B,C), number of range bins,
     laser polarization, high voltage, vertical resolution, 
     ADC range in mV (20,100,500), ADC bit used for the bit to mV conversion
     laser repetiotion rate, detected wavelength, channel polarization] """
     
    if file_format == 'new':
        cols = [
            'active', 
            'acquisition_mode', 
            'laser_id', 
            'bins', 
            'laser_polarization', 
            'pmt_high_voltage', 
            'range_resolution', 
            'wave_pol', 
            'unk1', 
            'unk2', 
            'unk3', 
            'unk4', 
            'analog_to_digital_resolution', 
            'shots', 
            'data_acquisition_range',
            'recorder_channel_id'
            ]
    else:
        cols = [
            'active', 
            'acquisition_mode', 
            'laser_id', 
            'bins', 
            'laser_polarization', 
            'pmt_high_voltage', 
            'range_resolution', 
            'wave_pol', 
            'analog_to_digital_resolution', 
            'shots', 
            'data_acquisition_range',
            'recorder_channel_id',
            'unk1', 
            'unk2', 
            'unk3'
            ]                

    channel_info = pd.DataFrame()
    
    temp_info = []
    
    # Store channel metadata in a nested list - trailing columns will be ignored 
    for i in range(num_channels):
        linevars = str(buffer.readline(), encoding="utf-8").split()
        channel = linevars[cols.index('recorder_channel_id')]
        
        if channel.startswith('BT') or channel.startswith('BC'): 
            temp_info.append(linevars[:len(cols)])

    temp_info = pd.DataFrame(
        temp_info, 
        index = np.arange(len(temp_info)), 
        columns = cols
        )

    header_channel_id = temp_info.loc[:,'recorder_channel_id'].values
    laser_id = temp_info.loc[:,'laser_id'].values
    
    # Combine from the recorder channel ID and the laser polarization    
    if len(header_channel_id) == len(set(header_channel_id)):
        recorder_channel_id = header_channel_id
    else:
        recorder_channel_id = [(f'{ch_id}_L{lr_id}') 
                               for ch_id, lr_id in zip(header_channel_id, laser_id)]
        
    # Check if the defined channels are unique (unique sets of licel id and laser number)   
    if len(recorder_channel_id) != len(set(recorder_channel_id)):
        raise FileReaderError('-- Error: At least two of the licel channels have both the same id and laser number. Please correct this in the recorder settings')

    channel_info.index = recorder_channel_id

    info_columns = [
        'acquisition_mode', 
        'laser_id', 
        'bins', 
        'range_resolution', 
        'shots',
        'data_acquisition_range', 
        'analog_to_digital_resolution', 
        'recorder_channel_id'
        ]

    channel_info.loc[:, info_columns] = temp_info.loc[:, info_columns].copy().values.astype(object)

    mask_an = channel_info.loc[:,'acquisition_mode'].values == "0"

    channel_info.loc[mask_an, 'data_acquisition_range'] = (1000. * channel_info.loc[mask_an, 'data_acquisition_range'].astype(float))
    channel_info.loc[~mask_an, 'data_acquisition_range'] = None
    
    wave = np.array(list(np.char.split(temp_info.wave_pol.values.astype('str'),
                                       sep='.')))[:,0].astype(float)
    
    channel_info.loc[:,'detected_wavelength'] = wave

    channel_info.loc[:,'dead_time_correction_type'] = 0. # default for Licel

    channel_info["acquisition_mode"] = channel_info["acquisition_mode"].astype("object")
    channel_info.loc[mask_an,'acquisition_mode'] = "a" # convert to atlas nomenclature
    channel_info.loc[~mask_an,'acquisition_mode'] = "p" # convert to atlas nomenclature
    
    return channel_info

def read_buffer(fname):
       
    """ Reads the binary file as a single byte sequence (buffer)"""
    
    with open(fname, 'rb') as f:
        buffer = f.read()
        
    return(buffer)

def unit_conv_bits_to_mV(channel_info, signal, shots):

    """Converts analog signals from bits to mV"""
    
    if len(signal) > 0:
        
        mask_an = channel_info.acquisition_mode.values == "a"
        
        channel_id_an = channel_info.index.values[mask_an]
        
        data_acquisition_range = channel_info.data_acquisition_range.astype(float)
        
        analog_to_digital_resolution = channel_info.analog_to_digital_resolution.astype(float)
        
        for ch in channel_id_an:
            ch_d = dict(channel = ch)
            # analog conversion (to mV)
            signal.loc[ch_d] = signal.loc[ch_d]*data_acquisition_range.loc[ch]/(shots.loc[ch_d]*(np.power(2,analog_to_digital_resolution.loc[ch])-1.))
    
    return(signal) 

# Read measurement
def read_dataset(dir_meas, meas_type = None):
    
    """ Reads information from the raw licel files"""
    
    # Setting sig, info, and time as empty lists in the beggining    
    sig_raw = []     
    shots = []
    start_time_arr = []
    end_time_arr = []
    filename = []
    
    system_info = []
    channel_info = []
    time_info = []
    
    if not(os.path.exists(dir_meas)):
        print('---- Warning : The folder for reading signals does not exist! '+\
              f'Check the input directory! \n Given folder: {dir_meas}')
    
    else:
        
        mfiles = glob.glob(os.path.join(dir_meas,'*.*'))
        
        mfiles = [file for file in mfiles if os.path.basename(file) != 'temp.dat']
        
        
        # for existing directory and files inside it, starts the reading of files     
        if len(mfiles) > 0:
        
            print(f'-- Reading {len(mfiles)} file(s)!')
        
            if not is_licel_header(mfiles[0]):
                raise FileReaderError(f"--QA test folder contains non Licel files: {dir_meas}")
        
            buffer = open(mfiles[0], 'rb')
            
            # Reading the licel file metadatas (header) - only for the first file
            system_info, stime, etime, num_channels, file_format = \
                read_header(buffer)
            
            channel_info = read_channels(
                buffer = buffer, 
                num_channels = num_channels, 
                file_format = file_format
                )
            
            num_channels = channel_info.index.size
            
            channels = channel_info.index.values
            
            bins = channel_info.loc[:,"bins"].values.astype(int)

            # Add the repetion rate, that was part of system info, to channel_info
            for ch in channels:
                if channel_info.loc[ch,"laser_id"] == "1" and system_info["laser_A_repetition_rate"] != None:
                    channel_info.loc[ch,"laser_repetition_rate"] = system_info["laser_A_repetition_rate"]
                if channel_info.loc[ch,"laser_id"] == "2" and system_info["laser_B_repetition_rate"] != None:
                    channel_info.loc[ch,"laser_repetition_rate"] = system_info["laser_B_repetition_rate"]
                if channel_info.loc[ch,"laser_id"] == "3" and system_info["laser_C_repetition_rate"] != None:
                    channel_info.loc[ch,"laser_repetition_rate"] = system_info["laser_C_repetition_rate"]

            # bins_arr = np.arange(max(bins))
            bins_arr = np.arange(1, max(bins) + 1, 1)

            # Creating empty signal, shots, and time arrays
            start_time_arr = np.nan*np.zeros(len(mfiles), dtype = object)
            end_time_arr = np.nan*np.zeros(len(mfiles), dtype = object)

            shots_arr = np.nan*np.zeros((len(mfiles), len(channels)), dtype = object)
            sig_arr = np.nan*np.zeros((len(mfiles), len(channels), len(bins_arr)), dtype = float)

            filename = np.empty(len(mfiles), dtype = object)
                
            buffer.close()
            
            # Iterate over the files
            for k in range(len(mfiles)):

                # Store filemname
                filename[k] = os.path.basename(mfiles[k])
                
                # Check if it is a Licel file
                if not is_licel_header(mfiles[k]):
                    raise FileReaderError(f"--The following file is not in raw Licel format: {mfiles[k]}")
                
                buffer = open(mfiles[k], 'rb')
                
                # Reading the licel file metadatas (header) - only for the first file
                system_info_i, stime_i, etime_i, num_channels_i, file_format_i = \
                    read_header(buffer)
                
                # Check if the files have the same licel format (new or old)
                if file_format_i != file_format:
                    raise FileReaderError(f"--Not all files have the same Licel format: {mfiles[k]}\nCompare with: {mfiles[0]}")
                
                channel_info_i = read_channels(
                    buffer = buffer, 
                    num_channels = num_channels,
                    file_format = file_format_i
                    )

                # Check if the number of channel is the same for all files
                if not channel_info_i.index.equals(channel_info.index):
                    raise FileReaderError(f"--Not all files have the same channels: {mfiles[k]}\nCompare with: {mfiles[0]}")

                # Check if the number of bins is the same for all files per channel
                if not channel_info_i.loc[:,'bins'].equals(channel_info.loc[:,'bins']):
                    raise FileReaderError(f"--Not all files include channels with the same number of bins: {mfiles[k]}\nCompare with: {mfiles[0]}")
                
                shots_arr[k,:] = channel_info.loc[:,"shots"].values
                
                sig = read_body(
                    buffer = buffer, 
                    num_channels = num_channels, 
                    bins = bins
                    )

                buffer.close()

                # Store signal, start and end time
                sig_arr[k, :, :] = sig

                start_time_arr[k] = stime_i
                
                if stime >= etime: #only possible if the files are different by only milliseconds 
                    end_time_arr[k] = etime_i + timedelta(milliseconds = 500)
                    # print(f'-- Warning! File {filename[k]} has the same start and end time reported (recording lasted < 1s). Please check it! ')
                else:
                    end_time_arr[k] = etime_i

            sig_raw = xr.DataArray(sig_arr, 
                                   coords=[start_time_arr, channels, bins_arr],
                                   dims=['time', 'channel', 'bins']).astype(float)
            
            shots = xr.DataArray(shots_arr,  
                                 coords=[start_time_arr, channels],
                                 dims=['time', 'channel']).astype(float)

            tdata = np.array([filename, 
                              datetimes_to_iso(start_time_arr), 
                              datetimes_to_iso(end_time_arr)], 
                             dtype = object)
 
            properties = ['filename', 'start_time', 'end_time']
            
            time_info = pd.DataFrame(tdata.T,  
                                     index = start_time_arr,
                                     columns = properties,
                                     dtype = object)  

            # Remove non BT or BC channels
            valid_channels = [c for c in sig_raw.channel.values if c.startswith(("BT", "BC"))]
            sig_raw = sig_raw.sel(channel = valid_channels)
            shots = shots.sel(channel = valid_channels)
            channel_info = channel_info.loc[valid_channels,:]
                         
            # Sort by time
            sig_raw = sig_raw.sortby('time').copy()
            shots = shots.sortby('time').copy()
            time_info = time_info.sort_index()

            # Convert bits to mV relying solely on the raw file header
            sig_raw = unit_conv_bits_to_mV(signal = sig_raw.copy(), shots = shots, channel_info = channel_info)
            
        else:
            CustomWarning(f"No files to read in: {dir_meas}") 

    return(system_info, channel_info, time_info, sig_raw, shots)

