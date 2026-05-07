#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 18:20:12 2025

@author: nikos
"""

import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime as dt
from datetime import timedelta
import xarray as xr
from scipy.io import loadmat
from utils.error_classes import FileReaderError
from utils.time_conversions import datetimes_to_iso
from utils.error_classes import CustomWarning

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
        
        mfiles = glob.glob(os.path.join(dir_meas,'*.mat'))
                
        # for existing directory and files inside it, starts the reading of files     
        if len(mfiles) > 0:
            print(f'-- Reading {len(mfiles)} file(s)!')
            
            data = loadmat(mfiles[0])
            
            # Reading the licel file metadatas (header) - only for the first file
            system_info = read_meas(buffer = data['header'])

            channel_info = read_channels(buffer = data['datasetinfo'])
            
            channels = channel_info.index.values
            
            # Add the repetion rate, that was part of system info, to channel_info
            for ch in channels:
                if channel_info.loc[ch,"laser"] == "1" and system_info["laser_A_repetition_rate"] != None:
                    channel_info.loc[ch,"laser_repetition_rate"]  = system_info["laser_A_repetition_rate"]
                if channel_info.loc[ch,"laser"] == "2" and system_info["laser_B_repetition_rate"] != None:
                    channel_info.loc[ch,"laser_repetition_rate"]  = system_info["laser_B_repetition_rate"]
                if channel_info.loc[ch,"laser"] == "3" and system_info["laser_C_repetition_rate"] != None:
                    channel_info.loc[ch,"laser_repetition_rate"]  = system_info["laser_C_repetition_rate"]

            # bins_arr = np.arange(1., channel_info.bins.max() + 1.)
            bins_arr = np.arange(0., channel_info.bins.max())

            # Creating empty signal, shots, and time arrays
            start_time_arr = np.nan*np.zeros(len(mfiles), dtype = object)
            end_time_arr = np.nan*np.zeros(len(mfiles), dtype = object)

            shots_arr = np.nan*np.zeros((len(mfiles), len(channels)), dtype = object)
            sig_arr = np.nan*np.zeros((len(mfiles), len(channels), len(bins_arr)), dtype = float)

            filename = np.empty(len(mfiles), dtype = object)

            # Iterate over the files
            for k in range(len(mfiles)):
                
                filename[k] = os.path.basename(mfiles[k])
                
                data = loadmat(mfiles[k])
                                
                stime, etime = read_time(buffer = data['header'])

                shots_arr[k,:] = read_shots(buffer = data['header'])
                                
                dataKeys = [k for k in data.keys() if k.startswith('set')]
                body = np.array([data[dataKeys[i]] for i in range(len(dataKeys))])

                # Store signal, start and end time
                sig_arr[k, :, :] = body[:,0,:]

                start_time_arr[k] = stime
                
                if stime == etime: #only possible if the files are different by only milliseconds 
                    end_time_arr[k] = etime + timedelta(milliseconds = 500)
                    print(f'-- Warning! File {filename[k]} has the same start and end time reported (recording lasted < 1s). Please check it! ')
                else:
                    end_time_arr[k] = etime
            
            sig_raw = xr.DataArray(sig_arr, 
                                   coords=[start_time_arr, channels, bins_arr],
                                   dims=['time', 'channel', 'bins']) 
                        
            shots = xr.DataArray(shots_arr,  
                                 coords=[start_time_arr, channels],
                                 dims=['time', 'channel'])
            
            properties = ['filename', 'start_time', 'end_time']
    
            tdata = np.array([filename, 
                              datetimes_to_iso(start_time_arr), 
                              datetimes_to_iso(end_time_arr)], 
                             dtype = object)
            
            time_info = pd.DataFrame(tdata.T,  
                                     index = start_time_arr,
                                     columns = properties)  
                        
            # Sort by time
            sig_raw = sig_raw.sortby('time').copy()
            shots = shots.sortby('time').copy()
            time_info = time_info.sort_index()
            
            # Convert MHz to summed counts          
            sig_raw = unit_conv_MHz_to_counts(signal = sig_raw.copy(), shots = shots, channel_info = channel_info)

        else:
            CustomWarning(f"No files to read in: {dir_meas}") 

    return(system_info, channel_info, time_info, sig_raw, shots)


def read_meas(buffer):

    """ Retrieves location and geometry relevant information from 
    the licel header [altitude, latitude, longitude, 
    zenith angle, azimuth angle] and laser relevant information from 
    the licel header [laser A repetion rate, laser B repetion rate if it exists
    laser C repetion rate if it exists]"""
    system_info = pd.Series()
     
    metadata = buffer[2][0][0].split()

    system_info['station_altitude'] = float(metadata[5])    
    system_info['station_latitude'] = np.round(float(metadata[6]), 4)
    system_info['station_longitude'] = np.round(float(metadata[7]), 4)
    
    if len(metadata) > 8:
        system_info['zenith_angle'] = float(metadata[8])

    if len(metadata) > 9:
        system_info['azimuth_angle'] = float(metadata[9])

    metadata = buffer[3][0][0].split()

    system_info['laser_A_repetition_rate'] = float(metadata[2])

    system_info['laser_B_repetition_rate'] = np.nan
        
    system_info['laser_C_repetition_rate'] = np.nan
    
    return(system_info)

def read_time(buffer):
    
    """ Retrieves temporal information from 
    the licel header [start time, stop time]"""

    metadata = buffer[2][0][0].split()

    start_date = metadata[1]
    start_time = metadata[2]
    end_date = metadata[3]
    end_time = metadata[4]

    stime = dt.strptime(start_date + ' ' + start_time, "%d/%m/%Y %H:%M:%S") # start meas
    etime = dt.strptime(end_date + ' ' + end_time, "%d/%m/%Y %H:%M:%S") # start meas
        
    return(stime, etime)

def read_channels(buffer):
    
    """ Collects channel specific information from the licel header
    [analog/photon mode (0/1), laser number (A,B,C), number of range bins,
     laser polarization, high voltage, vertical resolution, 
     ADC range in mV (20,100,500), ADC bit used for the bit to mV conversion
     laser repetiotion rate, detected wavelength, channel polarization] """

    channel_info = pd.DataFrame()
       
    cols = ['active', 
            'acquisition_mode', 
            'unk1', 
            'bins', 
            'pmt_high_voltage', 
            'range_resolution', 
            'wave_pol', 
            'analog_to_digital_resolution', 
            'data_acquisition_range',
            'recorder_channel_id',
            'unk2',
            'unk3',
            'unk4']
    
    # Convert header to text, parse metadata
    header = np.array([buffer[i][0][0] for i in range(len(buffer))])

    # Header rows
    header = np.array([line.split()[:len(cols)] for line in header], 
                      dtype = object)

    arr_head = pd.DataFrame(header, columns = cols, dtype = object)
    
    arr_head['laser'] = 1

    # Combine from the recorder channel ID and the laser polarization    
    header_channel_id = list(arr_head.recorder_channel_id.values)
    if len(header_channel_id) == len(set(header_channel_id)):
        recorder_channel_id = header_channel_id
    else:
        raise FileReaderError("read_licel_matlab: Duplicate header_channel_id values")

    channel_info.index = recorder_channel_id

    info_columns = ['acquisition_mode', 'laser', 'bins', 
                    'range_resolution', 'data_acquisition_range', 
                    'analog_to_digital_resolution']
    
    channel_info.loc[:, 'recorder_channel_id'] = arr_head.loc[:, 'recorder_channel_id'].copy().values
    channel_info.loc[:, info_columns] = arr_head.loc[:, info_columns].copy().values.astype(object)

    mask_an = channel_info.loc[:,'acquisition_mode'].values == "0"
    channel_info.loc[mask_an, 'data_acquisition_range'] = (1000. * channel_info.loc[mask_an, 'data_acquisition_range'].astype(float))
    channel_info.loc[~mask_an, 'data_acquisition_range'] = None

    wave = np.array(list(np.char.split(arr_head.wave_pol.values.astype('str'),
                                       sep='.')))[:,0].astype(float)

    channel_info.loc[:,'detected_wavelength'] = wave
    
    channel_info.loc[:,'channel_bandwidth'] = 1. # default when not available in the raw files
    
    channel_info.loc[:,'dead_time_correction_type'] = 0. # default for Licel

    channel_info.loc[:,'bins'] = channel_info.loc[:,'bins'].astype(int) # convert to int

    channel_info["acquisition_mode"] = channel_info["acquisition_mode"].astype("object")
    channel_info.loc[mask_an,'acquisition_mode'] = "a" # convert to atlas nomenclature
    channel_info.loc[~mask_an,'acquisition_mode'] = "p" # convert to atlas nomenclature

    return(channel_info)

def read_shots(buffer):

    """ Retrieves the laser shot information"""
     
    metadata = buffer[3][0][0].split()

    shots = float(metadata[1])
    
    if shots == 0:
        shots = np.nan

    return(shots)

def read_buffer(fname):
       
    """ Reads the binary file as a single byte sequence (buffer)"""
    
    with open(fname, 'rb') as f:
        buffer = f.read()
        
    return(buffer)


def unit_conv_MHz_to_counts(signal, shots, channel_info):
    
    """
    General:
        Converts photon units from MHz to summed counts (licel default)
        
    Input:
        signal: 
            A 2D or 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, ...). The units of the photon 
            channels must be raw counts
            
        shots: 
            A 2D xarray with the laser shots per channel and timeframe.
            It should include the following dimensions: (time, channel, ...) 
            The index should correspond to the channel dimension of sig            

        channel_info: 
            A 2D pandas Dataframe that includes the channel metadata per channel
            
    Returns:
        
        sig_out: 
            An xarray in the same shape as sig. The units of the photon 
            channels are converted from MHz to raw counts
            
    """
    
    if len(signal) > 0:
        
        signal_out = signal.copy()
                        
        mask_pc = channel_info.acquisition_mode.values == 1
        
        channel_id_pc = channel_info.index.values[mask_pc]
                
        range_resolution = channel_info.range_resolution
            
        for ch in channel_id_pc:
    
            ch_d = dict(channel = ch)
                        
            bin_temporal_resolution = range_resolution.loc[ch] / 150. # in μs

            signal_out.loc[ch_d] = signal_out.loc[ch_d].copy() * shots.loc[:,ch] * bin_temporal_resolution 
    
    return(signal_out) 