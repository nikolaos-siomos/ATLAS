"""
@author: Peristera

Read the raw PollyXT data
"""
import glob
import pandas as pd
import numpy as np
import os
import xarray as xr
from readers.check_file_format import detect_netcdf
from utils.error_classes import FileReaderError
from utils.time_conversions import datetimes_to_iso

def read_dataset(dir_meas, meas_type = None):
    
    # Setting sig, info, and time as empty lists in the beggining    
    sig_raw = []     
    shots = []
    time_info = []
    
    system_info = []
    channel_info = []
    
    list_sig = []
    list_time = []
    list_shots = []
    
    if not(os.path.exists(dir_meas)):
        print('---- Warning : The folder for reading signals does not exist! '+\
              f'Check the input directory! \n Given folder: {dir_meas}')
    
    else:
        mfiles = glob.glob(os.path.join(dir_meas,'*.nc'))
        
        mfiles = [file for file in mfiles if os.path.basename(file) != 'temp.dat']
        
        if detect_netcdf(mfiles[0]) == None:
            raise FileReaderError(f"--QA test folder contains non netcdf files: {dir_meas}")
    
        # for existing directory and files inside it, starts the reading of files     
        if len(mfiles) > 0:
            
            print(f'-- Reading {len(mfiles)} file(s)!')
                
            raw_data = xr.open_dataset(mfiles[0])
            
            if "raw_signal" not in raw_data:
                raise FileReaderError(f"raw_signal parameter not found in the netcdf file. This is not a polly_xt raw file")
            
            # Reading the polly_xt file metadatas (header) - only for the first file            
            system_info = read_meas(raw_data = raw_data)            
            channel_info = read_channels(raw_data = raw_data)
            
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
            
            for k in range(len(mfiles)):
                                
                if detect_netcdf(mfiles[k]) == None:
                    raise FileReaderError(f"--The following file is not in netcdf format: {mfiles[k]}")

                raw_data = xr.open_dataset(mfiles[k])
                
                filename = np.empty(raw_data.time.size, dtype = object)
                
                filename[:] = os.path.basename(mfiles[k])
            
                raw_signal = raw_data.raw_signal.values.astype(float)
                raw_shots = raw_data.measurement_shots.values
                
                # Convert the time to npdatetime format and mask them
                start_time = raw_data.measurement_time

                start_time_arr = convert_time_to_npdatetime(start_time) 
                
                end_time_arr = start_time_arr + (start_time_arr[1] - start_time_arr[0])
            
                # Define signal xarray
                sig_f = xr.DataArray(raw_signal, 
                                   coords=[start_time_arr, bins_arr, channels],
                                   dims=['time', 'bins', 'channel'])
                                
                shots_f = xr.DataArray(raw_shots, 
                                       coords=[start_time_arr, channels],
                                       dims=['time', 'channel']) 

                # Define temporal data pandas Dataframe
                properties = ['filename', 'start_time', 'end_time']

                tdata = np.array([filename, 
                                  datetimes_to_iso(start_time_arr), 
                                  datetimes_to_iso(end_time_arr)], dtype = object)
                   
                time_info_f = pd.DataFrame(tdata.T,  
                                           index = start_time_arr,
                                           columns = properties)  
                    
                # Append the arrays to list in order to concatenate later
                list_sig.append(sig_f)
                list_shots.append(shots_f)                
                list_time.append(time_info_f)
        
            # Append in the time dimension all the time frames
            sig_raw = xr.concat(list_sig, dim='time')
            shots = xr.concat(list_shots, dim='time')
            time_info = pd.concat(list_time)
            
            # Transpose the dims in [time, channel, bins]            
            sig_raw = sig_raw.transpose('time','channel','bins')
            shots = shots.transpose('time','channel')
            
            # Sort by time
            sig_raw = sig_raw.sortby('time').copy()
            shots = shots.sortby('time').copy()
            time_info = time_info.sort_index() 

        else:
            print('---- Warning! No files to read in: {dir_meas}') 

    return(system_info, channel_info, time_info, sig_raw, shots)

def convert_time_to_npdatetime(time):
    # meas_time: 2d array : [times, 2], 2: [yyyymmdd, second of day]
    frames = time.shape[0]
    time_dt = []
    
    for i in range(frames): 
        temp_str = str(time.values[i][0])
        hr = np.divide(time.values[i][1],(60*60)); minute = (hr - np.fix(hr))*60; sec = (minute - np.fix(minute))*60
        date = temp_str[0:4]+'-'+temp_str[4:6]+'-'+temp_str[6:8]+'T'+"%02d" %np.fix(hr)+':'+"%02d" %np.fix(minute)+':'+"%02d" %np.fix(sec)
        time_dt.append(np.datetime64(date).astype("datetime64[ns]"))
    
    return(np.array(time_dt))

def read_meas(raw_data):

    """ Retrieves location and geometry relevant information from 
    the polly_xt metadata [altitude, latitude, longitude, 
    zenith angle, azimuth angle] and laser relevant information from 
    the licel header [laser A repetion rate, laser B repetion rate if it exists
    laser C repetion rate if it exists]"""
    
    system_info = pd.Series()
    
    system_info['altitude'] = raw_data.location_height.values   
    
    system_info['zenith_angle'] = raw_data.zenithangle.values
    system_info['zenith_angle'] = np.round(system_info['zenith_angle'], decimals = 2)
    
    system_info['azimuth_angle'] = 0.
    
    system_info['laser_A_repetition_rate'] = raw_data.laser_rep_rate.values
    
    system_info['laser_B_repetition_rate'] = np.nan
        
    system_info['laser_C_repetition_rate'] = np.nan

    return(system_info)

def read_channels(raw_data):
    channel_info = pd.DataFrame(index = (raw_data.channel.values + 1).astype(str))
    
    channel_info.loc[:,'acquisition_mode'] = "p" # only photon channels
    channel_info.loc[:,'laser'] = "1" # only 1 laser
    channel_info.loc[:,'bins'] = int(raw_data.height[-1].values + 1)
    channel_info.loc[:,'range_resolution'] = raw_data.measurement_height_resolution.values * 0.15 # nanosecond to meters
    channel_info.loc[:,'data_acquisition_range'] = None # analog channels only
    channel_info.loc[:,'analog_to_digital_resolution'] = None # analog channels only
    channel_info.loc[:,'detected_wavelength'] = raw_data.if_center.values
    channel_info.loc[:,'channel_bandwidth'] = raw_data.if_fwhm.values
    channel_info.loc[:,'dead_time_correction_type'] = 1. # default for PollyXTs

    
    return(channel_info)


            