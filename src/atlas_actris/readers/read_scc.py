"""
@authors: Siomos and Paschou

================
Input:
    arg 1: full filename from current running directory, should be a string scalar
    arg 2: the code (first letters) of the measurements of the lidar 
Returns:
    arg 1-2: the signals if xarray format with dimensions (time, channel, bins)
    arg 3-4: info from header per channel (e.g. pmt type, laser polarization, bins, resolution, shots, channel type(p,s, total), ADC_range, ADC_bit, wavelength)
    arg 5-6: the ground altitude and the measurement angle(off-zenith)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import glob
import datetime as dt
import xarray as xr
from readers.check_file_format import detect_netcdf
from utils.error_classes import FileReaderError
from utils.time_conversions import datetimes_to_iso
from utils.error_classes import CustomWarning

# Read measurement
def read_dataset(dir_meas: str, meas_type: str):
    
    """
    dir_meas: Measurement folder (can contain more than one netcdf files)
    meas_type: If set to drk then the dark signals will be extracted"""
    
    # Setting sig, info, and time as empty lists in the beggining    
    sig_raw = []     
    shots = []
    time_info = []
    
    system_info = []
    channel_info = []
    
    list_sig = []
    list_time = []
    list_shots = []
    
    if not(Path(dir_meas).exists()):
        CustomWarning("The folder for reading signals does not exist! "+\
              f"Check the input directory! \n Given folder: {dir_meas}")
    
    else:
        
        mfiles = [p for p in Path(dir_meas).glob("*.*") if p.is_file()]
    
        if len(mfiles) > 0:
            
            if detect_netcdf(mfiles[0]) == None:
                raise FileReaderError(f"--QA test folder contains non netcdf files: {dir_meas}")
        
            print(f'-- Reading {len(mfiles)} file(s)!')

            # for existing directory and files inside it, starts the reading of files     
            for k in range(len(mfiles)):

                if detect_netcdf(mfiles[k]) == None:
                    raise FileReaderError(f"--The following file is not in netcdf format: {mfiles[k]}")
                
                raw_data = xr.open_dataset(mfiles[k])
                
                if "Measurement_ID" not in raw_data.attrs:
                    raise FileReaderError(f"Measurement_ID parameter not found in the netcdf file. This is not a scc raw file")
                    
                # Reading the scc file metadata
                time_info_f = get_time_info(raw_data, 
                                            meas_type = meas_type,
                                            filename = mfiles[k].name)
                if time_info_f.empty:
                    return(system_info, channel_info, time_info, sig_raw, shots)
                
                if k == 0:
                    system_info = read_meas(raw_data = raw_data)            
                    channel_info = read_channels(raw_data = raw_data)
                else:
                    continue
                       
                # Reading the licel signals
                sig_raw_f = read_signals(raw_data, 
                                         time = time_info_f.index, 
                                         channels = channel_info.index, 
                                         meas_type = meas_type)
                
                shots_f = read_shots(raw_data, 
                                     time = time_info_f.index, 
                                     channels = channel_info.index,
                                     meas_type = meas_type)
                
                channel_info["bins"] = channel_info.index.size * [int(sig_raw_f.shape[-1] + 1)]

                # Append the arrays to list in order to concatenate later
                list_sig.append(sig_raw_f)
                list_shots.append(shots_f)                
                list_time.append(time_info_f)
                
            if len(list_sig) > 0:
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
            CustomWarning(f"No files to read in: {dir_meas}")  

    return(system_info, channel_info, time_info, sig_raw, shots)

def read_meas(raw_data):
    
    system_info = pd.Series()
        
    return(system_info)

def read_channels(raw_data):
    
    ch_index = raw_data["channel_ID"].values.astype(str)

    channel_info = pd.DataFrame(index = ch_index)
    channel_info["data_acquisition_range"] = raw_data["DAQ_Range"].values
    
    return(channel_info)

def get_time_info(raw_data, meas_type, filename):

    time_info = pd.Series()
    
    if meas_type == 'drk':
        if "Raw_Bck_Start_Time" not in list(raw_data.variables):
            print("--Raw_Bck_Start_Time parameter not found. No dark profile embedded in SCC file -> skipping")
            return(time_info)
        if "Raw_Bck_Stop_Time" not in list(raw_data.variables):
            print("--Raw_Bck_Stop_Time parameter not found. No dark profile embedded in SCC file -> skipping")
            return(time_info)
        sdate = raw_data.RawBck_Start_Date
        stime = raw_data.RawBck_Start_Time_UT
        start_time_sec = raw_data.Raw_Bck_Start_Time[:,0].values.astype(float)
        stop_time_sec = raw_data.Raw_Bck_Stop_Time[:,0].values.astype(float)
        filenames = np.empty(raw_data.time_bck.size, dtype = object)

    else:
        if "Raw_Data_Start_Time" not in list(raw_data.variables):
            print("--Raw_Data_Start_Time parameter not found. No profile other than dark embedded in SCC file -> skipping")
            return(time_info)
        if "Raw_Data_Stop_Time" not in list(raw_data.variables):
            print("--Raw_Data_Stop_Time parameter not found. No profile other than dark embedded in SCC file -> skipping")
            return(time_info)
        sdate = raw_data.RawData_Start_Date
        stime = raw_data.RawData_Start_Time_UT
        start_time_sec = raw_data.Raw_Data_Start_Time[:,0].values.astype(float)
        stop_time_sec = raw_data.Raw_Data_Stop_Time[:,0].values.astype(float)
        filenames = np.empty(raw_data.time.size, dtype = object)

    # Convert and store start time
    sdt = dt.datetime.strptime(sdate + ' ' + stime, "%Y%m%d %H%M%S") # start meas
                    
    start_time_arr = np.array([sdt + dt.timedelta(seconds = t) for t in start_time_sec])
    end_time_arr = np.array([sdt + dt.timedelta(seconds = t) for t in stop_time_sec])
    
    filenames[:] = filename
    
    tdata = np.array([filenames, 
                      datetimes_to_iso(start_time_arr), 
                      datetimes_to_iso(end_time_arr)], 
                     dtype = object)
        
    time_info = pd.DataFrame(tdata.T,  
                             index = start_time_arr,
                             columns = ['filename', 'start_time', 'end_time'])  

    return(time_info)

def read_signals(raw_data, time, channels, meas_type):
    
    if meas_type == 'drk':
        if "Background_Profile" not in list(raw_data.variables):
            raise FileReaderError("--Background_Profile parameter not found. Is this really a dark measurement file?")
            
        sig_arr = raw_data["Background_Profile"].values
    else:
        if "Raw_Lidar_Data" not in list(raw_data.variables):
            raise FileReaderError("--Raw_Lidar_Data parameter not found. Is this really a non-dark measurement file?")
        
        sig_arr = raw_data["Raw_Lidar_Data"].values
    
    sig_arr[sig_arr >= 9.96e+36] = np.nan 
    
    bins = 1. + np.arange(0, sig_arr.shape[-1])

    sig_raw = xr.DataArray(sig_arr, 
                           coords=[time, channels, bins], #range_sig
                           dims=['time', 'channel', 'bins']) #'range' 

    # Sort by time
    sig_raw = sig_raw.copy().sortby('time')
    
    
    return(sig_raw)

def read_shots(raw_data, time, channels, meas_type):
    
    if meas_type == "drk":
        shots = np.tile(np.median(raw_data["Laser_Shots"].values, axis = 0), (len(time), 1))
    else:
        shots = raw_data["Laser_Shots"].values
        

    shots = xr.DataArray(shots, 
                         coords=[time, channels], #range_sig
                         dims=['time', 'channel']) #'range' 

    # Sort by time
    shots = shots.copy().sortby('time')
    
    return(shots)
