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
import os
import numpy as np
import pandas as pd
import glob
import datetime as dt
import xarray as xr
from lidar_processing import signal
import timeit, sys

# Read measurement
def dtfs(dir_meas, cfg, meas_type):
    
    # Setting sig, info, and time as empty lists in the beggining    
    sig_raw= []    
    info = []
    ground_alt, SZA, azimuth = [], [], []
    
    if not(os.path.exists(dir_meas)):
        sys.exit('---- Error : The folder for reading signals does not exist! \r\n' +\
                 f'---- Check the input directory! \n Given folder: {dir_meas}')
    
    else:
        
        mfile = glob.glob(os.path.join(dir_meas,'*.nc'))
    
    
        # for existing directory and files inside it, starts the reading of files     
        if len(mfile) == 1:
            print('-- SCC input file found!')
            
            file = xr.open_dataset(mfile[0])
            
            # Read header once
            info = read_scc_info(file, cfg)
            
            # Reading the scc file metadatas
            sdt, edt, time, ground_alt, SZA, azimuth = read_metas(file, cfg = cfg, meas_type = meas_type)
                   
            # Reading the licel signals
            sig_raw = read_scc_signals(file, info = info, time = time, channels = info.index, meas_type = meas_type)
            
        if len(mfile) == 0:
            print('---- Warning! Folder empty \n'+\
                  f'---> !! Skip reading measurement files from folder {dir_meas}')  

        if len(mfile) > 1:
            sys.exit('---- Error! More than one SCC input files found in folder \n'+\
                  f'{dir_meas} ---> !! Aborting ')  


    return(sig_raw, [], info, [], ground_alt, SZA, azimuth)

def read_scc_info(file, cfg):

    # Header rows
    scc_channels = file.channel_ID.values
  
    cfg_channels = cfg.ch_list
    
    # Link channel ID from config file to the SCC IDs, row-wise 
    channels = scc_to_cfg(scc_channels, cfg_channels)

    info = pd.DataFrame(index = channels, dtype = object)
    
    info.loc[:,'ch_mode'] = cfg.ch_mode.values
    info.loc[:,'bins'] = cfg.bins.values
    info.loc[:,'laser_pol'] = cfg.laser_pol.values
    info.loc[:,'resol'] = cfg.raw_resol
    info.loc[:,'shots'] = np.mean(file.Laser_Shots,axis = 0)
    info.loc[:,'ADC_range'] = cfg.ADC_range
    info.loc[:,'ADC_bit'] = np.nan
    info.loc[:,'wave'] = cfg.wave
    info.loc[:,'ch_pol'] = cfg.ch_pol
    info.loc[:, 'sampl_rate'] =  150. / info.resol.values  # MHz
    
    return(info)

def read_metas(file, cfg, meas_type):

    if meas_type == 'normal':
        sdate = file.RawData_Start_Date
        stime = file.RawData_Start_Time_UT
        etime = file.RawData_Stop_Time_UT
        tscale = file.Raw_Data_Start_Time[:,0].values.astype(float)
    elif meas_type == 'dark':
        sdate = file.RawBck_Start_Date
        stime = file.RawBck_Start_Time_UT
        etime = file.RawBck_Stop_Time_UT
        tscale = file.Raw_Bck_Start_Time[:,0].values.astype(float)
    else:
        sys.exit('---- Error: Measurement type other than normal, or dark. '+\
                 '---- Please correct and reprocess ---> Aborting !!')

    # Convert and store start time
    sdt = dt.datetime.strptime(sdate + ' ' + stime, "%Y%m%d %H%M%S") # start meas
    edt = dt.datetime.strptime(sdate + ' ' + etime, "%Y%m%d %H%M%S") # start meas
                    
    time = np.array([sdt + dt.timedelta(seconds = t) for t in tscale])
    

    ground_alt = cfg.lidar['altitude']
    SZA = file.Laser_Pointing_Angle.values[0]
    azimuth = cfg.angles['azimuth']
    
    return(sdt, edt, time, ground_alt, SZA, azimuth)

def read_scc_signals(file, info, time, channels, meas_type):
    
    sig_arr = np.nan*np.zeros(file.Raw_Lidar_Data.shape, dtype = float)

    if meas_type == 'normal' or meas_type == 'calibration':
        sig_arr = file.Raw_Lidar_Data.values
    elif meas_type == 'dark':
        sig_arr = file.Background_Profile.values
    else:
        sys.exit('---- Error: Measurement type other than normal, dark, or calibration. '+\
                 '---- Please correct and reprocess ---> Aborting !!')
                 
    bins = 1. + np.arange(0, sig_arr.shape[-1])

    sig_raw = xr.DataArray(sig_arr, 
                           coords=[time, channels, bins], #range_sig
                           dims=['time', 'channel', 'bins']) #'range' 

    # Sort by time
    sig_raw = sig_raw.copy().sortby('time')
    
    # Photon units conversion from summed counts to MHz
    sig_raw = signal.unit_conv_counts_to_MHz(sig_raw.copy(), info)

    return(sig_raw)


def scc_to_cfg(scc_channels, cfg_channels):
        
    if len(cfg_channels)>0:
        
        if len(cfg_channels) != len(scc_channels):
            raise ValueError('config_file channels ID and LICEL channels ID '+\
                             'have not the same length! Please revise the config_file.ini!')
        
        #Replace the LICEL channels ID with the IDs given by the user in config_file
        scc_channels = cfg_channels
        
    return(scc_channels)
