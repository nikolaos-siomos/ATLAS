"""
@author: Peristera

Read the raw PollyXT data
"""
import glob
import pandas as pd
import numpy as np
import os
import xarray as xr


def dtfs(dir_meas, meas_type):
    
    # Setting sig, info, and time as empty lists in the beggining    
    sig_raw = []     
    shots = []
    time_info = []
    
    meas_info = []
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
        
        # for existing directory and files inside it, starts the reading of files     
        if len(mfiles) > 0:
            print(f'-- Folder contains {len(mfiles)} file(s)!')
                
            raw_data = xr.open_dataset(mfiles[0])
            
            # Reading the polly_xt file metadatas (header) - only for the first file            
            meas_info = read_meas(raw_data = raw_data)            
            channel_info = read_channels(raw_data = raw_data)
            
            channels = channel_info.index.values
            bins_arr = np.arange(1., channel_info.bins.max() + 1.)  
            
            for k in range(len(mfiles)):
            
                raw_data = xr.open_dataset(mfiles[k])

                # Mask measurements (normal, +45, -45) according to cal_angle (True:norm or +45 or -45 / False:other)     
                if meas_type == 'pcb':
                    mask, mask_p45, mask_m45 = get_cal_info(raw_data)    
                else:
                    mask = (raw_data.depol_cal_angle.values == 0.)

                    
                raw_signal = raw_data.raw_signal[mask,:,:].values.astype(float)
                start_time = raw_data.measurement_time
                raw_shots = raw_data.measurement_shots.values[mask,:]
                
                filename = np.empty(raw_data.time.size, dtype = object)
                folder = np.empty(raw_data.time.size, dtype = object)
                
                position = np.empty(raw_data.time.size, dtype = object)
                
                filename[:] = os.path.basename(mfiles[k])
                
                if meas_type == 'pcb':
                    folder[mask_p45] = '+45'
                    folder[mask_m45] = '-45'
                    position[mask_p45] = 2
                    position[mask_m45] = 1
                    
                elif meas_type == 'tlc':
                    folder[:] = (mfiles[k]).split(os.sep)[-2] 

                position = position[mask]
                filename = filename[mask]
                folder = folder[mask]
                
                # Convert the time to npdatetime format and mask them
                start_time_arr = convert_time_to_npdatetime(start_time) 

                start_time_arr = start_time_arr[mask]
                end_time_arr = start_time_arr + (start_time_arr[1] - start_time_arr[0])
            
                # Define signal xarray
                sig_f = xr.DataArray(raw_signal, 
                                   coords=[end_time_arr, bins_arr, channels],
                                   dims=['time', 'bins', 'channel'])
                                
                shots_f = xr.DataArray(raw_shots, 
                                       coords=[end_time_arr, channels],
                                       dims=['time', 'channel']) 

                # Define temporal data pandas Dataframe
                if meas_type == 'pcb':
                    properties = ['folder', 'filename', 
                                  'start_time', 'end_time', 'position']
                    tdata = np.array([folder, filename, 
                                      start_time_arr, end_time_arr, position], 
                                     dtype = object)
                    
                else:
                    properties = ['folder', 'filename', 'start_time', 'end_time']
                    tdata = np.array([folder, filename, 
                                      start_time_arr, end_time_arr], 
                                     dtype = object)
                    
                time_info_f = pd.DataFrame(tdata.T,  
                                           index = end_time_arr,
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
            print('---- Warning! Folder empty \n'+\
                  f'---> !! Skip reading measurement files from folder {dir_meas}')  


    return(meas_info, channel_info, time_info, sig_raw, shots)

def convert_time_to_npdatetime(time):
    # meas_time: 2d array : [times, 2], 2: [yyyymmdd, second of day]
    frames = time.shape[0]
    time_dt = []
    
    for i in range(frames): 
        temp_str = str(time.values[i][0])
        hr = np.divide(time.values[i][1],(60*60)); minute = (hr - np.fix(hr))*60; sec = (minute - np.fix(minute))*60
        date = temp_str[0:4]+'-'+temp_str[4:6]+'-'+temp_str[6:8]+'T'+"%02d" %np.fix(hr)+':'+"%02d" %np.fix(minute)+':'+"%02d" %np.fix(sec)
        time_dt.append(np.datetime64(date))
    
    return(np.array(time_dt))

def read_meas(raw_data):

    """ Retrieves location and geometry relevant information from 
    the polly_xt metadata [location, altitude, latitude, longitude, 
    zenith angle, azimuth angle] and laser relevant information from 
    the licel header [laser A repetion rate, laser B repetion rate if it exists
    laser C repetion rate if it exists]"""
    
    meas_info = pd.Series()

    meas_info['location'] = raw_data.location
    
    meas_info['altitude'] = raw_data.location_height.values   
    meas_info['latitude'] = np.round(raw_data.location_coordinates.values[0], 4)
    meas_info['longitude'] = np.round(raw_data.location_coordinates.values[1], 4)
    
    meas_info['zenith_angle'] = raw_data.zenithangle.values
    meas_info['azimuth_angle'] = 0.
    
    meas_info['laser_A_repetition_rate'] = raw_data.laser_rep_rate.values

    return(meas_info)

def read_channels(raw_data):
    channel_info = pd.DataFrame(index = raw_data.channel)
    
    channel_info.loc[:,'acquisition_mode'] = 1 # only photon channels
    channel_info.loc[:,'laser'] = 1 # only 1 laser
    channel_info.loc[:,'bins'] = raw_data.height[-1] + 1
    channel_info.loc[:,'laser_polarization'] = 1 # only linear polarization
    channel_info.loc[:,'pmt_high_voltage'] = raw_data.pm_voltage.values
    channel_info.loc[:,'range_resolution'] = raw_data.measurement_height_resolution.values * 0.15 # nanosecond to meters
    channel_info.loc[:,'data_acquisition_range'] = np.nan # analog channels only
    channel_info.loc[:,'analog_to_digital_resolution'] = np.nan # analog channels only
    channel_info.loc[:,'detected_wavelength'] = raw_data.if_center.values
    channel_info.loc[:,'channel_bandwidth'] = raw_data.if_fwhm.values
    
    # Fill in default values that depend on the raw file metadata
    channel_info.loc[:,'background_low_bin'] = channel_info.loc[:,'bins'] - 600
    channel_info.loc[:,'background_high_bin'] = channel_info.loc[:,'bins'] - 100
    channel_info.loc[:,'emitted_wavelength'] = np.nan * channel_info.loc[:,'detected_wavelength'] 
    
    for i in range(channel_info.loc[:,'detected_wavelength'] .size):
        if channel_info.loc[i,'detected_wavelength'] >= 340. and \
            channel_info.loc[i,'detected_wavelength'] < 520.:
                channel_info.loc[i,'emitted_wavelength'] = 355.
                
        if channel_info.loc[i,'detected_wavelength'] >= 520. and \
            channel_info.loc[i,'detected_wavelength'] < 1000.:
                channel_info.loc[i,'emitted_wavelength'] = 532.
                
        if channel_info.loc[i,'detected_wavelength'] >= 1000.:
            channel_info.loc[i,'emitted_wavelength'] = 1064.
                
    return(channel_info)


def get_cal_info(raw_data):
        
    pre_mask = (raw_data.depol_cal_angle.values != 0.)

    angles = raw_data.depol_cal_angle.values
    
    cal_angle = np.bincount(np.ceil(angles[pre_mask]).astype(int)).argmax()
    
    perfect_angles = np.array([45, 225, 135, 315])
    idx = (np.abs(perfect_angles - cal_angle)).argmin()

    if perfect_angles[idx] in [45, 225]:
        mask_p45 = (np.round(angles - cal_angle, decimals = 0) == 0)
        mask_m45 = ((np.round(angles - cal_angle - 90., decimals = 0) == 0.) |\
                    (np.round(angles - cal_angle + 90., decimals = 0) == 0.))

    if perfect_angles[idx] in [135, 315]:
        mask_m45 = (np.round(angles - cal_angle, decimals = 0) == 0)
        mask_p45 =((np.round(angles - cal_angle - 90., decimals = 0) == 0.) |\
                    (np.round(angles - cal_angle + 90., decimals = 0) == 0.))        

    mask = mask_p45 | mask_m45
        
    return(mask, mask_p45, mask_m45)