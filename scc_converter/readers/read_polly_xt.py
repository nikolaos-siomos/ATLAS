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
        
        # for existing directory and files inside it, starts the reading of files     
        if len(mfiles) > 0:
            print(f'-- Folder contains {len(mfiles)} file(s)!')
                
            raw_data = xr.open_dataset(mfiles[0])
            
            # Reading the polly_xt file metadatas (header) - only for the first file            
            system_info = read_meas(raw_data = raw_data)            
            channel_info = read_channels(raw_data = raw_data)
            
            channels = channel_info.index.values
            # bins_arr = np.arange(1., channel_info.bins.max() + 1.)  
            bins_arr = np.arange(0., channel_info.bins.max())  
            
            for k in range(len(mfiles)):
            
                raw_data = xr.open_dataset(mfiles[k])

                # Mask measurements (normal, +45, -45) according to cal_angle (True:norm or +45 or -45 / False:other)     
                mask_zer, mask_p45, mask_m45 = \
                    get_cal_info(pol_cal_angle = raw_data.depol_cal_angle.values, 
                                 meas_type = meas_type)    
                
                filename = np.empty(raw_data.time.size, dtype = object)
                folder = np.empty(raw_data.time.size, dtype = object)
                position = np.empty(raw_data.time.size, dtype = object)
                
                filename[:] = os.path.basename(mfiles[k])
                
                if meas_type == 'pcb':
                    mask = (mask_p45) | (mask_m45)
                    folder[mask_p45] = '+45'
                    folder[mask_m45] = '-45'
                    position[mask_p45] = 2
                    position[mask_m45] = 1                
                elif meas_type == 'nrm':
                    mask = mask_zer                    
                    position[mask_zer] = 0
                    folder[mask_zer] = 'nrm'
                else:
                    mask = mask_zer                    
                    position[mask_zer] = 0   
                    
                if meas_type == 'tlc':
                    folder[:] = (mfiles[k]).split(os.sep)[-2] 
                
                if mask.any():
                    raw_signal = raw_data.raw_signal[mask,:,:].values.astype(float)
                    raw_shots = raw_data.measurement_shots.values[mask,:]
                    
                    position = position[mask]
                    filename = filename[mask]
                    folder = folder[mask]
                    
                    # Convert the time to npdatetime format and mask them
                    start_time = raw_data.measurement_time
    
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
            print('---- Warning! Folder empty \n'+\
                  f'---> !! Skip reading measurement files from folder {dir_meas}')  

    return(system_info, channel_info, time_info, sig_raw, shots)

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
    the polly_xt metadata [altitude, latitude, longitude, 
    zenith angle, azimuth angle] and laser relevant information from 
    the licel header [laser A repetion rate, laser B repetion rate if it exists
    laser C repetion rate if it exists]"""
    
    system_info = pd.Series()
    
    system_info['altitude'] = raw_data.location_height.values   
    system_info['latitude'] = raw_data.location_coordinates.values[0]
    system_info['longitude'] = raw_data.location_coordinates.values[1]
    
    system_info['latitude'] = np.round(system_info['latitude'], decimals = 4)
    system_info['longitude'] = np.round(system_info['longitude'], decimals = 4)
    
    system_info['zenith_angle'] = raw_data.zenithangle.values
    system_info['zenith_angle'] = np.round(system_info['zenith_angle'], decimals = 2)
    
    system_info['azimuth_angle'] = 0.
    
    system_info['laser_A_repetition_rate'] = raw_data.laser_rep_rate.values
    
    system_info['laser_B_repetition_rate'] = np.nan
        
    system_info['laser_C_repetition_rate'] = np.nan

    return(system_info)

def read_channels(raw_data):
    channel_info = pd.DataFrame(index = (raw_data.channel.values + 1).astype(str))
    
    channel_info.loc[:,'acquisition_mode'] = 1 # only photon channels
    channel_info.loc[:,'laser'] = 1 # only 1 laser
    channel_info.loc[:,'bins'] = raw_data.height[-1].values + 1
    channel_info.loc[:,'range_resolution'] = raw_data.measurement_height_resolution.values * 0.15 # nanosecond to meters
    channel_info.loc[:,'data_acquisition_range'] = np.nan # analog channels only
    channel_info.loc[:,'analog_to_digital_resolution'] = np.nan # analog channels only
    channel_info.loc[:,'detected_wavelength'] = raw_data.if_center.values
    channel_info.loc[:,'channel_bandwidth'] = raw_data.if_fwhm.values
    channel_info.loc[:,'dead_time_correction_type'] = 1. # default for PollyXTs

    
    return(channel_info)


def get_cal_info(pol_cal_angle, meas_type):
    
    if meas_type != 'ray' and meas_type != 'pcb':
        
        mask_zer = np.ones(pol_cal_angle.shape, dtype = bool)
        
        mask_p45 = np.zeros(pol_cal_angle.shape, dtype = bool)
        
        mask_m45 = np.zeros(pol_cal_angle.shape, dtype = bool)
  
    else:
        
        check_cal_angle(pol_cal_angle)
        
        fix_pos = find_fixed_pos(pol_cal_angle)
        
        zer_pos, p45_pos, m45_pos = assign_cal_positions(fix_pos)
        
        mask_zer, mask_p45, mask_m45 = \
            get_masks(pol_cal_angle, 
                      zer_pos = zer_pos, 
                      p45_pos = p45_pos, 
                      m45_pos = m45_pos)
        
    return(mask_zer, mask_p45, mask_m45)
            

def check_cal_angle(pol_cal_angle):
    
    mask_fin = np.isfinite(pol_cal_angle)
    
    if not mask_fin.all():
        print('-- Warning: The pol. calibrator angle provided in the raw files contains non-finite values. This can lead to problems with the separation of the pol. cal. measurement from the normal measurment. Please make sure that the separation works as expected')
    
    mask_err = (np.abs(pol_cal_angle) <= 360) | (pol_cal_angle == 999.) | (~np.isfinite(pol_cal_angle))
    
    if not mask_err.all():
        raise Exception(f'-- Error: The absolute pol. calibrator angle provided in the raw files above 360 degrees without being 999. For example {pol_cal_angle[~mask_err]}. The value 999 is allowed in some old polly_xt systems showing that the calibrator is at the zero position')
    
    return()

def find_fixed_pos(pol_cal_angle):
    
    mask_dif = np.ones(pol_cal_angle.shape, dtype = bool)

    mask_dif[1:] = (np.abs(pol_cal_angle[1:] - pol_cal_angle[:-1]) <= 0.5) & \
        (np.isfinite(pol_cal_angle[1:])) & (np.isfinite(pol_cal_angle[:-1]))
    
    mask_dif[1:-1] = (mask_dif[1:-1]) & (mask_dif[2:])

    uniques, counts =np.unique(pol_cal_angle[mask_dif], return_counts=True)
    
    sorted_uniques = uniques[np.argsort(counts)]
    
    if sorted_uniques.size == 1:
        print('-- Warning: The provided normal measurement does not seem to include a calibration measurement. Please make sure that this is really the case')
    
    elif sorted_uniques.size == 2:
        raise Exception('-- Error: The provided normal measurement does not seem to include a normal part. Please correct the input file or report this issue to CARS if not necessary')
    
    fix_pos = uniques[np.argsort(counts)][:3]
    
    return(fix_pos)

def assign_cal_positions(fix_pos):
    
    mask_zer_pos = (np.abs(fix_pos) < 40.) | (np.abs(fix_pos -360.) < 40.) | (fix_pos == 999.)
    
    if np.sum(mask_zer_pos) == 0:
        raise Exception(f'-- Error: It is not possible to assume the zero calibrator position because the absolute calibrator angle in any of the 3 fixed positions {fix_pos} is neither less than 40 degrees nor between 320 and 360 degrees nor is it exactly 999. Please correct the input file or report this issue to CARS if not necessary ')
    elif np.sum(mask_zer_pos) > 1:
        raise Exception(f'-- Error: It is not possible to assume the zero calibrator position because the absolute calibrator angle is less than 40 degrees or between 320 and 360 degrees or exactly 999 in more than one fixed positions {fix_pos}. Please correct the input file or report this issue to CARS if not necessary ')
    else:
        zer_pos = fix_pos[mask_zer_pos][0]

    if len(fix_pos) == 3:
        mask_p45_pos = (fix_pos != zer_pos) & (((fix_pos > 0.) & (fix_pos < 90.)) | ((fix_pos > 180.) & (fix_pos < 270.)))
                                               
        if np.sum(mask_p45_pos) == 0:
            raise Exception(f'-- Error: It is not possible to assume the +45 calibrator position because the absolute calibrator angle in any of the 2 fixed non zero positions {fix_pos} is not within 0 and 90 degrees nor within 180 and 270 degrees. Please correct the input file or report this issue to CARS if not necessary ')
        elif np.sum(mask_p45_pos) > 1:
            raise Exception(f'-- Error: It is not possible to assume the +45 calibrator position because the absolute calibrator angle is within 0 and 90 degrees or within 180 and 270 degrees in both fixed non zero positions {fix_pos}. Please correct the input file or report this issue to CARS if not necessary ')
        else:
            p45_pos = fix_pos[mask_p45_pos][0]
    else:
        p45_pos = np.nan
    
    if len(fix_pos) == 3:

        mask_m45_pos = (fix_pos != zer_pos) & (((fix_pos < 0.) & (fix_pos > -90.)) | ((fix_pos > 90.) & (fix_pos < 180.)) | ((fix_pos > 270.) & (fix_pos < 360.)))
        
        if np.sum(mask_m45_pos) == 0:
            raise Exception(f'-- Error: It is not possible to assume the -45 calibrator position because the absolute calibrator angle in any of the 2 fixed non zero positions {fix_pos} is not within 0 and -90 degrees nor within 90 and 180 degrees nor within 270 and 360 degrees. Please correct the input file or report this issue to CARS if not necessary ')
        elif np.sum(mask_p45_pos) > 1:
            raise Exception(f'-- Error: It is not possible to assume the -45 calibrator position because the absolute calibrator angle is within 0 and -90 degrees or within 90 and 180 degrees or within 270 and 360 degrees in both fixed non zero positions {fix_pos}. Please correct the input file or report this issue to CARS if not necessary ')
        else:
            m45_pos = fix_pos[mask_m45_pos][0]
    else:
        m45_pos = np.nan
        
    return(zer_pos, p45_pos, m45_pos)

def get_masks(pol_cal_angle, zer_pos, p45_pos, m45_pos):
    
    mask_zer = ((np.abs(pol_cal_angle - zer_pos)) <= 0.05) & \
        (np.isfinite(pol_cal_angle))
    
    mask_zer = trim_edges(mask_zer)
    
    if p45_pos == p45_pos:
        mask_p45 = (np.abs(pol_cal_angle - p45_pos) <= 0.05) & \
            (np.isfinite(pol_cal_angle))
        mask_p45 = trim_edges(mask_p45)

    else:
       mask_p45 = np.zeros(pol_cal_angle.size, dtype = bool)
       
    if m45_pos == m45_pos:
        mask_m45 = (np.abs(pol_cal_angle - m45_pos) <= 0.05) & \
            (np.isfinite(pol_cal_angle))
        mask_m45 = trim_edges(mask_m45)
    else:
        mask_m45 = np.zeros(pol_cal_angle.size, dtype = bool)
      
    return(mask_zer, mask_p45, mask_m45)

def trim_edges(mask):
    
    mask[1:-1] = (mask[2:].astype(int) - mask[1:-1].astype(int) == 0) & \
        (mask[1:-1].astype(int) - mask[:-2].astype(int) == 0) & \
            (mask[1:-1] == True)
    mask[0] = (mask[1].astype(int) - mask[0].astype(int) == 0) & (mask[0] == True)
    mask[-1] = (mask[-1].astype(int) - mask[-2].astype(int) == 0) & (mask[-1] == True)
    
    return(mask)
    
            