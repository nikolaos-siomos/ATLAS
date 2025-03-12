import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime as dt
from datetime import timedelta
import xarray as xr
from scipy.io import loadmat

# Read measurement
def dtfs(dir_meas):
    
    """ Reads information from the raw licel files"""
    
    # Setting sig, info, and time as empty lists in the beggining    
    sig_raw = []     
    shots = []
    folder = []
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
            print(f'-- Folder contains {len(mfiles)} file(s)!')
            
            data = loadmat(mfiles[0])
            keys = [key for key in data.keys()]
            
            # Reading the licel file metadatas (header) - only for the first file
            system_info = read_meas(buffer = data[keys[5]])

            channel_info = read_channels(buffer = data[keys[4]])
            
            channels = channel_info.index.values
            # bins_arr = np.arange(1., channel_info.bins.max() + 1.)
            bins_arr = np.arange(0., channel_info.bins.max())

            # Creating empty signal, shots, and time arrays
            start_time_arr = np.nan*np.zeros(len(mfiles), dtype = object)
            end_time_arr = np.nan*np.zeros(len(mfiles), dtype = object)

            shots_arr = np.nan*np.zeros((len(mfiles), len(channels)), dtype = object)
            sig_arr = np.nan*np.zeros((len(mfiles), len(channels), len(bins_arr)), dtype = float)

            filename = np.empty(len(mfiles), dtype = object)
            folder = np.empty(len(mfiles), dtype = object)

            # Iterate over the files
            for k in range(len(mfiles)):
                
                filename[k] = os.path.basename(mfiles[k])
                
                data = loadmat(mfiles[k])
                
                keys = [key for key in data.keys()]
                
                stime, etime = read_time(buffer = data[keys[5]])

                shots_arr[k,:] = read_shots(buffer = data[keys[5]])
                
                body = np.array([data[keys[6+i]] for i in range(len(data)-6)])

                # Store signal, start and end time
                sig_arr[k, :, :] = body[:,0,:]

                start_time_arr[k] = stime
                
                if stime == etime: #only possible if the files are different by only milliseconds 
                    end_time_arr[k] = etime + timedelta(milliseconds = 500)
                    print(f'-- Warning! File {filename[k]} has the same start and end time reported (recording lasted < 1s). Please check it! ')
                else:
                    end_time_arr[k] = etime

                     
                if (mfiles[k]).split(os.sep)[-2] in ['north', 'east', 'south', 'west', 'inner', 'outer', '+45', '-45', 'static']:
                    folder[k] = (mfiles[k]).split(os.sep)[-2]
            
            sig_raw = xr.DataArray(sig_arr, 
                                   coords=[end_time_arr, channels, bins_arr],
                                   dims=['time', 'channel', 'bins']) 
            
            shots = xr.DataArray(shots_arr,  
                                 coords=[end_time_arr, channels],
                                 dims=['time', 'channel'])
            
            tdata = np.array([folder, filename, 
                              start_time_arr, end_time_arr], dtype = object)

            properties = ['folder', 'filename', 'start_time', 'end_time']
            
            time_info = pd.DataFrame(tdata.T,  
                                     index = end_time_arr,
                                     columns = properties)  
                        
            # Sort by time
            sig_raw = sig_raw.sortby('time').copy()
            shots = shots.sortby('time').copy()
            time_info = time_info.sort_index()
            
        else:
            print('---- Warning! Folder empty \n'+\
                  f'---> !! Skip reading measurement files from folder {dir_meas}')  

    return(system_info, channel_info, time_info, sig_raw, shots)


def read_meas(buffer):

    """ Retrieves location and geometry relevant information from 
    the licel header [altitude, latitude, longitude, 
    zenith angle, azimuth angle] and laser relevant information from 
    the licel header [laser A repetion rate, laser B repetion rate if it exists
    laser C repetion rate if it exists]"""
    system_info = pd.Series()
     
    metadata = buffer[2][0][0].split()

    system_info['altitude'] = float(metadata[5])    
    system_info['latitude'] = np.round(float(metadata[6]), 4)
    system_info['longitude'] = np.round(float(metadata[7]), 4)
    
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
    recorder_channel_id = []
    for i in range(len(header[:,-1])):
        recorder_channel_id.append(f'{arr_head.recorder_channel_id.iloc[i]}_L{str(int(arr_head.laser.iloc[i]))}')  
    # Check if the defined channels are unique (unique sets of licel id and laser number)   
        if len(set(recorder_channel_id)) < len(recorder_channel_id):
            raise Exception('-- Error: At least two of the licel channels have both the same id and laser number. Please correct this in the recorder settings')

    channel_info.index = recorder_channel_id

    info_columns = ['acquisition_mode', 'laser', 'bins', 
                    'range_resolution', 'data_acquisition_range', 
                    'analog_to_digital_resolution']
    
    channel_info.loc[:, 'recorder_channel_id'] = arr_head.loc[:, 'recorder_channel_id'].copy().values
    channel_info.loc[:, info_columns] = arr_head.loc[:, info_columns].copy().values.astype(float)

    mask_an = channel_info.loc[:,'acquisition_mode'].values == 0
    channel_info.loc[:,'data_acquisition_range'][mask_an] = 1000.*channel_info.loc[:,'data_acquisition_range'][mask_an]
    channel_info.loc[:,'data_acquisition_range'][~mask_an] = np.nan

    wave = np.array(list(np.char.split(arr_head.wave_pol.values.astype('str'),
                                       sep='.')))[:,0].astype(float)

    channel_info.loc[:,'detected_wavelength'] = wave
    
    channel_info.loc[:,'channel_bandwidth'] = 1. # default when not available in the raw files
    
    channel_info.loc[:,'dead_time_correction_type'] = 0. # default for Licel

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