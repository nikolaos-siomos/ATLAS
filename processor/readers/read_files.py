"""
@authors: Nikolaos Siomos (nikolaos.siomos@lmu.de)

================
"""
import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import netCDF4 as nc

# Read measurement
def short_reader(fpath, exclude_telescope_type, exclude_channel_type, 
                 exclude_acquisition_mode, exclude_channel_subtype, 
                 use_channels):
    
    """
    General:
        Collects all the necessary information from the QA file
        
    Input:
        fpath: 
            The path to the input QA netcdf file (string)
        
        exclude_telescope_type:
            Channel with a field type (1rd letter of the ID) in
            exclude_telescope_type will be excluded 
                    
        exclude_channel_type:
            Channel with a channel type (2rd letter of the ID) in
            exclude_channel will be excluded 
                    
        exclude_acquisition_mode:
            Channel with a channel type (3rd letter of the ID) in
            exclude_mode will be excluded 
            
        exclude_channel_subtype:
            Channel with a subtype (4rd letter of the ID) in
            exclude_subtype will be excluded 
            
        use_channels:
            Channel IDs to include for the calculations
            
            
    Returns:
        
        system_info: 
            A pandas series with all the necessary system metadata fetched 
            from the QA file
            
        channel_info : 
            A pandas Dataframe with all the necessary metadata for each 
            channel fetched from the QA file
            
        time_info : 
            A pandas Dataframe with all the necessary metadata for each 
            temporal frame fetched from the QA file
            
        time_info_d : 
            A pandas Dataframe with all the necessary metadata for each
            temporal frame of the dark measurements fetched from the QA file
            Defaults to an empty list if a dark measurement is not provided
    
        signal : 
            A 3D xarray Dataarray with the measured signals, extracted 
            from the QA file. The dimensions should be (time, channel, bins/range)
    
        signal_d : 
            A 3D xarray Dataarray with the measured dark signals, extracted 
            from the QA file. The dimensions should be (time, channel, bins/range)
            Defaults to an empty list if a dark measurement is not provided
    
        shots : 
            A 2D xarray Dataarray with the shots per timeframe and channel, 
            extracted from the QA file. The dimensions should be (time, channel)
    
        shots_d : 
            A 2D xarray Dataarray with the dark shots per timeframe and channel, 
            extracted from the QA file. The dimensions should be (time, channel)
            Defaults to an empty list if a dark measurement is not provided
            
    """


    # Setting signal, info, and time as empty lists in the beggining    
    signal = []
    signal_d = []
    
    shots = []
    shots_d = []   
    
    time_info = [] 
    time_info_d = []
        
    # Read QA netcdf file
    file = xr.open_dataset(fpath)
        
    # Read channel metadata
    channel_info = \
        channel_metadata(file)
    
    # Mask out channels based on their characteristics
    channels = channel_info.index.values
    
    if use_channels == None:
        use_channels = channels
    
    mask = np.array([ch[4] not in exclude_telescope_type and
                     ch[5] not  in exclude_channel_type and
                     ch[6] not  in exclude_acquisition_mode and
                     ch[7] not  in exclude_channel_subtype and
                     ch in use_channels for ch in channels])
    
    if all(~mask):
        raise Exception('-- Error: The provided channel filtering arguments are too strict and exclude all channels. Please revide the following arguments: exclude_telescope_type, exclude_channel_type, exclude_acquisition_mode, exclude_channel_subtype, channels')
    
    valid_channels = channels[mask]

    # Read lidar metadata
    system_info = lidar_metadata(file)
    
    # Read time 
    time_info = time_metadata(file, system_info)
    if 'Background_Profile' in file.data_vars: 
        time_info_d = time_metadata_d(file, system_info)

    # Reading the licel signals
    signal = signals(file, time_info = time_info, channel_info = channel_info)
    signal = signal.copy().loc[dict(channel = valid_channels)]
    if 'Background_Profile' in file.data_vars: 
        signal_d = signals(file, time_info = time_info_d, 
                           channel_info = channel_info, isdark = True)
        signal_d = signal_d.copy().loc[dict(channel = valid_channels)]


    # Reading the laser shots
    shots = laser_shots(file, time_info = time_info, channel_info = channel_info)
    shots = shots.copy().loc[dict(channel = valid_channels)]
    if 'Background_Profile' in file.data_vars: 
        shots_d = laser_shots(file, time_info = time_info_d, 
                              channel_info = channel_info, isdark = True)
        shots_d = shots_d.copy().loc[dict(channel = valid_channels)]

    channel_info = channel_info.loc[valid_channels,:]

    return(system_info, channel_info, time_info, time_info_d,
           signal, signal_d, shots, shots_d)

def lidar_metadata(file):
    
    """
    General:
        Extracts the general lidar metadata from the QA file
        
    Input:
        fpath: 
            The path to the input QA netcdf file (string)
            
    Returns:
        
        system_info: 
            A pandas series with all the necessary lidar metadata fetched 
            from the QA file
          
    """

    keys_attrs = ['Sounding_File_Name',
                  'Rayleigh_File_Name',
                  'Lidar_Name',
                  'Lidar_ID',
                  'Station_Name',
                  'Station_ID',
                  'Version_Name',
                  'Version_ID',
                  'Configuration_Name',
                  'Configuration_ID',
                  'Altitude_meter_asl',
                  'Latitude_degrees_north',
                  'Longitude_degrees_east',
                  'Measurement_ID',
                  'Measurement_type',
                  'RawBck_Start_Date',
                  'RawBck_Start_Time_UT',
                  'RawBck_Stop_Time_UT',
                  'RawData_Start_Date',
                  'RawData_Start_Time_UT',
                  'RawData_Stop_Time_UT']
    
    data_attrs = [file.attrs[key] for key in keys_attrs if key in file.attrs]
    keys_attrs = [key for key in keys_attrs if key in file.attrs]
    
    keys_vars = ['Molecular_Calc',
                 'Pressure_at_Lidar_Station',
                 'Temperature_at_Lidar_Station']

    data_vars = [file.variables[key].values for key in keys_vars 
                 if key in file.variables]
    keys_vars = [key for key in keys_vars if key in file.variables]
    
    keys_extras = ['Laser_Pointing_Angle',
                   'Laser_Pointing_Azimuth_Angle']

    data_extras = [file.variables[key].values[0] for key in keys_extras 
                   if key in file.variables]

    keys_extras = [key for key in keys_extras if key in file.variables]

    data = []
    data.extend(data_attrs)
    data.extend(data_vars)
    data.extend(data_extras)
    
    keys = []
    keys.extend(keys_attrs)
    keys.extend(keys_vars)
    keys.extend(keys_extras)
    
    system_info = pd.Series(data = data, index = keys, dtype = object)
    
    return(system_info)

def channel_metadata(file):
    
    """
    General:
        Extracts the metadata for each channel from the QA file
        
    Input:
        fpath: 
            The path to the input QA netcdf file (string)
            
    Returns:
            
        channel_info : 
            A pandas Dataframe with all the necessary metadata for each 
            channel fetched from the QA file

    """

    channels = np.array([np.round(file.Detected_Wavelength.values,decimals=0)\
                         .astype('int').astype(str)[i].zfill(4) +\
                         file.atlas_channel_label.values.astype(str)[i] \
                         for i in range(file.atlas_channel_label.size)])

    keys = ['ADC_resolution',
            'Background_Low_Bin',
            'Background_High_Bin',
            'Channel_Bandwidth',
            'channel_ID',
            'DAQ_Range',
            'Dead_Time',
            'Dead_Time_Correction_Type',
            'Detected_Wavelength',
            'Emitted_Wavelength',
            'Laser_Repetition_Rate',
            'Raw_Data_Range_Resolution',
            'DAQ_Trigger_Offset']
    
    data = [file.variables[key].values for key in keys if key in file.variables]
    keys = [key for key in keys if key in file.variables]

    channel_info = pd.DataFrame(data = np.array(data, dtype = object).T, 
                                index = channels, 
                                columns = keys, dtype = object)
    
    if (channel_info.values == nc.default_fillvals['f8']).any():
        channel_info[channel_info.values == nc.default_fillvals['f8']] = np.nan
    if (channel_info.values == nc.default_fillvals['i4']).any():
        channel_info[channel_info.values == nc.default_fillvals['i4']] = np.nan

    return(channel_info)

def time_metadata(file, system_info):
    
    """
    General:
        Collects all the metadata for each timeframe from the QA file
        
    Input:
        file: 
            An xarray dataset with all the QA file information

        system_info: 
            A pandas series with all the necessary lidar metadata fetched 
            from the QA file
            
    Returns:
        
        time_info : 
            A pandas Dataframe with all the necessary metadata for each 
            temporal frame fetched from the QA file
            
    """
    
    all_keys = ['Raw_Data_Start_Time',
                'Raw_Data_Stop_Time',
                'Filename',
                'sector',
                'calibrator_position']
    
    types = [int, int, str, str, str]
    
    data = []
    keys = []
    
    for i in range(len(all_keys)):
        
        if all_keys[i] in file.variables:

            var = file.variables[all_keys[i]].values.astype(types[i])
            
            if len(var.shape) == 2:
                data_i = np.squeeze(var, axis = 1)
    
            else:
                data_i = var

            data.append(data_i)
            keys.append(all_keys[i])

    data = np.array(data, dtype = object)    
    keys = np.array(keys, dtype = object) 
    
    sdate = system_info['RawData_Start_Date']

    stime = system_info['RawData_Start_Time_UT']

    start_t = file['Raw_Data_Start_Time'].values[:,0].astype(float)

    timeframes = time_array(sdate = sdate, stime = stime, start_t = start_t)
    
    time_info = pd.DataFrame(data = np.array(data, dtype = object).T, 
                             index = timeframes, 
                             columns = keys, dtype = object)
    
    return(time_info)

def time_metadata_d(file, system_info):
    
    """
    General:
        Collects all the metadata for each dark measurment 
        timeframe from the QA file
        
    Input:
        file: 
            An xarray dataset with all the QA file information

        system_info: 
            A pandas series with all the necessary lidar metadata fetched 
            from the QA file
            
    Returns:
        
        time_info_d : 
            A pandas Dataframe with all the necessary metadata for each 
            dark measurement timeframe fetched from the QA file
            
    """
    
    all_keys = ['Bck_Data_Start_Time',
                'Bck_Data_Stop_Time',
                'Filename_Bck']

    types = [int, int, str]
    
    data = []
    keys = []
    
    for i in range(len(all_keys)):
        
        if all_keys[i] in file.variables:

            var = file.variables[all_keys[i]].values.astype(types[i])
            
            if len(var.shape) == 2:
                data_i = np.squeeze(var, axis = 1)
    
            else:
                data_i = var

            data.append(data_i)
            keys.append(all_keys[i])

    data = np.array(data, dtype = object)    
    keys = np.array(keys, dtype = object)  

    sdate = system_info['RawBck_Start_Date']

    stime = system_info['RawBck_Start_Time_UT']

    start_t = file['Bck_Data_Start_Time'].values[:,0].astype(float)

    timeframes = time_array(sdate = sdate, stime = stime, start_t = start_t)

    time_info = pd.DataFrame(data = np.array(data, dtype = object).T, 
                             index = timeframes, 
                             columns = keys, dtype = object)
    
    return(time_info)

def time_array(sdate, stime, start_t):
    
    """
    General:
        Calculates the times for each timeframe in a datetime64[ns] format
        
    Input:
        sdate: 
            A string with the measurement starting date information (yyyymmdd)

        stime: 
            A string with the measurement starting time information (hhmmss)

        start_t:
            An 1D numpy array with the starting time of each measurement 
            time frame - in seconds after the measurement start  
            
    Returns:
        
        timeframes : 
            A 1D numpy array with the datetime objects for each measurement
            starting time
            
    """
    
    # Convert and store start time
    sdt = dt.datetime.strptime(sdate + ' ' + stime, "%Y%m%d %H%M%S")

    start_dt = np.array([sdt + dt.timedelta(seconds = t) for t in start_t])
    
    timeframes = start_dt 
    
    return(timeframes)

def laser_shots(file, time_info, channel_info, isdark = False):

    """
    General:
        Collects the number of shots per timeframe and channel
        from the QA file
        
    Input:
        file: 
            An xarray dataset with all the QA file information

        time_info : 
            A pandas Dataframe with all the necessary metadata for each 
            measurement timeframe fetched from the QA file

        channel_info : 
            A pandas Dataframe with all the necessary metadata for each 
            channel fetched from the QA file

        isdark :
            A boolean value. If set to Treu, the dark signals are read instead
            
    Returns:
        
        shots : 
            A 2D xarray Dataarray with the shots per timeframe and channel, 
            extracted from the QA file. The dimensions should be (time, channel)
    
    """

    timeframes =  time_info.index.values
    
    channels = channel_info.index.values
    
    if not isdark:
        shots_arr = file.Laser_Shots.values
    else:
        shots_arr = file.Background_Shots.values
        
    
    shots = xr.DataArray(shots_arr, 
                         dims = ['time', 'channel'],
                         coords = [timeframes, channels])

    # Sort by time
    shots = shots.copy().sortby('time')

    return(shots)

def signals(file, time_info, channel_info, isdark = False):
 
    """
    General:
        Collects the signal information per timeframe, channel, and rangebin
        from the QA file
        
    Input:
        file: 
            An xarray dataset with all the QA file information

        time_info : 
            A pandas Dataframe with all the necessary metadata for each 
            measurement timeframe fetched from the QA file

        channel_info : 
            A pandas Dataframe with all the necessary metadata for each 
            channel fetched from the QA file
        
        isdark :
            A boolean value. If set to Treu, the dark signals are read instead
            
    Returns:
        
        signal : 
            A 3D xarray Dataarray with the measured signals, extracted 
            from the QA file. The dimensions should be (time, channel, bins/range)
    
    """
    if not isdark:
        signal_arr = file.Raw_Lidar_Data.values
    else:
        signal_arr = file.Background_Profile.values
        
    bins = np.arange(0., signal_arr.shape[-1])
    
    channels = channel_info.index.values

    timeframes =  time_info.index.values
    
    signal = xr.DataArray(data = signal_arr, 
                       dims = ['time', 'channel', 'bins'],
                       coords = [timeframes, channels, bins])
    
    # Sort by time
    signal = signal.copy().sortby('time')
    signal = signal.copy().where(signal != nc.default_fillvals['f8'])
    
    return(signal)

def dark(file, time_info, channel_info):
    
    """
    General:
        Collects the signal information per timeframe, channel, and rangebin
        from the QA file
        
    Input:
        file: 
            An xarray dataset with all the QA file information

        time_info : 
            A pandas Dataframe with all the necessary metadata for each 
            measurement timeframe fetched from the QA file

        channel_info : 
            A pandas Dataframe with all the necessary metadata for each 
            channel fetched from the QA file
            
    Returns:
        
        signal_d : 
            A 3D xarray Dataarray with the measured dark signals, extracted 
            from the QA file. The dimensions should be (time, channel, bins/range)
            Defaults to an empty list if a dark measurement is not provided
    
    """
    
    signal_arr = file.Background_Profile.values

    bins = np.arange(0., signal_arr.shape[-1])
    
    channels = channel_info.index.values
    
    timeframes =  time_info.index.values
    
    signal = xr.DataArray(data = signal_arr, 
                       dims = ['time', 'channel', 'bins'],
                       coords = [timeframes, channels, bins])
    
    # Sort by time
    signal = signal.copy().sortby('time')
    signal = signal.copy().where(signal != nc.default_fillvals['f8'])

    return(signal)
