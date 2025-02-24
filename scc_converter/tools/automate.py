"""
@author: Peristera
"""
import os, sys
import numpy as np
import xarray as xr
import pandas as pd

def get_meas_type(args):

    """Identifies if the measurement being processed is a Rayleigh, 
    a Telecover, a Depolarization Calibration, or a standalone Dark."""
    
    meas_type = []
    
    dirs = {'drk' : args['dark_folder'],
            'ray' : args['rayleigh_folder'],
            'tlc_sec' : args['telecover_sectors_folder'],
            'tlc_rin' : args['telecover_rings_folder'],
            'pcb_stc' : args['pcb_cal_stc_folder'],
            'pcb_p45' : args['pcb_cal_p45_folder'],
            'pcb_m45' : args['pcb_cal_m45_folder']}                                        
        
    if os.path.exists(dirs['ray']): 
        meas_type.append('rayleigh')

    if os.path.exists(dirs['tlc_sec']) or os.path.exists(dirs['tlc_rin']): 
        meas_type.append('telecover')

    if (os.path.exists(dirs['pcb_p45']) and os.path.exists(dirs['pcb_m45'])) \
        or os.path.exists(dirs['pcb_stc']):
        meas_type.append('polarization_calibration')

    if os.path.exists(dirs['drk']):
        meas_type.append('dark')
    
    if len(meas_type) == 0:
        raise Exception('-- Error: None of the input folders exists!')

    return(meas_type)


def check_telecover_sec(path, files_per_sector):

    """Ensures that the telecover folder is properly set, otherwise it raises 
    an error """
        
    if os.path.exists(path) and files_per_sector == None:
    
        allowed_folders = np.array(['north', 'east', 'south', 'west'])

        dirs = np.array([dir_i for dir_i in os.listdir(path)
                         if os.path.isdir(os.path.join(path,dir_i))])

        unk_folder = np.array([dir_i not in allowed_folders 
                               for dir_i in dirs])

        mis_folder = np.array([dir_i not in dirs 
                               for dir_i in allowed_folders])
        
        if any(unk_folder):
                
            raise Exception('-- Error: Unknown subfolder(s) detected  ' +
                f'in the telecover sectors folder: {dirs[unk_folder]} ' +
                f'in {path} Please remove in order to proceed!')

        if any(mis_folder):
                    
            raise Exception('-- Error: Subfolder(s) missing ' +
                f'in the telecover sectors folder: {allowed_folders[mis_folder]} ' +
                f'Please either provide all: {allowed_folders} ' +
                'folders or use the files_per_sector '+
                'argument as: --files_per_sector <number_of_files>')

    elif os.path.exists(path) and files_per_sector != None:
        dirs = np.array([dir_i for dir_i in os.listdir(path)
                         if os.path.isdir(os.path.join(path,dir_i))])
        if len(dirs) > 0:
            raise Exception(f'Directories where detected inside the telecover folder ({path}) even though files_per_sector option was provided. Please either don;t provide files_per_sector at all(default) and create the necessary telecover subfolder inside the telecover folder or put all telecover files directly insided the telcover folder without any subfolders ')

    return()

def check_telecover_rin(path, files_per_ring):

    """Ensures that the telecover folder is properly set, otherwise it raises 
    an error """
        
    if os.path.exists(path) and files_per_ring == None:
    
        allowed_folders = ['inner', 'outer']

        dirs = np.array([dir_i for dir_i in os.listdir(path)
                         if os.path.isdir(os.path.join(path,dir_i))])

        unk_folder = np.array([dir_i not in allowed_folders 
                               for dir_i in dirs])

        mis_folder = np.array([dir_i not in dirs 
                               for dir_i in allowed_folders])
        
        if any(unk_folder):
                
            raise Exception('-- Error: Unknown subfolder(s) detected  ' +
                f'in the telecover rings folder: {dirs[unk_folder]} ' +
                f'in {path} Please remove in order to proceed!')

        if any(mis_folder):
                    
            raise Exception('-- Error: Subfolder(s) missing ' +
                'in the telecover rings folder: allowed_folders{mis_folder} ' +
                f'Please either provide all: {allowed_folders} ' +
                'folders or use the files_per_ring '+
                'argument as: --files_per_ring <number_of_files>')
            
    elif os.path.exists(path) and files_per_ring != None:
        dirs = np.array([dir_i for dir_i in os.listdir(path)
                         if os.path.isdir(os.path.join(path,dir_i))])
        if len(dirs) > 0:
            raise Exception(f'Directories where detected inside the telecover folder ({path}) even though files_per_sector option was provided. Please either don;t provide files_per_sector at all(default) and create the necessary telecover subfolder inside the telecover folder or put all telecover files directly insided the telcover folder without any subfolders ')

    return()


def detect_overflows(sig, shots, channel_info, time_info, meas_type, method = 0):

    """
    General:
        Detects and removes lidar profiles with photon values above the
        maximum allowed countrate or analog values above the data 
        acquisition range
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, bins). 

        shots:
            A 2D xarray with the laser shots per timeframe and channel, 
            it should include the following dimensions: (time, channel).  
            
        channel_info:
            A pandas Dataframe with exactly one index entry per lidar channel.
            the following columns must be included:
                
            acquisition_mode: 
                The acquision_mode values per channel 
                (0 for analog, 1 for photon)
                
            data_acquisition_range: 
                The data acquisition range of the analog channels.

        time_info: 
            A pandas dataframe with exactly one index entry per 
            lidar timeframe. The index should correspond to the time 
            dimension of sig. The following column must be included         
            
            filename: 
                The raw file filename. 
        
        meas_type: 
            A 3 letter identifier that specifies the measurement type,
            it can be one of ray, tlc, pcb, drk        
        
        method:
            An integer. If set to 0 only the check for overflows will be 
            performed and an error will be raised if a single value is 
            encountered. If set to 1 all files with overflows will be
            discarded. If set to 2 overflows will be removed and 
            interpolated instead
            
            
    Returns:
        
        sig, shots, time:
            If no overflows are detected then the 3 variables reuturn intact
            
            If overflows are detected:
                method = 0 : the code exits with a diagnostic error
            
                method = 1 : timeframes with at least one overflowed bin are
                             removed
                
                method = 2 : only sig changes, overflowed values are replaced
                             by interpolated values from the surrounding bins

                method = 3 : do nothing about it, use only while debugging            
    """
    
    meas_label = {'ray': 'Rayleigh', 'tlc' : 'Telecover',
                  'pcb': 'Polarization Calibration', 'drk': 'Dark'}
    
    print('-----------------------------------------')
    print(f'Handling {meas_label[meas_type]} overflow values...')
    print('-----------------------------------------')
    
    acquisition_mode = channel_info.acquisition_mode
    
    daq_range = channel_info.data_acquisition_range

    filename = time_info.filename
        
    # Get an overflow mask for each bin
    mask = get_overflow_mask(sig, acquisition_mode, daq_range)

    if method == 0 and mask.any(): # Detect the problematic profiles and raise error

        overflow_method_0(mask = mask, filename = filename)
            
    elif method == 1 and mask.any(): # Remove the problematic profiles
    
        sig, shots, time_info = overflow_method_1(sig = sig.copy(), 
                                                  shots = shots.copy(),
                                                  time_info = time_info, 
                                                  mask = mask, 
                                                  filename = filename)
        print("Profiles with overflowed bin were succesfully removed!")
        
    elif method == 2 and mask.any(): # Replace the overflowed values with interpolated ones from the nearby bins
        
        sig = overflow_method_2(sig = sig.copy(), 
                                mask = mask, 
                                filename = filename)
        print("Overflowed were succesfully replaced!")
    
    elif method == 3 and mask.any():

        print("-- Warning: Overflows were detected but no action has been performed! Use trim_overflows to 3 only when debugging!")

    elif not mask.any():
        
        print("No bins with overflows have been encountered!")
    
    print('-----------------------------------------')
    print("")
     

    return(sig, shots, time_info)

def get_overflow_mask(sig, acquisition_mode, daq_range):

    channels = sig.channel.values
    
    time = sig.time.values

    bins = sig.bins.values
    
    mask = xr.DataArray(np.zeros([time.size, channels.size, bins.size], 
                                 dtype = bool),
                        dims = ['time','channel', 'bins'],
                        coords = [time, channels, bins])
    
    for ch in channels:

        ch_d = dict(channel = ch)

        if acquisition_mode.loc[ch] == 1: #3rd digit of channel name is the acquisition mode (a or p)
        
            max_count = np.power(2.,15)

            crit = (sig.loc[ch_d] >= max_count)
            
            if crit.any():
                print(f"-- Warning: Channel {ch} - Photon signal count values above the maximum allowed summed counts were detected! ")

                mask.loc[ch_d] = crit.values                                           

        if acquisition_mode.loc[ch] == 0: #3rd digit of channel name is the acquisition mode (a or p)

            max_mV = daq_range.loc[ch]
        
            crit = (sig.loc[ch_d] >= max_mV) | (sig.loc[ch_d] <= 0.)
            
            if crit.any():
                print(f"-- Warning: Channel {ch} - Analog signal mV values above the data acqusition range or below 0. were detected! ")

                mask.loc[ch_d] = crit.values
    
    return(mask)

def overflow_method_0(mask, filename):
    
    time = mask.time.values
    
    channels = mask.channel.values
    
    bins = mask.bins.values
    
    # use the index and not the time....
    if mask.any():
        
        print("-- Error at least one bin with an overflow was detected. ")
        print("-- Please revise the following bins: ")
        
        mask_t = mask.any(dim = 'bins').any(dim = 'channel').values
                    
        time_ovf = time[mask_t]
                    
        for t in time_ovf:
            filename_ovf = filename[t]
            
            t_d = dict(time = t)
                            
            ch_ovf = channels[mask.any(dim = 'bins').loc[t_d].values]
            
            for ch in ch_ovf:
            
                ch_d = dict(channel = ch)
                
                bins_ovf = bins[mask.loc[t_d].loc[ch_d].values].astype(int)
                
                print(f"    file: {filename_ovf} | ch: {ch} | bins: {bins_ovf}")
        
        raise Exception("-- Error: In order to continue with an automated overflow removal use the trim_overflow argument with value 1 or 2 (default is 0) ")
        
    return()

def overflow_method_1(sig, shots, time_info, mask, filename):

    time = mask.time.values
    
    mask_t = mask.any(dim = 'bins').any(dim = 'channel').values

    print(f"-- Warning: {np.sum(mask_t)} profiles with at least one overflowed bin have been detected ")
    print("-- Warning: trim_overflows = 1: The following profiles have been removed ")
                
    time_ovf = time[mask_t]

    time_cor = time[~mask_t]
                
    for t in time_ovf:
        filename_ovf = filename[t]
        
        print(f"    {filename_ovf} ")

    time_d = dict(time = time_cor)
    
    sig_out = sig.copy().loc[time_d] 

    shots_out = shots.copy().loc[time_d] 
    
    time_info = time_info.loc[time_cor,:]
    
    return(sig_out, shots_out, time_info)

def overflow_method_2(sig, mask, filename):

    mask_t = mask.any(dim = 'bins').any(dim = 'channel').values

    ovfs = mask.sum(dim = 'bins')

    print(f"-- Warning: {np.sum(mask_t)} profiles with at least one overflowed bin have been detected ")
    print("-- Warning: trim_overflows = 2: Will attempt to replace the overflows ")
          
    if (ovfs > 100).any():
        print("-- Warning: More that 100 overflowed bins encountered in single profiles... Interpolation is too risky, please revise the input files or consider switching to --trim_overflows 2 ")
    
    sig = sig.copy().where(~mask)
    
    sig_out = sig.copy().interpolate_na(dim = "bins", method = "linear")

    print(f"-- Warning: {np.sum(mask).values} overflows have been replaced by interpolating across the bins ")
    
    return(sig_out)
  