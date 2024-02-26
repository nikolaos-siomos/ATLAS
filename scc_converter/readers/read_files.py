"""
@author: Peristera
"""
import os, sys, glob
import numpy as np
from ..readers import read_licel, read_polly_xt, read_licel_matlab, read_polly_xt_first, read_licel_old2rack
import xarray as xr
import pandas as pd
from datetime import datetime
import re

def rayleigh(finput_ray, file_format):
    
    """Extracts the raw signal, shots, and rest metadata information out of the 
    raw input files. The default format is currently licel. The signal units
    are always mV for analog and counts for photon channels"""
    
    # Reading
    print('-----------------------------------------')
    print('Start reading Rayleigh signals...')
    print('-----------------------------------------')
    
    # Select reader based on the file format
    if file_format == 'polly_xt':
        system_info, channel_info, time_info, sig, shots = \
            read_polly_xt.dtfs(dir_meas = finput_ray, meas_type = 'ray')
            
    elif file_format == 'polly_xt_first':
        system_info, channel_info, time_info, sig, shots = \
            read_polly_xt_first.dtfs(dir_meas = finput_ray)
            
    elif file_format == 'licel':
        system_info, channel_info, time_info, sig, shots = \
            read_licel.dtfs(dir_meas = finput_ray)
            
    elif file_format == 'licel_matlab':
        system_info, channel_info, time_info, sig, shots = \
            read_licel_matlab.dtfs(dir_meas = finput_ray)

    elif file_format == 'licel_old2rack':
        system_info, channel_info, time_info, sig, shots = \
            read_licel_old2rack.dtfs(dir_meas = finput_ray)


    print('Reading Rayleigh signals complete!')
    print('-----------------------------------------')
    print('')

    return(sig, shots, system_info, channel_info, time_info)

def telecover(finput_sec, finput_rin, file_format, 
              files_per_sector = None, files_per_ring = None):
    
    """Extracts the raw signal, shots, and rest metadata information out of the 
    raw input files. The default format is currently licel. The signal units
    are always mV for analog and counts for photon channels"""
    
    # Reading
    print('-----------------------------------------')
    print('Start reading telecover signals...')
    print('-----------------------------------------')
    
    sig_sec = []
    sector_sec = []
    shots_sec = []
    time_info_sec = []
    
    sig = []
    shots = []
    time_info = []
    channel_info = []
    system_info = []
    
    if files_per_sector == None:
        if os.path.exists(finput_sec): 
            if os.listdir(finput_sec):
                for sector in ['north', 'east', 'south', 'west']:
                    path = os.path.join(finput_sec, sector)
                    print(f'-- Reading {sector} sector..')           

                    # Select reader based on the file format
                    if file_format == 'polly_xt':
                        system_info, channel_info, time_info, sig, shots = \
                            read_polly_xt.dtfs(dir_meas = path, meas_type = 'tlc')
                    
                    elif file_format == 'polly_xt_first':
                        system_info, channel_info, time_info, sig, shots = \
                            read_polly_xt_first.dtfs(dir_meas = path)
                          
                    elif file_format == 'licel':
                        system_info, channel_info, time_info, sig, shots = \
                            read_licel.dtfs(dir_meas = path)
                            
                    elif file_format == 'licel_matlab':
                        system_info, channel_info, time_info, sig, shots = \
                            read_licel_matlab.dtfs(dir_meas = path)

                    elif file_format == 'licel_old2rack':
                        system_info, channel_info, time_info, sig, shots = \
                            read_licel_old2rack.dtfs(dir_meas = path)

                    sector = folder_to_sector(folder = time_info['folder'].values)
    
                    time_info['sector'] = sector
                        
                    sig_sec.append(sig)
                    shots_sec.append(shots)
                    sector_sec.append(sector)
                    time_info_sec.append(time_info)
                
                sig = xr.concat(sig_sec, dim = 'time').sortby('time')
                shots = xr.concat(shots_sec, dim = 'time').sortby('time')
                time_info = pd.concat(time_info_sec).sort_index()
            else:
                print(f'-- Warning: Folder {finput_sec} is empty! No files to read ')
                
    else:
        if os.path.exists(finput_sec):
            if os.listdir(finput_sec):
                print('-- Reading sectors..')           
                
                # Select reader based on the file format
                if file_format == 'polly_xt':
                    system_info, channel_info, time_info, sig, shots = \
                        read_polly_xt.dtfs(dir_meas = finput_sec, meas_type = 'tlc')
                          
                elif file_format == 'polly_xt_first':
                    system_info, channel_info, time_info, sig, shots = \
                        read_polly_xt_first.dtfs(dir_meas = finput_sec)
                        
                elif file_format == 'licel':
                    system_info, channel_info, time_info, sig, shots = \
                        read_licel.dtfs(dir_meas = finput_sec)
                
                elif file_format == 'licel_matlab':
                    system_info, channel_info, time_info, sig, shots = \
                        read_licel_matlab.dtfs(dir_meas = finput_sec)
                        
                elif file_format == 'licel_old2rack':
                    system_info, channel_info, time_info, sig, shots = \
                        read_licel_old2rack.dtfs(dir_meas = finput_sec)
                    
                sector = time_to_sector(folder = time_info['folder'], 
                                        files_per_sector = files_per_sector)
                
                time_info['sector'] = sector
    
                sig_sec.append(sig)
                shots_sec.append(shots)
                sector_sec.append(sector)
                time_info_sec.append(time_info)
    
                sig = xr.concat(sig_sec, dim = 'time').sortby('time')
                shots = xr.concat(shots_sec, dim = 'time').sortby('time')
                time_info = pd.concat(time_info_sec).sort_index()
            else:
                print(f'-- Warning: Folder {finput_sec} is empty! No files to read ')
                

    if files_per_ring == None:
        if os.path.exists(finput_rin):
            if os.listdir(finput_rin):
                for ring in ['inner', 'outer']:
                    path = os.path.join(finput_rin, ring)
                    print(f'-- Reading {ring} rings..')           
                    # Select reader based on the file format
                    if file_format == 'polly_xt':
                        system_info, channel_info, time_info, sig, shots = \
                            read_polly_xt.dtfs(dir_meas = path, meas_type = 'tlc')
                          
                    elif file_format == 'polly_xt_first':
                        system_info, channel_info, time_info, sig, shots = \
                            read_polly_xt_first.dtfs(dir_meas = path)
                            
                    elif file_format == 'licel':
                        system_info, channel_info, time_info, sig, shots = \
                            read_licel.dtfs(dir_meas = path)
                            
                    elif file_format == 'licel_matlab':
                        system_info, channel_info, time_info, sig, shots = \
                            read_licel_matlab.dtfs(dir_meas = path)
                            
                    elif file_format == 'licel_old2rack':
                        system_info, channel_info, time_info, sig, shots = \
                            read_licel_old2rack.dtfs(dir_meas = path)
                
                    ring = folder_to_sector(folder = time_info['folder'].values)
                    time_info['sector'] = ring
                        
                    sig_sec.append(sig)
                    shots_sec.append(shots)
                    sector_sec.append(ring)
                    time_info_sec.append(time_info)
        
                sig = xr.concat(sig_sec, dim = 'time').sortby('time')
                shots = xr.concat(shots_sec, dim = 'time').sortby('time')
                time_info = pd.concat(time_info_sec).sort_index()
            else:
                print(f'-- Warning: Folder {finput_rin} is empty! No files to read ')
        
    else:
        if os.path.exists(finput_rin):
            if os.listdir(finput_rin):
                print('-- Reading rings..')    
                
                if file_format == 'polly_xt':
                    system_info, channel_info, time_info, sig, shots = \
                        read_polly_xt.dtfs(dir_meas = finput_rin, meas_type = 'tlc')

                elif file_format == 'polly_xt_first':
                    system_info, channel_info, time_info, sig, shots = \
                        read_polly_xt_first.dtfs(dir_meas = finput_rin)
                        
                elif file_format == 'licel':
                    system_info, channel_info, time_info, sig, shots = \
                        read_licel.dtfs(dir_meas = finput_rin)
                
                elif file_format == 'licel_matlab':
                    system_info, channel_info, time_info, sig, shots = \
                        read_licel_matlab.dtfs(dir_meas = finput_rin)
                        
                elif file_format == 'licel_old2rack':
                    system_info, channel_info, time_info, sig, shots = \
                        read_licel_old2rack.dtfs(dir_meas = finput_rin)                        
                        
                ring = time_to_ring(folder = time_info['folder'], 
                                    files_per_ring = files_per_ring)
            
                time_info['sector'] = ring
        
                sig_sec.append(sig)
                shots_sec.append(shots)
                sector_sec.append(ring)
                time_info_sec.append(time_info)
    
                sig = xr.concat(sig_sec, dim = 'time').sortby('time')
                shots = xr.concat(shots_sec, dim = 'time').sortby('time')
                time_info = pd.concat(time_info_sec).sort_index()

            else:
                print(f'-- Warning: Folder {finput_rin} is empty! No files to read ')
        
    print('Reading telecover signals complete!')
    print('-----------------------------------------')
    print('')

    return(sig, shots, system_info, channel_info, time_info)

def polarization_calibration(finput_ray, finput_p45, finput_m45, finput_stc, file_format):
    
    """Extracts the raw signal, shots, and rest metadata information out of the 
    raw input files. The default format is currently licel. The signal units
    are always mV for analog and counts for photon channels"""
    
    # Reading
    print('-----------------------------------------')
    print('Start reading Polarization Calibration signals...')
    print('-----------------------------------------')
            
    sig_pos = []
    shots_pos = []
    position_pos = []
    time_info_pos = []

    sig = []
    shots = []
    time_info = []
    channel_info = []
    system_info = []
    
    if file_format == 'licel':
        if os.path.exists(finput_stc):
            if os.listdir(finput_stc):
                print('-- Reading static calibration files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_licel.dtfs(dir_meas = finput_stc)
            
                position = np.array(time_info.index.size * [0])
                time_info['position'] = position
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_stc} is empty! No files to read ')
                       
        if os.path.exists(finput_m45):
            if os.listdir(finput_m45):
                print('-- Reading -45 files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_licel.dtfs(dir_meas = finput_m45)
    
                position = np.array(time_info.index.size * [1])
                time_info['position'] = position
        
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_m45} is empty! No files to read ')
                    

        if os.path.exists(finput_p45):
            if os.listdir(finput_p45):
                print('-- Reading +45 files..')  
  
                system_info, channel_info, time_info, sig, shots = \
                    read_licel.dtfs(dir_meas = finput_p45)
    
                position = np.array(time_info.index.size * [2])
                time_info['position'] = position
                
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_p45} is empty! No files to read ')
                 
        sig = xr.concat(sig_pos, dim = 'time').sortby('time')
        shots = xr.concat(shots_pos, dim = 'time').sortby('time')
        time_info = pd.concat(time_info_pos).sort_index()

    elif file_format == 'licel_matlab':
        if os.path.exists(finput_stc):
            if os.listdir(finput_stc):
                print('-- Reading static calibration files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_licel_matlab.dtfs(dir_meas = finput_stc)
            
                position = np.array(time_info.index.size * [0])
                time_info['position'] = position
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_stc} is empty! No files to read ')
                       
        if os.path.exists(finput_m45):
            if os.listdir(finput_m45):
                print('-- Reading -45 files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_licel_matlab.dtfs(dir_meas = finput_m45)
    
                position = np.array(time_info.index.size * [1])
                time_info['position'] = position
        
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_m45} is empty! No files to read ')
                    

        if os.path.exists(finput_p45):
            if os.listdir(finput_p45):
                print('-- Reading +45 files..')  
  
                system_info, channel_info, time_info, sig, shots = \
                    read_licel_matlab.dtfs(dir_meas = finput_p45)
    
                position = np.array(time_info.index.size * [2])
                time_info['position'] = position
                
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_p45} is empty! No files to read ')
        
        if len(sig_pos) != 0:        
            sig = xr.concat(sig_pos, dim = 'time').sortby('time')
            shots = xr.concat(shots_pos, dim = 'time').sortby('time')
            time_info = pd.concat(time_info_pos).sort_index()

    elif file_format == 'licel_old2rack':
        if os.path.exists(finput_stc):
            if os.listdir(finput_stc):
                print('-- Reading static calibration files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_licel_old2rack.dtfs(dir_meas = finput_stc)
            
                position = np.array(time_info.index.size * [0])
                time_info['position'] = position
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_stc} is empty! No files to read ')
                       
        if os.path.exists(finput_m45):
            if os.listdir(finput_m45):
                print('-- Reading -45 files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_licel_old2rack.dtfs(dir_meas = finput_m45)
    
                position = np.array(time_info.index.size * [1])
                time_info['position'] = position
        
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_m45} is empty! No files to read ')
                    
    
        if os.path.exists(finput_p45):
            if os.listdir(finput_p45):
                print('-- Reading +45 files..')  
    
                system_info, channel_info, time_info, sig, shots = \
                    read_licel_old2rack.dtfs(dir_meas = finput_p45)
    
                position = np.array(time_info.index.size * [2])
                time_info['position'] = position
                
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_p45} is empty! No files to read ')
        
        if len(sig_pos) != 0:        
            sig = xr.concat(sig_pos, dim = 'time').sortby('time')
            shots = xr.concat(shots_pos, dim = 'time').sortby('time')
            time_info = pd.concat(time_info_pos).sort_index()

    elif file_format == 'polly_xt_first':
        if os.path.exists(finput_stc):
            if os.listdir(finput_stc):
                print('-- Reading static calibration files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_polly_xt_first.dtfs(dir_meas = finput_stc)
            
                position = np.array(time_info.index.size * [0])
                time_info['position'] = position
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_stc} is empty! No files to read ')
                       
        if os.path.exists(finput_m45):
            if os.listdir(finput_m45):
                print('-- Reading -45 files..')  
                
                system_info, channel_info, time_info, sig, shots = \
                    read_polly_xt_first.dtfs(dir_meas = finput_m45)
    
                position = np.array(time_info.index.size * [1])
                time_info['position'] = position
        
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_m45} is empty! No files to read ')
                    

        if os.path.exists(finput_p45):
            if os.listdir(finput_p45):
                print('-- Reading +45 files..')  
  
                system_info, channel_info, time_info, sig, shots = \
                    read_polly_xt_first.dtfs(dir_meas = finput_p45)
    
                position = np.array(time_info.index.size * [2])
                time_info['position'] = position
                
                sig_pos.append(sig)
                shots_pos.append(shots)
                position_pos.append(position)
                time_info_pos.append(time_info)
            else:
                print(f'-- Warning: Folder {finput_p45} is empty! No files to read ')
        
        if len(sig_pos) != 0:        
            sig = xr.concat(sig_pos, dim = 'time').sortby('time')
            shots = xr.concat(shots_pos, dim = 'time').sortby('time')
            time_info = pd.concat(time_info_pos).sort_index()

    # Select reader based on the file format
    elif file_format == 'polly_xt':
        if os.path.exists(finput_ray):
            if os.listdir(finput_ray):
                system_info, channel_info, time_info, sig, shots = \
                    read_polly_xt.dtfs(dir_meas = finput_ray, meas_type = 'pcb')
            else:
                print(f'-- Warning: Folder {finput_ray} is empty! No files to read ')
                         
            
    print('Reading Polarization Calibration signals complete!')
    print('-----------------------------------------')
    print('')

    return(sig, shots, system_info, channel_info, time_info)


def dark(finput_drk, file_format):
    
    """Extracts the raw signal, shots, and rest metadata information out of the 
    raw input files. The default format is currently licel. The signal units
    are always mV for analog and counts for photon channels"""
    
    # Reading
    print('-----------------------------------------')
    print('Start reading dark signals...')
    print('-----------------------------------------')
    
    # Select reader based on the file format
    if file_format == 'polly_xt':
        system_info, channel_info, time_info, sig, shots = \
            read_polly_xt.dtfs(dir_meas = finput_drk, meas_type = 'drk')
    elif file_format == 'polly_xt_first':
        system_info, channel_info, time_info, sig, shots = \
            read_polly_xt_first.dtfs(dir_meas = finput_drk)
    elif file_format == 'licel':
        system_info, channel_info, time_info, sig, shots = \
            read_licel.dtfs(dir_meas = finput_drk)
    elif file_format == 'licel_matlab':
        system_info, channel_info, time_info, sig, shots = \
            read_licel_matlab.dtfs(dir_meas = finput_drk)
    elif file_format == 'licel_old2rack':
        system_info, channel_info, time_info, sig, shots = \
            read_licel_old2rack.dtfs(dir_meas = finput_drk)
            
    print('Reading dark signals complete!')
    print('-----------------------------------------')
    print('')
  
    return(sig, shots, system_info, channel_info, time_info)


def radiosonde(finput_rs, delimiter, skip_header, skip_footer, 
               usecols, units, mtime, ground):

    """Extracts the meteorological information out of the 
    raw radiosonde file."""
    
    # Reading
    print('-----------------------------------------')
    print('Start reading radiosonde file...')
    print('-----------------------------------------')
    
    paths = glob.glob(os.path.join(finput_rs,'*_*.txt'))
    
    if len(paths) == 0:
        raise Exception(f"-- Error: No txt file was found in the radiosonde folder: {finput_rs} Please make sure that the radiosonde files are in txt format")
    
    lib_delimiter =  {"S": "",
                      "C": ",",
                      "T": "\t"}
    
    # Unit conversion functions
    def km_asl_to_m_asl(x):
        return(1E3 * x)

    def m_agl_to_m_asl(x, ground = 0.):
        return(x + ground)

    def km_agl_to_m_asl(x, ground = 0.):
        return(1E3 * x + ground)
    
    def geo_to_asl(x):
        Re = 6.371E6
        return(x * Re / (Re - x))
    
    def Pa_to_hPa(x):
        return(1E-2 * x)
    
    def atm_to_hPa(x):
        return(x * 1013.25)
    
    def C_to_K(x):
        return(x + 273.16)

    def Cx10_to_K(x):
        return(x/10. + 273.16)
    
    def fraction_to_percent(x):
        return(100. * x)
    
    if len(paths) == 0 :
        raise Exception("-- Error: No txt file provided in the radiosonde folder! Please provide a single file with the radiosonde data with a filename that starts with 'yyyymmdd_hhmm' and ends with '.txt' ")

    bname = [os.path.basename(path) for path in paths]

    bad_length = [len(name) < 14 for name in bname]
    
    if any(bad_length) :
        raise Exception(f"-- Error: Radiosonde filename with wrong length detected! Please revise the following files: {bname[bad_length]}. They should start with 'yyyymmdd_hhmm' and end with '.txt' ")
    else:
        pattern = "20[0-9]{2}[0-1][0-9][0-3][0-9]_[0-2][0-9][0-5][0-9]"
        bad_format = [not(bool(re.search(pattern,name))) for name in bname]

        if any(bad_format):
            raise Exception(f"-- Error: Radiosonde filename with wrong format detected! Please revise the following files: {bname[bad_format]}. They should start with 'yyyymmdd_hhmm' and end with '.txt' ")
        else:
            dates = [name[:13].split('_')[0] for name in bname]
            times = [name[:13].split('_')[1] for name in bname]

    bad_dates = [int(date[:4]) not in np.arange(1960,9999,1) or \
                 int(date[4:6]) not in np.arange(1,13,1) or \
                 int(date[6:8]) not in np.arange(1,32,1) for date in dates]

    bad_times = [int(time[:2]) not in np.arange(0,24,1) or \
                 int(time[2:4]) not in np.arange(0,60,1) for time in times]
        
    if any(bad_dates):
        raise Exception(f"-- Error: The date provided in at least one radiosonde filename is not correct. Please revise the following files: {np.array(bname)[bad_dates]}. It should start with 'yyyymmdd_hhmm' and end with '.txt' ")

    if any(bad_times):
        raise Exception(f"-- Error: The time provided in at least radiosond filename is not correct. Please revise the following files: {np.array(bname)[bad_times]}. It should start with 'yyyymmdd_hhmm' and end with '.txt' ")
        
    date_dt = np.array([datetime.strptime(date + time,'%Y%m%d%H%M') 
                        for date, time in zip(dates, times)])
    
    delta_t = np.array([(dt - mtime).total_seconds() /3600. for dt in date_dt])

    ind_rs = np.argmin(np.abs(delta_t))
    
    if not any(np.abs(delta_t) < 24):
        raise Exception(f"-- Error: The nearest radiosonde in time {bname[ind_rs]} was launched with a time difference of {np.round(delta_t[ind_rs],decimals=1)} hours with respect to the middle time of the measurement! Please provide a radiosond file with less than 18 hours temporal difference")
    else:
        print(f'-- Selected radiosonde file: {bname[ind_rs]}')
        
    if usecols[3] == None:
        parameters = ['P', 'T']
        usecols = usecols[:3]

    else:
        parameters = ['P', 'T', 'RH']
        
    data = np.genfromtxt(paths[ind_rs],skip_header = skip_header, 
                         skip_footer = skip_footer,
                         delimiter = lib_delimiter[delimiter], 
                         autostrip = True,
                         usecols = np.array(usecols) - 1, dtype = float)

    if units[0] == 'km_asl':
        data[:,0] = km_asl_to_m_asl(data[:,0])

    if units[0] in ['m_agl', 'km_agl']:
        if ground == None:
            raise Exception("-- Error: The altitude parameter of the rsonde_geodata field is mandatory when the radiosond height is in agl (altitude above ground units). Please provide at least 1 float corresponding to the station altitude: --rsonde_geodata 60.0")
        else:
            if units[0] == 'm_agl':
                data[:,0] = m_agl_to_m_asl(data[:,0], ground = ground)
            else:
                data[:,0] = km_agl_to_m_asl(data[:,0], ground = ground)

    # if 'geo' in units[0]:
    #     data[:,0] = geo_to_asl(data[:,0])
        
    if units[1] == 'Pa':
        data[:,1] = Pa_to_hPa(data[:,1])    

    if units[1] == 'atm':
        data[:,1] = atm_to_hPa(data[:,1])   
        
    if units[2] == 'C':
        data[:,2] = C_to_K(data[:,2]) 

    if units[2] == 'Cx10':
        data[:,2] = Cx10_to_K(data[:,2]) 
        
    if units[3] == 'fraction' and len(usecols) == 4:
        data[:,3] = fraction_to_percent(data[:,3])     
        
    alt = data[:,0]
    
    atmo = xr.DataArray(data[:,1:], 
                        coords = [alt, parameters], 
                        dims = ['height', 'parameters'] )
    
    return(dates[ind_rs], times[ind_rs], atmo)

def folder_to_sector(folder):

    fld = ['north','east','south','west','outer','inner']
    sec = [1,2,3,4,5,6]
            
    sector = np.nan * np.zeros(folder.shape)
    
    for i in range(len(fld)):
        sector[folder == fld[i]] = sec[i]
        
    return(sector)

def folder_to_position(folder):

    fld = ['static', '-45', '+45']
    sec = [0, 1, 2]
            
    position = np.nan * np.zeros(folder.shape)
    
    for i in range(len(fld)):
        position[folder == fld[i]] = sec[i]
        
    return(position)

def time_to_sector(folder, files_per_sector):
    
    blocks = folder.size / files_per_sector
    
    if blocks - np.floor(blocks) > 0.:
        raise Exception("-- Error: The files_per_sector argument was provided but " +
                 "the number of telecover files cannot be evenly divided by it! " +
                 "Please revise the files_per_sector value. If the number of " +
                 "files per sector was not constant during measurements then " +
                 "provide the telecover in individual folders per sector.")
    
    sec = [1, 2, 3, 4]
    
    sec_list = int(np.floor(blocks / 4.)) * sec
    sec_list.extend(sec[:((blocks - np.floor(blocks / 4.)) * 4.).astype(int)])
            
    sector = np.nan * np.zeros(folder.shape)
    
    for i in range(len(sec_list)):
        sector[i*files_per_sector:(i+1)*files_per_sector] = sec_list[i]
    
    return(sector)

def time_to_ring(folder, files_per_ring):
    
    blocks = folder.size / files_per_ring
    
    if blocks - np.floor(blocks) > 0.:
        raise Exception("-- Error: The files_per_ring argument was provided but " +
                 "the number of telecover files cannot be evenly divided by it! " +
                 "Please revise the files_per_ring value. If the number of " +
                 "files per ring was not constant during measurements then " +
                 "provide the telecover in individual folders per ring.")
    
    sec = [5, 6]
    
    sec_list = int(np.floor(blocks / 2.)) * sec
    sec_list.extend(sec[:((blocks - np.floor(blocks / 2.)) * 2.).astype(int)])
            
    ring = np.nan * np.zeros(folder.shape)
    
    for i in range(len(sec_list)):
        ring[i*files_per_ring:(i+1)*files_per_ring] = sec_list[i]
    
    return(ring)