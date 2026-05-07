#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 11:17:36 2025

@author: nikos
"""

import os, sys, glob
import numpy as np
import xarray as xr
import re
import netCDF4 as nc
from utils.toolbox import get_mid_time, find_nearest_file
from molecular.utilities import number_density_at_pt, saturation_vapour_pressure
from utils.unit_conversions import km_asl_to_m_asl, m_agl_to_m_asl, \
    km_agl_to_m_asl, hPa_to_Pa, atm_to_Pa, C_to_K, Cx10_to_K, percent_to_fraction


default_parsing_options = {'delimiter' : 'C', 
                          'skip_header' : 4, 
                          'skip_footer' : 1, 
                          'usecols' : [2,1,3,5],
                          'units' : ["m_asl", "hPa", "C" , "percent"]}

lib_delimiter =  {"S": "",
                  "C": ",",
                  "T": "\t"}


def date_from_filename(input_file):
    
    rsond_start_date = input_file[:13].split('_')[0]
    rsonde_start_time = input_file[:13].split('_')[1]
    
    return(rsond_start_date, rsonde_start_time)


def get_filename_list(input_folder, pattern):
    
    paths = glob.glob(os.path.join(input_folder,pattern))

    if len(paths) == 0:
   
        file_validity_flag = "empty"
        
        print(f"-- Warning: No txt file was found in the wyoming radiosonde folder: {input_folder} Please make sure that the radiosonde files are in txt format")
    
    else:
        bnames = [os.path.basename(path) for path in paths]
        file_validity_flag = 'usable'

    return(bnames, file_validity_flag)
    
def convert_units(data, parsing_options = None, ground = None):
    
    units = default_parsing_options["units"]
    
    if units[0] == 'km_asl':
        data.height.values = km_asl_to_m_asl(data.height.values)

    if units[0] in ['m_agl', 'km_agl']:
        if ground == None:
            raise Exception("-- Error: The altitude parameter of the rsonde_geodata field is mandatory when the radiosond height is in agl (altitude above ground units). Please provide at least 1 float corresponding to the station altitude: --rsonde_geodata 60.0")
        else:
            if units[0] == 'm_agl':
                data.height.values = m_agl_to_m_asl(
                    data.height.values, 
                    ground = ground
                    )
            else:
                data[:,0] = km_agl_to_m_asl(
                    data.height.values, 
                    ground = ground
                    )

    # if 'geo' in units[0]:
    #     data[:,0] = geo_to_asl(data[:,0])
        
    if units[1] == 'hPa':
        data.loc['P',:] = hPa_to_Pa(data.loc['P',:])    

    if units[1] == 'atm':
        data.loc['P',:] = atm_to_Pa(data.loc['P',:])   
        
    if units[2] == 'C':
        data.loc['T',:] = C_to_K(data.loc['T',:]) 

    if units[2] == 'Cx10':
        data.loc['T',:] = Cx10_to_K(data.loc['T',:]) 
        
    if units[3] == 'percent':
        data.loc['RH',:] = percent_to_fraction(data.loc['RH',:])  
        
    return(data)
        
def read_radiosonde_ascii(input_file, mtime, parsing_options = None):

    if parsing_options is None:
        parsing_options = default_parsing_options
        
    # Unpack parsing options
    usecols = parsing_options['usecols']
    skip_header = parsing_options['skip_header']
    skip_footer = parsing_options['skip_footer']
    delimiter = parsing_options['delimiter']

    data = np.genfromtxt(input_file,
                         skip_header = skip_header, 
                         skip_footer = skip_footer,
                         delimiter = lib_delimiter[delimiter], 
                         autostrip = True,
                         usecols = np.array(usecols) - 1, dtype = float)
    
    parameters = ['P', 'T', 'RH']         
        
    height = data[:,0]
    
    meteo = xr.DataArray(data[:,1:].T, 
                        coords = [parameters, height], 
                        dims = ['parameters', 'height'])
    
    return(meteo)

def export_path(output_folder, meas_ID):
    """Creates a sting variable with the full path to the output QA file: 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html"""
    
    nc_path = os.path.join(output_folder,f'rs_{meas_ID[:-4]}00.nc')
    
    if os.path.exists(output_folder) == False:
        os.makedirs(nc_path)
        
    if os.path.exists(nc_path):
        os.unlink(nc_path)
        
    return(nc_path)

def make_nc_var(ds, name, value, dtype, dims = []):  
    """Function called by the *_file functions in order to fascilitate variable
    creation in the netcdf"""
    
    if dtype == 'int':
        func = np.int32
        default_val = nc.default_fillvals['i4']
        
    if dtype == 'float':
        func = np.double
        default_val = nc.default_fillvals['f8']
      
    if len(dims) == 0:
        value = func(value)
    else:
        value[value != value] = default_val
        value = value.astype(dtype)

    var = ds.createVariable(name, func, dims)
    
    if len(dims) == 0:
        var[:] = value

    elif len(dims) == 1:
        var[:] = value

    elif len(dims) == 2:
        var[:,:] = value
        
    elif len(dims) == 3:
        var[:,:,:] = value
        
    else:
        sys.exit('-- Error: 4 or higher dimensional arrays not supported in function make_nc_var')
        
    return()

def export_radiosonde_scc(nc_path, st_name, wmo_id, wban_id, 
                          date, time, ground, lat, lon, meteo):

    print('-----------------------------------------')
    print('Start exporting to a radiosonde QA file...')
    print('-----------------------------------------')
    
    """Creates the radiosonde netcdf file according to the SCC format 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html
    and exports it to nc_path"""

    n_points = meteo.height.size

    ds = nc.Dataset(nc_path,mode='w')

# Adding Dimensions
    ds.createDimension('points', n_points)
        
# Adding Global Parameters    
    ds.Altitude_meter_asl = ground;

    if lat != None:
        ds.Latitude_degrees_north = lat;

    if lon != None:
        ds.Longitude_degrees_east = lon;

    ds.Measurement_type = 'rs';

    ds.Sounding_Start_Date = date;

    ds.Sounding_Start_Time_UT = f'{time}00';

    if st_name != None:
        ds.Sounding_Station_Name = st_name;
        
    if wmo_id != None:
        ds.WMO_Station_Number = str(wmo_id);

    if wban_id != None:
        ds.WBAN_Station_Number = str(wban_id);

# Adding Variables
    make_nc_var(ds, name = 'Altitude', value = meteo.height.values, dtype = 'float', dims = ('points',))

    make_nc_var(ds, name = 'Pressure', value = meteo.loc[dict(parameters = 'P')].values, dtype = 'float', dims = ('points',))

    make_nc_var(ds, name = 'Temperature', value = meteo.loc[dict(parameters = 'T')].values, dtype = 'float', dims = ('points',))

    if 'RH' in meteo.parameters.values:
        make_nc_var(ds, name = 'RelativeHumidity', value = meteo.loc[dict(parameters = 'RH')].values, dtype = 'float', dims = ('points',))

    ds.close()
    
    return()



def export_to_scc(output_folder, meteo, metadata):
    
    date = metadata['rsonde_start_date']
    time = metadata['rsonde_start_time']
    station_id = metadata['station_id']
    
    # Creating the radiosonde ID
    rsonde_ID = f"{date}{station_id}{time[:4]}"
    
    # Creating the paths and folders
    nc_path = export_path(output_folder = output_folder, meas_ID = rsonde_ID)

    # Making the raw SCC file
    export_radiosonde_scc(nc_path = nc_path, 
                          date = metadata['rsonde_start_date'], 
                          time = metadata['rsonde_start_time'], 
                          ground = metadata['rsonde_altitude'], 
                          lat = metadata['rsonde_latitude'], 
                          lon = metadata['rsonde_longitude'], 
                          st_name = metadata['rsonde_station_name'], 
                          wmo_id = metadata['rsonde_wmo_number'], 
                          wban_id = metadata['rsonde_wban_number'], 
                          meteo = meteo)
    
    print('--Succesfully generated an SCC radiosonde file!')
    print('')
    
def evaluate_flag(file_validity_flag, label):
    
    if file_validity_flag == 'usable':
        print(f'-- Succesfully parsed {label} meteorological file!')
    elif file_validity_flag == 'unusable':
        print(f'-- The {label} meteorological file is not usable!')
    elif file_validity_flag == 'empty':
        pass
    else:
        print(f'-- Warning: File validity flag {file_validity_flag} not foreseen!')
        
def wyoming_filename_checks(bnames):
    
    bad_length = [len(name) < 14 for name in bnames]
    
    if any(bad_length) :
        raise Exception(f"-- Error: Radiosonde filename with wrong length detected! Please revise the following files: {bnames[bad_length]}. They should start with 'yyyymmdd_hhmm' and end with '.txt' ")
    else:
        pattern = "20[0-9]{2}[0-1][0-9][0-3][0-9]_[0-2][0-9][0-5][0-9]"
        bad_format = [not(bool(re.search(pattern,name))) for name in bnames]

        if any(bad_format):
            raise Exception(f"-- Error: Radiosonde filename with wrong format detected! Please revise the following files: {bnames[bad_format]}. They should start with 'yyyymmdd_hhmm' and end with '.txt' ")
        else:
            dates = [name[:13].split('_')[0] for name in bnames]
            times = [name[:13].split('_')[1] for name in bnames]

    bad_dates = [int(date[:4]) not in np.arange(1960,9999,1) or \
                 int(date[4:6]) not in np.arange(1,13,1) or \
                 int(date[6:8]) not in np.arange(1,32,1) for date in dates]

    bad_times = [int(time[:2]) not in np.arange(0,24,1) or \
                 int(time[2:4]) not in np.arange(0,60,1) for time in times]
        
    if any(bad_dates):
        raise Exception(f"-- Error: The date provided in at least one wyoming radiosonde filename is not correct. Please revise the following files: {np.array(bnames)[bad_dates]}. It should start with 'yyyymmdd_hhmm' and end with '.txt' ")

    if any(bad_times):
        raise Exception(f"-- Error: The time provided in at least one wyoming radiosond filename is not correct. Please revise the following files: {np.array(bnames)[bad_times]}. It should start with 'yyyymmdd_hhmm' and end with '.txt' ")
            
def add_number_density(meteo):
    
    height = meteo.height.values
    pressure = meteo.loc['P',:].values
    temperature = meteo.loc['T',:].values
    relative_humidity = meteo.loc['RH',:].values
    
    number_density = number_density_at_pt(
        pressure = pressure, 
        temperature = temperature, 
        relative_humidity = relative_humidity, 
        ideal=False
        )
    
    number_density = xr.DataArray([number_density], 
                                  coords = [['N'], height], 
                                  dims = ['parameters', 'height'])
    
    meteo = xr.concat([meteo, number_density], dim="parameters")
    
    return(meteo)

def add_saturation_vapour_pressure(meteo):
    
    height = meteo.height.values
    pressure = meteo.loc['P',:].values
    temperature = meteo.loc['T',:].values
    
    sat_vap_pressure = saturation_vapour_pressure(
        pressure = pressure, 
        temperature = temperature, 
        )
    
    sat_vap_pressure = xr.DataArray([sat_vap_pressure], 
                                  coords = [['e_s'], height], 
                                  dims = ['parameters', 'height'])
    
    meteo = xr.concat([meteo, sat_vap_pressure], dim="parameters")
    
    return(meteo)

def select_radiosonde_wyoming(input_folder, signal_times, metadata):
    
    """Extracts the meteorological information out of the 
    wyoming radiosonde file."""
    
    # Reading
    print('-----------------------------------------')
    print('Parsing Wyoming meteorological files...')
    print('-----------------------------------------')
    
    mtime = get_mid_time(signal_times)
    
    bnames, file_validity_flag = \
        get_filename_list(input_folder, pattern = '*_*.txt')
    
    if file_validity_flag != "empty":
        
        wyoming_filename_checks(bnames)
        
        input_file, file_validity_flag = find_nearest_file(
            mtime = mtime, 
            bnames = bnames,
            filetype = 'wyoming'
            )
        
        meteo = \
            read_radiosonde_ascii(
                input_file = os.path.join(input_folder, input_file), 
                mtime = mtime
                )
        
        meteo = convert_units(data = meteo)
        
        rsonde_start_date, rsonde_start_time = date_from_filename(input_file)
        metadata["rsonde_start_date"] = rsonde_start_date
        metadata["rsonde_start_time"] = rsonde_start_time
        
        meteo = add_number_density(meteo)

        meteo = add_saturation_vapour_pressure(meteo)

        evaluate_flag(file_validity_flag, label = 'wyoming')
        
        metadata['meteo_file'] = 'wyoming'
        metadata['meteo_validity_flag'] = file_validity_flag

    else:
        meteo = None
    
    return(meteo, metadata)


def select_radiosonde_ecmwf(input_folder, signal_times, metadata):

    mtime = get_mid_time(signal_times)
    
    bnames, file_validity_flag = \
        get_filename_list(input_folder, pattern = '*_ecmwf.nc')

    if file_validity_flag != "empty":
        input_file, file_validity_flag = find_nearest_file(
            mtime = mtime, 
            bnames = bnames,
            filetype = 'ecmwf'
            )
        
        data = xr.open_dataset(os.path.join(input_folder, input_file))
        
        metadata['latitude'] = data.latitude.values
        metadata['longitude'] = data.longitude.values
        
        height = data.height.interp(time = mtime)
        pressure = data.pressure.interp(time = mtime)
        temperature = data.temperature.interp(time = mtime)
        relative_humidity = data.rh.interp(time = mtime)
        
        meteo = xr.concat([pressure, temperature, relative_humidity],
                         dim="parameters")
        
        meteo.attrs = {}
        
        meteo = meteo.assign_coords(parameters = ['P', 'T', 'RH'])
        meteo = meteo.rename(level='height')
        meteo = meteo.assign_coords(height = height.values + metadata['station_altitude']).sortby("height")
        
        meteo = add_number_density(meteo)

        meteo = add_saturation_vapour_pressure(meteo)

        metadata["rsonde_start_date"] = meteo.time.dt.strftime("%Y%m%d").values
        metadata["rsonde_start_time"] = meteo.time.dt.strftime("%H%M").values

        evaluate_flag(file_validity_flag, label = 'ecmwf')
        
        metadata['meteo_file'] = 'ecmwf'
        metadata['meteo_validity_flag'] = file_validity_flag

    else:
        meteo = None

    return(meteo, metadata)
    