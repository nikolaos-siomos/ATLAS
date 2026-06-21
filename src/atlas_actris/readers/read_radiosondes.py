#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 11:17:36 2025

@author: nikos
"""

import re
import numpy as np
import xarray as xr
import os, sys, glob
import netCDF4 as nc

from utils.error_classes import CustomWarning
from utils.printouts import print_header, print_entry
from molecular.utilities import number_density_at_pt, saturation_vapour_pressure

from utils.unit_conversions import (
    km_asl_to_m_asl, 
    m_agl_to_m_asl, 
    km_agl_to_m_asl, 
    hPa_to_Pa, 
    atm_to_Pa, 
    C_to_K, 
    Cx10_to_K, 
    percent_to_fraction
    )

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
    
def convert_units(data, units = None, ground = None):
        
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

    make_nc_var(ds, name = 'Pressure', value = meteo.loc[dict(atmo_parameters = 'P')].values, dtype = 'float', dims = ('points',))

    make_nc_var(ds, name = 'Temperature', value = meteo.loc[dict(atmo_parameters = 'T')].values, dtype = 'float', dims = ('points',))

    if 'RH' in meteo.atmo_parameters.values:
        make_nc_var(ds, name = 'RelativeHumidity', value = meteo.loc[dict(atmo_parameters = 'RH')].values, dtype = 'float', dims = ('points',))

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
            
def clean_meteo_height_index(meteo):
    """Sort radiosonde profiles and remove invalid or duplicate height levels.

    xarray interpolation requires the interpolation coordinate to be unique.
    Radiosonde files can contain repeated height levels, especially after
    rounding or unit conversion. Keep the first occurrence of each valid
    height level and preserve the original atmospheric-parameter dimension.
    """

    if "height_asl" not in meteo.dims:
        return meteo

    meteo = meteo.sortby("height_asl")

    height = meteo.height_asl.values
    valid = np.isfinite(height)

    if not np.all(valid):
        CustomWarning(
            "Invalid radiosonde height levels detected. "
            "Removing NaN or infinite height values."
        )
        meteo = meteo.isel(height_asl=valid)
        height = meteo.height_asl.values

    if height.size == 0:
        raise Exception(
            "-- Error: No valid radiosonde height levels remain after cleaning."
        )

    _, unique_ind = np.unique(height, return_index=True)
    unique_ind = np.sort(unique_ind)

    if unique_ind.size < height.size:
        CustomWarning(
            "Duplicate radiosonde height levels detected. "
            "Keeping the first occurrence of each height level."
        )
        meteo = meteo.isel(height_asl=unique_ind)

    return meteo

def add_number_density(meteo):
    
    height_asl = meteo.height_asl.values
    pressure = meteo.loc['P',:].values
    temperature = meteo.loc['T',:].values
    relative_humidity = meteo.loc['RH',:].fillna(0.).values
        
    number_density = number_density_at_pt(
        pressure = pressure, 
        temperature = temperature, 
        relative_humidity = relative_humidity, 
        ideal=False
        )
    
    number_density = xr.DataArray([number_density], 
                                  coords = [['N'], height_asl], 
                                  dims = ['atmo_parameters', 'height_asl'])
    
    meteo = xr.concat([meteo, number_density], dim="atmo_parameters")
    
    return(meteo)

def add_saturation_vapour_pressure(meteo):
    
    height_asl = meteo.height_asl.values
    pressure = meteo.loc['P',:].values
    temperature = meteo.loc['T',:].values
    
    sat_vap_pressure = saturation_vapour_pressure(
        pressure = pressure, 
        temperature = temperature, 
        )
    
    sat_vap_pressure = xr.DataArray([sat_vap_pressure], 
                                  coords = [['e_s'], height_asl], 
                                  dims = ['atmo_parameters', 'height_asl'])
    
    meteo = xr.concat([meteo, sat_vap_pressure], dim="atmo_parameters")
    
    return(meteo)

def read_radiosonde_ecmwf(station_altitude, radiosonde_info):
    
    radiosonde_file = radiosonde_info.loc['radiosonde_file'].item()
    mtime = radiosonde_info.loc['measurement_time'].item()
    
    data = xr.open_dataset(radiosonde_file)

    # Convert to numpy datetime64 for robust comparison
    mtime64 = np.datetime64(mtime)
    tmin = data.time.values.min()
    tmax = data.time.values.max()

    # If measurement time is outside the ECMWF time window,
    # use the nearest available boundary profile.
    # Otherwise, interpolate normally in time.
    if mtime64 < tmin:
        CustomWarning(
            "Measurement time is before the first ECMWF profile. "
            f"Using first ECMWF profile instead. "
            f"measurement_time={mtime64}, first_ecmwf_time={tmin}"
        )

        height = data.height.sel(time=tmin)
        pressure = data.pressure.sel(time=tmin)
        temperature = data.temperature.sel(time=tmin)
        relative_humidity = data.rh.sel(time=tmin)

    elif mtime64 > tmax:
        CustomWarning(
            "Measurement time is after the last ECMWF profile. "
            f"Using last ECMWF profile instead. "
            f"measurement_time={mtime64}, last_ecmwf_time={tmax}"
        )

        height = data.height.sel(time=tmax)
        pressure = data.pressure.sel(time=tmax)
        temperature = data.temperature.sel(time=tmax)
        relative_humidity = data.rh.sel(time=tmax)

    else:
        height = data.height.interp(time=mtime)
        pressure = data.pressure.interp(time=mtime)
        temperature = data.temperature.interp(time=mtime)
        relative_humidity = data.rh.interp(time=mtime)
    
    meteo = xr.concat(
        [pressure, temperature, relative_humidity],
        dim="atmo_parameters"
    )
    
    meteo.attrs = {}
    
    meteo = meteo.assign_coords(atmo_parameters=['P', 'T', 'RH'])
    meteo = meteo.rename(level='height_asl')
    meteo = meteo.assign_coords(
        height_asl=height.values + station_altitude
    )

    meteo = clean_meteo_height_index(meteo)
    
    meteo = add_number_density(meteo)

    meteo = add_saturation_vapour_pressure(meteo)

    return meteo

def read_radiosonde_wyoming(radiosonde_info):

    radiosonde_file = radiosonde_info.loc['radiosonde_file'].item()
    
    data = np.genfromtxt(radiosonde_file,
                         skip_header = 1, 
                         skip_footer = 0,
                         delimiter = ',', 
                         autostrip = True,
                         usecols = np.array([4, 3, 5, 8]), dtype = float)
    
    atmo_parameters = ['P', 'T', 'RH']         
        
    height = data[:,0]
    
    meteo = xr.DataArray(data[:,1:].T, 
                        coords = [atmo_parameters, height], 
                        dims = ['atmo_parameters', 'height_asl'])
    
        
    meteo = convert_units(meteo, units = ['m_asl','hPa','C','percent'])

    meteo = clean_meteo_height_index(meteo)
    
    meteo = add_number_density(meteo)

    meteo = add_saturation_vapour_pressure(meteo)
    
    return meteo

def read_radiosonde_scc(radiosonde_info):
    
    radiosonde_file = radiosonde_info.loc['radiosonde_file']
    
    file = xr.open_dataset(radiosonde_file)
    
    P = file.Pressure.values
    T = file.Temperature.values
    height = file.Altitude.values
        
    if 'RelativeHumidity' in file.variables:
        RH = file.RelativeHumidity
    
    else:
        RH = np.nan * np.zeros(height.size)
        
    atmo_parameters = ['P', 'T', 'RH']         

    meteo = xr.DataArray(
        np.stack([P, T, RH], axis = 0), 
        coords = [atmo_parameters, height], 
        dims = ['atmo_parameters', 'height_asl']
        )
    
    meteo = convert_units(meteo, units = ['m_asl','hPa','K','percent'])

    meteo = clean_meteo_height_index(meteo)
    
    meteo = add_number_density(meteo)

    meteo = add_saturation_vapour_pressure(meteo)
    
    # P_i = logarithmic_pressure_interpolation(
    #     H = H,
    #     P = P,
    #     heights = height_arr
    #     )
    
    return meteo

def read_radiosonde_ascii(caller_info, radiosonde_info):
    
    radiosonde_file = radiosonde_info.loc['radiosonde_file'].item()
    
    # Unpack parsing options
    rsonde_skip_header = caller_info['rsonde_skip_header']
    rsonde_skip_footer = caller_info['rsonde_skip_footer']
    rsonde_delimiter = caller_info['rsonde_delimiter']
    rsonde_column_index = caller_info['rsonde_column_index']
    rsonde_column_units = caller_info['rsonde_column_units']
    rsonde_station_altitude = caller_info['rsonde_station_altitude']
    
    if rsonde_delimiter == 'S': rsonde_delimiter = None
    if rsonde_delimiter == 'C': rsonde_delimiter = ','
    
    data = np.genfromtxt(radiosonde_file,
                         skip_header = rsonde_skip_header, 
                         skip_footer = rsonde_skip_footer,
                         delimiter = rsonde_delimiter, 
                         autostrip = True,
                         usecols = np.array(rsonde_column_index) - 1, dtype = float)
    
    if len(rsonde_column_units) == 3:
        data = np.vstack([data, np.nan * np.zeros(data.shape[0])])
        
    atmo_parameters = ['P', 'T', 'RH']         
        
    height = data[:,0]
    
    meteo = xr.DataArray(data[:,1:].T, 
                        coords = [atmo_parameters, height], 
                        dims = ['atmo_parameters', 'height_asl'])
    
    meteo = convert_units(
        meteo, 
        units = rsonde_column_units, 
        ground = rsonde_station_altitude
        )

    meteo = clean_meteo_height_index(meteo)
    
    meteo = add_number_density(meteo)

    meteo = add_saturation_vapour_pressure(meteo)
    
    return meteo

def load_radiosonde(caller_info, metadata):
    
    meteo = {}
                
    print_header("Parsing radiosonde file")

    radiosonde_info = metadata['radiosonde_info']
                
    for key in radiosonde_info:
        
        radiosonde_status = radiosonde_info[key].loc['radiosonde_status'].item()
        
        if radiosonde_status == 0:
                    
            radiosonde_format = radiosonde_info[key].loc['radiosonde_format'].item()
            
            if radiosonde_format == 'ecmwf':
                print_entry("Parsing downloaded Cloudnet meteorological files")
                station_altitude = float(metadata['system_info'][key].loc['station_altitude'].values)
                da = read_radiosonde_ecmwf(station_altitude, radiosonde_info[key])
    
            elif radiosonde_format == 'wyoming':
                print_entry("Parsing downloaded Wyoming radiosonde")
                da = read_radiosonde_wyoming(radiosonde_info[key])
                
            elif radiosonde_format == 'ascii':
                print_entry("Parsing manually provided ASCII radiosonde")
                da = read_radiosonde_ascii(caller_info, radiosonde_info[key])
                
            elif radiosonde_format == 'scc':
                print_entry("Parsing manually provided SCC radiosonde")
                da = read_radiosonde_scc(radiosonde_info[key])
            
            else:
                da = None
                CustomWarning(f"Radiosonde format {radiosonde_format} not understood. Radiosonde was not parsed")
            
            if da is not None:
                meteo[key] = da.chunk({
                    "atmo_parameters": -1,
                    "height_asl": -1,
                })
            
    return(meteo, metadata)

    
    