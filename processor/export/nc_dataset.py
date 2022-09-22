#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 12:50:13 2022

@author: nick
"""

import netCDF4
import os
import numpy as np
import xarray as xr

def rayleigh(sig, molec, meteo, meas_info, channel_info, time_info, 
             molec_info,  heights, ranges, version, dir_out):

    print('-----------------------------------------')
    print('Start exporting to a Rayleigh ATLAS file...')
    print('-----------------------------------------')
    
    nc_file = sig[dict(time = 0)]\
        .to_dataset(name = 'Range_Corrected_Signals')\
            .fillna(netCDF4.default_fillvals['f8'])
        
    for idx in meas_info.index.values:
        nc_file.attrs[idx] = meas_info[idx]
    
    for clm in channel_info.columns.values:
        nc_file[clm] = (['channel'], channel_info.loc[:,clm]\
                        .fillna(netCDF4.default_fillvals['f8']))
    
    for clm in time_info.columns.values:
        nc_file[clm] = time_info.loc[:,clm].iloc[0]
    
    for prop in meteo.properties.values:
        nc_file[prop] = (['channel', 'bins'], 
                         meteo.loc[dict(properties = prop)].values)

    for prop in molec.properties.values:
        nc_file[prop] = (['channel', 'bins'], 
                         molec.loc[dict(properties = prop)].values)

    nc_file.attrs['Molecular_atmosphere_method'] = \
        molec_info['Molecular_atmosphere_method']
        
    if molec_info['Molecular_atmosphere_method'] == 'Radiosonde':
        discard = ['Molecular_atmosphere_method', 'Measurement_type']
        keys = [key for key in molec_info.index.values if key not in discard]
        
        for idx in keys:
            nc_file.attrs[f'Radiosonde_{idx}'] = molec_info[idx]
    
            
    nc_file['Height_levels'] = heights   

    nc_file['Range_levels'] = ranges
    
    nc_file.attrs['version'] = version

    nc_file.attrs['processing_software'] = 'ATLAS'
    
    fname = f'ray_{meas_info.Measurement_ID}_ATLAS_{version}.nc'
    
    nc_file.to_netcdf(path = os.path.join(dir_out, fname))
    
    nc_file.close()
    
    print('-- Rayleigh ATLAS file succesfully created!')
    print('-----------------------------------------')
    print('')

    return(nc_file)

def telecover(sig, meas_info, channel_info, time_info, 
              heights, ranges, version, dir_out):

    print('-----------------------------------------')
    print('Start exporting to a Telecover ATLAS file...')
    print('-----------------------------------------')
    
    sector = time_info.sector.values.astype(int)
    
    time = sig.time.values
    
    if any(sector == 1) == True and any(sector == 5) == False:

        idx_north = np.where(sector == 1)[0]
        idx_east  = np.where(sector == 2)[0]
        idx_south = np.where(sector == 3)[0]
        idx_west  = np.where(sector == 4)[0]

        nc_file = xr.Dataset(coords = {'time_n' : time[idx_north],
                                       'time_e' : time[idx_east],
                                       'time_s' : time[idx_south],
                                       'time_w' : time[idx_west],
                                       'channel' : sig.channel.values, 
                                       'bins': sig.bins.values})

        nc_file['Range_Corrected_Signals_North_Sector'] = \
            (['time_n','channel','bins'], sig[dict(time = idx_north)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
                
        nc_file['Range_Corrected_Signals_East_Sector'] = \
            (['time_e','channel','bins'], sig[dict(time = idx_east)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
                
        nc_file['Range_Corrected_Signals_South_Sector'] = \
            (['time_s','channel','bins'], sig[dict(time = idx_south)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
    
        nc_file['Range_Corrected_Signals_West_Sector'] = \
            (['time_w','channel','bins'], sig[dict(time = idx_west)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
        
        for clm in time_info.columns.values:
            if clm not in ['sector']:
                nc_file[f'{clm}_North_Sector'] = \
                    (['time_n'], time_info.loc[:,clm].iloc[idx_north].values)
                    
                nc_file[f'{clm}_East_Sector'] = \
                    (['time_e'], time_info.loc[:,clm].iloc[idx_east].values)
                    
                nc_file[f'{clm}_South_Sector'] = \
                    (['time_s'], time_info.loc[:,clm].iloc[idx_south].values)
                    
                nc_file[f'{clm}_West_Sector'] = \
                    (['time_w'], time_info.loc[:,clm].iloc[idx_west].values)
    

    if any(sector == 1) == False and any(sector == 5) == True:
        
        idx_inner = np.where(sector == 5)[0]
        idx_outer = np.where(sector == 6)[0]
        
        nc_file = xr.Dataset(coords = {'time_i' : time[idx_inner],
                                       'time_o' : time[idx_outer],
                                       'channel' : sig.channel.values, 
                                       'bins': sig.bins.values})
        

        nc_file['Range_Corrected_Signals_Inner_Ring'] = \
            (['time','channel','bins'], sig[dict(time = idx_inner)]\
             .fillna(netCDF4.default_fillvals['f8']).values)  

        nc_file['Range_Corrected_Signals_Outer_Ring'] = \
            (['time','channel','bins'], sig[dict(time = idx_outer)]\
             .fillna(netCDF4.default_fillvals['f8']).values)     

        for clm in time_info.columns.values:
            if clm not in ['sector']:
                nc_file[f'{clm}_Inner_Ring'] = \
                    (['time_i'], time_info.loc[:,clm].iloc[idx_inner].values)
                    
                nc_file[f'{clm}_Outer_Ring'] = \
                    (['time_o'], time_info.loc[:,clm].iloc[idx_outer].values)
        

    if any(sector == 1) == True and any(sector == 5) == True:

        idx_north = np.where(sector == 1)[0]
        idx_east  = np.where(sector == 2)[0]
        idx_south = np.where(sector == 3)[0]
        idx_west  = np.where(sector == 4)[0]
        idx_outer  = np.where(sector == 5)[0]
        idx_inner  = np.where(sector == 6)[0]

        nc_file = xr.Dataset(coords = {'time_n' : time[idx_north],
                                       'time_e' : time[idx_east],
                                       'time_s' : time[idx_south],
                                       'time_w' : time[idx_west],
                                       'time_i' : time[idx_inner],
                                       'time_o' : time[idx_outer],
                                       'channel' : sig.channel.values, 
                                       'bins': sig.bins.values})

        nc_file['Range_Corrected_Signals_North_Sector'] = \
            (['time_n','channel','bins'], sig[dict(time = idx_north)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
                
        nc_file['Range_Corrected_Signals_East_Sector'] = \
            (['time_e','channel','bins'], sig[dict(time = idx_east)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
                
        nc_file['Range_Corrected_Signals_South_Sector'] = \
            (['time_s','channel','bins'], sig[dict(time = idx_south)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
    
        nc_file['Range_Corrected_Signals_West_Sector'] = \
            (['time_w','channel','bins'], sig[dict(time = idx_west)]\
             .fillna(netCDF4.default_fillvals['f8']).values)
                
        nc_file['Range_Corrected_Signals_Inner_Ring'] = \
            (['time_i','channel','bins'], sig[dict(time = idx_inner)]\
             .fillna(netCDF4.default_fillvals['f8']).values)  

        nc_file['Range_Corrected_Signals_Outer_Ring'] = \
            (['time_o','channel','bins'], sig[dict(time = idx_outer)]\
             .fillna(netCDF4.default_fillvals['f8']).values)     

        
        for clm in time_info.columns.values:
            if clm not in ['sector']:
                nc_file[f'{clm}_North_Sector'] = \
                    (['time_n'], time_info.loc[:,clm].iloc[idx_north].values)
                    
                nc_file[f'{clm}_East_Sector'] = \
                    (['time_e'], time_info.loc[:,clm].iloc[idx_east].values)
                    
                nc_file[f'{clm}_South_Sector'] = \
                    (['time_s'], time_info.loc[:,clm].iloc[idx_south].values)
                    
                nc_file[f'{clm}_West_Sector'] = \
                    (['time_w'], time_info.loc[:,clm].iloc[idx_west].values)
    
                nc_file[f'{clm}_Inner_Ring'] = \
                    (['time_i'], time_info.loc[:,clm].iloc[idx_inner].values)
                    
                nc_file[f'{clm}_Outer_Ring'] = \
                    (['time_o'], time_info.loc[:,clm].iloc[idx_outer].values)
        
    for idx in meas_info.index.values:
        nc_file.attrs[idx] = meas_info[idx]
    
    for clm in channel_info.columns.values:
        nc_file[clm] = (['channel'], channel_info.loc[:,clm]\
                        .fillna(netCDF4.default_fillvals['f8']))

            
    nc_file['Height_levels'] = heights   

    nc_file['Range_levels'] = ranges
    
    nc_file.attrs['version'] = version

    nc_file.attrs['processing_software'] = 'ATLAS'
    
    fname = f'tlc_{meas_info.Measurement_ID}_ATLAS_{version}.nc'
    
    nc_file.to_netcdf(path = os.path.join(dir_out, fname))
    
    nc_file.close()

    print('-- Telecover ATLAS file succesfully created!')
    print('-----------------------------------------')
    print('')

    return(nc_file)

def calibration(sig, sig_ray, meteo, molec, meas_info, meas_info_ray, 
                channel_info, channel_info_ray, time_info, time_info_ray, 
                molec_info, heights, ranges, ranges_ray, heights_ray, 
                version, dir_out):

    print('-----------------------------------------')
    print('Start exporting to a Calibration ATLAS file...')
    print('-----------------------------------------')
    
    pos = time_info.calibrator_position.astype(int)
    time = sig.time.values
    
    idx_m45 = np.where(pos == 1)[0]
    idx_p45 = np.where(pos == 2)[0]
    
    nc_file = xr.Dataset(coords = {'time_m45' : time[idx_m45],
                                   'time_p45' : time[idx_p45],
                                   'channel' : sig.channel.values, 
                                   'bins': sig.bins.values})

    nc_file['Range_Corrected_Signals_minus_45'] = \
        (['time_m45','channel','bins'], sig[dict(time = idx_m45)]\
         .fillna(netCDF4.default_fillvals['f8']).values)
            
    nc_file['Range_Corrected_Signals_plus_45'] = \
        (['time_p45','channel','bins'], sig[dict(time = idx_p45)]\
         .fillna(netCDF4.default_fillvals['f8']).values)
        
    for idx in meas_info.index.values:
        nc_file.attrs[f'{idx}_Calibration'] = meas_info[idx]
    
    for clm in channel_info.columns.values:
        nc_file[f'{clm}_Calibration'] = (['channel'], channel_info.loc[:,clm]\
                                         .fillna(netCDF4.default_fillvals['f8']))
    
    for clm in time_info.columns.values:
        if clm not in ['calibrator_position']:
            nc_file[f'{clm}_minus_45'] = (['time_m45'], time_info.loc[:,clm]\
                                          .iloc[idx_m45])
                            
            nc_file[f'{clm}_plus_45'] = (['time_p45'], time_info.loc[:,clm]\
                                         .iloc[idx_p45])
            
    nc_file['Height_levels_Calibration'] = heights   

    nc_file['Range_levels_Calibration'] = ranges
    
    nc_file.attrs['bins_r'] = sig_ray.bins.values
    
    nc_file['Range_Corrected_Signals_Rayleigh'] = (['channel','bins_r'], 
                                                   sig_ray[dict(time = 0)].values)

    for idx in meas_info_ray.index.values:
        nc_file.attrs[f'{idx}_Rayleigh'] = meas_info_ray[idx]
    
    for clm in channel_info_ray.columns.values:
        nc_file[f'{clm}_Rayleigh'] = (['channel'], channel_info_ray.loc[:,clm]\
                                      .fillna(netCDF4.default_fillvals['f8']))
    
    for clm in time_info_ray.columns.values:
        nc_file[f'{clm}_Rayleigh'] = time_info_ray.loc[:,clm].iloc[0]
            
    nc_file['Height_levels_Rayleigh'] = (['channel','bins_r'], heights.values)

    nc_file['Range_levels_Rayleigh'] = (['channel','bins_r'], ranges.values)
    
    for prop in meteo.properties.values:
        nc_file[prop] = (['channel', 'bins'], 
                         meteo.loc[dict(properties = prop)].values)

    for prop in molec.properties.values:
        nc_file[prop] = (['channel', 'bins'], 
                         molec.loc[dict(properties = prop)].values)
    
    nc_file.attrs['Molecular_atmosphere_method'] = \
        molec_info['Molecular_atmosphere_method']
        
    if molec_info['Molecular_atmosphere_method'] == 'Radiosonde':
        discard = ['Molecular_atmosphere_method', 'Measurement_type']
        keys = [key for key in molec_info.index.values if key not in discard]
        
        for idx in keys:
            nc_file.attrs[f'Radiosonde_{idx}'] = molec_info[idx]
    

    nc_file.attrs['version'] = version

    nc_file.attrs['processing_software'] = 'ATLAS'
    
    fname = f'pcl_{meas_info.Measurement_ID}_ATLAS_{version}.nc'
    
    nc_file.to_netcdf(path = os.path.join(dir_out, fname))
    
    nc_file.close()
    
    print('-- Calibration ATLAS file succesfully created!')
    print('-----------------------------------------')
    print('')

    return(nc_file)

def quicklook(sig, meas_info, channel_info, time_info, 
              heights, ranges, version, dir_out):

    print('-----------------------------------------')
    print('Start exporting to a Quicklook ATLAS file...')
    print('-----------------------------------------')
    
    nc_file = sig.to_dataset(name = 'Range_Corrected_Signals')\
        .fillna(netCDF4.default_fillvals['f8'])
        
    for idx in meas_info.index.values:
        nc_file.attrs[idx] = meas_info[idx]
    
    for clm in channel_info.columns.values:
        nc_file[clm] = (['channel'], channel_info.loc[:,clm]\
                        .fillna(netCDF4.default_fillvals['f8']))
    
    for clm in time_info.columns.values:
        nc_file[clm] = (['time'], time_info.loc[:,clm]\
                        .fillna(netCDF4.default_fillvals['f8']))
            
    nc_file['Height_levels'] = heights   

    nc_file['Range_levels'] = ranges
    
    nc_file.attrs['version'] = version

    nc_file.attrs['processing_software'] = 'ATLAS'
    
    fname = f'qck_{meas_info.Measurement_ID}_ATLAS_{version}.nc'
    
    nc_file.to_netcdf(path = os.path.join(dir_out, fname))
    
    print('-- Quicklook ATLAS file succesfully created!')
    print('-----------------------------------------')
    print('')

    return(nc_file)
