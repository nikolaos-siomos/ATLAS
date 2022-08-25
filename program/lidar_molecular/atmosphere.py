#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 20:11:12 2022

@author: nick
"""

from molecular import us_std
import xarray as xr
import os 
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from lidar_molecular.raman_scattering import GaussianFilter
from lidar_molecular.raman_scattering import RotationalRamanBackscatter
from lidar_molecular import make_gas, rayleigh_scattering 

def short_molec(heights, meas_info, channel_info, time_info, external_info):
    
    if 'Sounding_File_Name' in meas_info.keys() :
        rsonde_path = os.path.join(os.path.dirname(external_info['input_file']),  
                                   meas_info['Sounding_File_Name'])
        
        meteo, molec_info = from_rsonde(path = rsonde_path, 
                                        heights = heights)
        
    elif 'Pressure_at_Lidar_Station' in meas_info.keys() and \
        'Temperature_at_Lidar_Station' in meas_info.keys():
                        
        meteo, molec_info = \
            from_us_std(heights = heights,
                        pressure = meas_info.Pressure_at_Lidar_Station, 
                        temperature = meas_info.Temperature_at_Lidar_Station,
                        elevation = meas_info.Altitude_meter_asl)
    
    
    channels = heights.channel.values

    bins = heights.bins.values
    
    ewl = channel_info.Emitted_Wavelength
    
    dwl = channel_info.Detected_Wavelength
    
    bdw = channel_info.Channel_Bandwidth
    
    gaussian_iff = GaussianFilter(dwl,bdw)    
    
    c_N2 = 0.780796
    c_O2 = 0.209448
    c_Ar = 0.009339
    c_CO2 = 0.000416
    
    c_td = c_N2 + c_O2 + c_Ar + c_CO2

    def saturation_vapour_pressure(P, T):
        ewp = (1.0007 + 3.46 * 1E-6 * P) * 6.1121 * \
            np.exp(17.502 * (T - 273.15)/(240.97 + T - 273.15))
        return(ewp)

    sel_bins = np.linspace(bins[0], bins[-1], 100).astype(int)

    ldr_flt = np.nan * np.zeros(size = sel_bins.size)

    bxs_flt = np.nan * np.zeros(size = sel_bins.size)
    
    for ch in channels:
        
        sl_d = dict(channel = ch, bins = sel_bins)
        
        P = meteo.loc[sl_d].loc[dict(properties  = 'pressure')].values
        T = meteo.loc[sl_d].loc[dict(properties  = 'temperature')].values
        RH = meteo.loc[sl_d].loc[dict(properties  = 'humidity')] .values
        
        ewp = saturation_vapour_pressure(P = P, T = T)
        
        c_H2O = ewp * (RH / 100.) / P 

        for i in range(sel_bins.size):
            
            c_N2_i = c_N2/(c_td + c_H2O[i]) 
            c_O2_i = c_O2/(c_td + c_H2O[i]) 
            c_Ar_i = c_Ar/(c_td + c_H2O[i])
            c_CO2_i = c_CO2/(c_td + c_H2O[i])
            c_H2O_i = c_H2O[i]
        
            N2 = make_gas.N2(ewl,relative_concentration = c_N2_i)
            O2 = make_gas.O2(ewl,relative_concentration = c_O2_i)
            Ar = make_gas.Ar(ewl,relative_concentration = c_Ar_i)
            CO2 = make_gas.CO2(ewl,relative_concentration = c_CO2_i)
            H2O = make_gas.H2O(ewl,relative_concentration = c_H2O_i)
        
            rrb = RotationalRamanBackscatter(wavelength = ewl, 
                                             max_J = 101, 
                                             temperature = T[i], 
                                             optical_filter = gaussian_iff,
                                             N2_parameters = N2, 
                                             O2_parameters = O2, 
                                             Ar_parameters = Ar, 
                                             CO2_parameters = CO2, 
                                             H2O_parameters = H2O) 
            
            ldr_flt[i] = rrb.delta_mol_rayleigh(method = 'line_summation')

            bxs_flt[i], _, _ = rrb.quantum_mechanic_rayleigh_cross_section()
            
        

            
    
    
    molec = meteo.copy()
    
    return(molec, molec_info)
        
def from_us_std(pressure, temperature, elevation, heights):
    
    channels = heights.channel
    
    bins = heights.bins
    
    properties = ['pressure','temperature','humidity']
    
    meteo = xr.DataArray(coords = [properties, channels, bins],
                         dims = ['properties', 'channel', 'bins'])
    
    molec_info = pd.Series()
    
    molec_info['method'] = 'US_Standard_Atmosphere'
    
    atm = us_std.Atmosphere(t_r = temperature, p_r = pressure, alt = elevation)

    for ch in channels:
    
        ch_d = dict(channel = ch)
        
        P = np.array([atm.pressure(h) for h in heights.loc[ch_d].values]) 
        
        T = np.array([atm.temperature(h) for h in heights.loc[ch_d].values]) 
        
        meteo.loc[ch_d].loc[dict(properties = 'pressure')] = P
        
        meteo.loc[ch_d].loc[dict(properties = 'temperature')] = T
        
        meteo.loc[ch_d].loc[dict(properties = 'humidity')] = np.zeros(P.size)
    
    return(meteo, molec_info)

def from_rsonde(path, heights):

    channels = heights.channel
    
    bins = heights.bins
    
    file = xr.open_dataset(path)
    
    P = file.Pressure.values
    T = file.Temperature.values
    H = file.Altitude.values
    
    if 'RelativeHumidity' in file.variables:
        RH = file.RelativeHumidity
    
    else:
        RH = np.zeros(H.size)
    
    molec_info = pd.Series()
    
    molec_info['method'] = 'Radiosonde'

    attrs = file.attrs
    
    for key in attrs.keys():
        molec_info[key] = attrs[key]
    
    properties = ['pressure','temperature','humidity']
    
    meteo = xr.DataArray(coords = [properties, channels, bins],
                         dims = ['properties', 'channel', 'bins'])
    
    for ch in channels:
        P_f = interp1d(H, P, bounds_error = False, fill_value = np.nan)
        
        T_f = interp1d(H, T, bounds_error = False, fill_value = np.nan)
        
        RH_f = interp1d(H, RH, bounds_error = False, fill_value = np.nan)
        
        ch_d = dict(channel = ch)
        
        height_arr = heights.loc[ch].values
        
        meteo.loc[ch_d].loc[dict(properties = 'pressure')] = P_f(height_arr)
        
        meteo.loc[ch_d].loc[dict(properties = 'temperature')] = T_f(height_arr)
        
        meteo.loc[ch_d].loc[dict(properties = 'humidity')] = RH_f(height_arr)
    
    return(meteo, molec_info)
    
    