#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 20:11:12 2022

@author: nick
"""

from ..arc import us_std
import xarray as xr
import os 
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from ..arc.raman_scattering import GaussianFilter, SquareFilter
from ..arc.raman_scattering import RotationalRaman
from ..arc import make_gas, rayleigh_scattering 
from ..arc.utilities import number_density_at_pt

def short_molec(heights, ranges, system_info, channel_info, 
                time_info, external_info):

    print('-----------------------------------------')
    print('Start making molecular calculations...')
    print('-----------------------------------------')
    
    ldr_types = ['p', 'c', 't']
    bsc_stypes = ['r', 't', 'n', 'o', 'h', 'l', 'm', 'x']
    
    if 'Sounding_File_Name' in system_info.keys() :
        rsonde_path = os.path.join(os.path.dirname(external_info['input_file']),  
                                   system_info['Sounding_File_Name'])
        
        meteo, molec_info = from_rsonde(path = rsonde_path, 
                                        heights = heights,
                                        elevation = system_info.Altitude_meter_asl)
        
        print('-- Radiosond file succesfully parsed!')

        
    elif 'Pressure_at_Lidar_Station' in system_info.keys() and \
        'Temperature_at_Lidar_Station' in system_info.keys():
                        
        meteo, molec_info = \
            from_us_std(heights = heights,
                        pressure = system_info.Pressure_at_Lidar_Station, 
                        temperature = system_info.Temperature_at_Lidar_Station,
                        elevation = system_info.Altitude_meter_asl)

        print('-- US Standard Atmosphere succesfully constructed!')
    
    
    channels = ranges.channel.values

    bins = ranges.bins.values
    
    ewl = channel_info.Emitted_Wavelength
    
    dwl = channel_info.Detected_Wavelength
    
    bdw = channel_info.Channel_Bandwidth
    
    ch_type = pd.Series([ch[5] for ch in channel_info.index], index = channels)

    ch_stype = pd.Series([ch[7] for ch in channel_info.index], index = channels)
    
    properties = ['Backscatter_Cross_Section', 
                  'Extinction_Cross_Section_Forward', 
                  'Extinction_Cross_Section_Backward',  
                  'Backscatter_Coefficient', 
                  'Extinction_Coefficient_Forward', 
                  'Extinction_Coefficient_Backward', 
                  'Molecular_Linear_Depolarization_Ratio', 
                  'Attenuated_Backscatter']
    
    molec = xr.DataArray(dims = ['properties', 'channel', 'bins'],
                         coords = [properties, channels, bins])
    
    
    c_N2_default = 0.780796
    c_O2_default = 0.209448
    c_Ar_default = 0.009339
    c_CO2_default = 0.000416
    
    c_total = c_N2_default + c_O2_default + c_Ar_default + c_CO2_default
    
    for ch in channels:
        print(f'-- Processing channel {ch}')

        ch_d = dict(channel = ch)
        
        ewl_ch = ewl.loc[ch]
        
        dwl_ch = dwl.loc[ch]

        bdw_ch = bdw.loc[ch]
        
        mdr = np.nan * np.zeros(bins.size)
        bxs_tot = np.nan * np.zeros(bins.size)
        exs_fth = np.nan * np.zeros(bins.size)
        exs_bck = np.nan * np.zeros(bins.size)
        
        wns_ch = 1E7 * (1. / ewl_ch - 1. / dwl_ch )
        
        if ch_type.loc[ch] in ['p', 'c', 't', 'a'] and np.abs(wns_ch) > 50:
            raise Exception(f"-- Error: Channel {ch_d['channel']} has type '{ch_type.loc[ch]}' (not Raman or fluorescence) but the wavenumber shift between the provided detected {dwl_ch} and emitted {ewl_ch} wavelengths is too large. Please correct the corresponding fields in the configuration file and/or in the raw file header")

        if ch_type.loc[ch] == 'r' and np.abs(wns_ch) > 500:
            raise Exception(f"-- Error: Channel {ch_d['channel']} has type '{ch_type.loc[ch]}' (rotational Ramn) but the wavenumber shift between the provided detected {dwl_ch} and emitted {ewl_ch} wavelengths is too large. Please correct the corresponding fields in the configuration file and/or in the raw file header")

        if ch_type.loc[ch] == 'v' and (wns_ch > 9000 or wns_ch < -50):
            raise Exception(f"-- Error: Channel {ch_d['channel']} has type '{ch_type.loc[ch]}' (vibrational Ramn) but the wavenumber shift between the provided detected {dwl_ch} and emitted {ewl_ch} wavelengths is too large. Please correct the corresponding fields in the configuration file and/or in the raw file header")

        max_bin = meteo.loc[ch_d].dropna(dim = 'bins', how = 'any')\
            .loc[dict(properties  = 'Number_Density')].bins.values[-1]
        sel_bins = np.linspace(bins[0], max_bin, 10).astype(int)
        sl_d = dict(channel = ch, bins = sel_bins)
            
        # Defining a gaussian filter for the IFF
        if ch_type.loc[ch] == 'r':
            filter_trans = SquareFilter(dwl_ch, bdw_ch) 
        elif ch_type.loc[ch] == 'v':
            filter_trans = None
        else:
            filter_trans = GaussianFilter(dwl_ch, bdw_ch)    

        # Extracting molecular parameters in specific range bins
        N = meteo.loc[ch_d].loc[dict(properties  = 'Number_Density')] .values
        P = meteo.loc[sl_d].loc[dict(properties  = 'Pressure')].values
        T = meteo.loc[sl_d].loc[dict(properties  = 'Temperature')].values
        RH = meteo.loc[sl_d].loc[dict(properties  = 'Relative_Humidity')].values

        # Convert Relative humidity to water vapor molar fraction
        ewp = saturation_vapour_pressure(P = P, T = T)
        c_H2O = ewp * (RH / 100.) / P         
        c_H2O[np.isnan(c_H2O)] = 0.

        # Calculate the molar fraction profiles for the other gases taking into account the water vapor profile
        c_N2  = c_N2_default  / (c_total + c_H2O) 
        c_O2  = c_O2_default  / (c_total + c_H2O) 
        c_Ar  = c_Ar_default  / (c_total + c_H2O)
        c_CO2 = c_CO2_default / (c_total + c_H2O)

        mdr_i = np.nan * np.zeros(sel_bins.size)
        bxs_tot_i = np.nan * np.zeros(sel_bins.size)
        exs_fth_i = np.nan * np.zeros(sel_bins.size)
        exs_bck_i = np.nan * np.zeros(sel_bins.size)

        # Automatically calculate the Raman vibrational shift from the emitted wavelength
        if ch_type.loc[ch] not in ['r', 'v', 'f']:
            swl_ch = ewl_ch
        elif ch_type.loc[ch] == 'v':
            if ch_stype.loc[ch] == 'n':
                swl_ch = 1./(1./ewl.loc[ch] - 2331.*1E-7) # N2 wavenumber shift is 2331 cm-1
            else:
                print(f'-- Warning: Molecular atmosphere calculations are only supported for N2. No molecular calculations were performed for channel {ch}')
        else:
            swl_ch = dwl_ch    
               
        for i in range(sel_bins.size):
            
            # Create the gas dictonaries needed for the Raman calculations
            N2 = make_gas.N2(ewl_ch,relative_concentration = c_N2[i])
            O2 = make_gas.O2(ewl_ch,relative_concentration = c_O2[i])
            Ar = make_gas.Ar(ewl_ch,relative_concentration = c_Ar[i])
            CO2 = make_gas.CO2(ewl_ch,relative_concentration = c_CO2[i])
            H2O = make_gas.H2O(ewl_ch,relative_concentration = c_H2O[i])

            # Create the RR object for the backscatter cross sections
            rrb = RotationalRaman(wavelength = ewl_ch, 
                                  max_J = 51, 
                                  temperature = T[i], 
                                  optical_filter = filter_trans,
                                  N2_parameters = N2, 
                                  O2_parameters = O2, 
                                  Ar_parameters = Ar, 
                                  CO2_parameters = CO2, 
                                  H2O_parameters = H2O) 

            # Create the RR object for the forward scattering/extinction cross sections            
            rre_f = RotationalRaman(wavelength = ewl_ch, 
                                    max_J = 51, 
                                    temperature = T[i], 
                                    N2_parameters = N2, 
                                    O2_parameters = O2, 
                                    Ar_parameters = Ar, 
                                    CO2_parameters = CO2, 
                                    H2O_parameters = H2O,
                                    istotal = True) 

            # Create the gas dictonaries needed for the Raman calculations
            N2 = make_gas.N2(swl_ch,relative_concentration = c_N2[i])
            O2 = make_gas.O2(swl_ch,relative_concentration = c_O2[i])
            Ar = make_gas.Ar(swl_ch,relative_concentration = c_Ar[i])
            CO2 = make_gas.CO2(swl_ch,relative_concentration = c_CO2[i])
            H2O = make_gas.H2O(swl_ch,relative_concentration = c_H2O[i])
            
            # Create the RR object for the backward scattering/extinction cross sections            
            rre_b = RotationalRaman(wavelength = swl_ch, 
                                    max_J = 51, 
                                    temperature = T[i], 
                                    N2_parameters = N2, 
                                    O2_parameters = O2, 
                                    Ar_parameters = Ar, 
                                    CO2_parameters = CO2, 
                                    H2O_parameters = H2O,
                                    istotal = True) 


            bxs_tot_i[i], _, _ = rrb.rayleigh_cross_section()
            exs_fth_i[i], _, _ = rre_f.rayleigh_cross_section()
            exs_bck_i[i], _, _ = rre_b.rayleigh_cross_section()
            
            if ch_stype.loc[ch] == 'n':
                bxs_tot_i[i] = c_N2[i] * bxs_tot_i[i]
            if ch_stype.loc[ch] == 'o':
                bxs_tot_i[i] = c_O2[i] * bxs_tot_i[i]

            mdr_i[i] = rrb.delta_mol_rayleigh(method = 'line_summation')

        # Interpolate for every range bin (calculations are performed only on selected bins)
        bxs_tot_f = interp1d(sel_bins, bxs_tot_i, bounds_error = False, fill_value = np.nan)
        exs_fth_f = interp1d(sel_bins, exs_fth_i, bounds_error = False, fill_value = np.nan)
        exs_bck_f = interp1d(sel_bins, exs_bck_i, bounds_error = False, fill_value = np.nan)
        mdr_f = interp1d(sel_bins, mdr_i, bounds_error = False, fill_value = np.nan)
        
        if ch_stype.loc[ch] in bsc_stypes:
            bxs_tot = bxs_tot_f(bins)  
        exs_bck = exs_bck_f(bins)
        exs_fth = exs_fth_f(bins)
        
        bcf_tot = N * bxs_tot 
        ecf_bck = N * exs_bck
        ecf_fth = N * exs_fth

        if ch_type.loc[ch] in ldr_types:
            mdr = mdr_f(bins)
            
        rng = ranges.loc[ch_d].values
        trn_bck = transmittance(x_range = rng, extinction = ecf_bck)
        trn_fth = transmittance(x_range = rng, extinction = ecf_fth)
        
        if ch_type.loc[ch] not in ['p', 'c']:
            atb = bcf_tot * trn_fth * trn_bck
        elif ch_type.loc[ch] in ['p']:
            atb = bcf_tot * trn_fth * trn_bck / (1. + mdr) 
        elif ch_type.loc[ch] in ['c']:
            atb = bcf_tot * trn_fth * trn_bck * mdr / (1. + mdr) 
                            
        # Pack into a single xarray object     
        molec.loc[dict(properties = 'Extinction_Cross_Section_Backward')]\
            .loc[ch_d] = exs_bck
        molec.loc[dict(properties = 'Extinction_Cross_Section_Forward')]\
            .loc[ch_d] = exs_fth
        molec.loc[dict(properties = 'Extinction_Coefficient_Backward')]\
            .loc[ch_d] = ecf_bck   
        molec.loc[dict(properties = 'Extinction_Coefficient_Forward')]\
            .loc[ch_d] = ecf_fth   
        molec.loc[dict(properties = 'Molecular_Linear_Depolarization_Ratio')]\
            .loc[ch_d] = mdr
        molec.loc[dict(properties = 'Backscatter_Cross_Section')]\
            .loc[ch_d] = bxs_tot       
        molec.loc[dict(properties = 'Backscatter_Coefficient')]\
            .loc[ch_d] = bcf_tot    
        molec.loc[dict(properties = 'Attenuated_Backscatter')]\
            .loc[ch_d] = atb

        # if ch_type.loc[ch] not in ldr_types:
        #     print(f"-- Warning: No molec. LDR calculation. Channel type {ch_type.loc[ch]} not in {ldr_types}")
        # if ch_stype.loc[ch] not in bsc_stypes:
        #     print(f"-- Warning: No molec. bsc calculation. Channel type {ch_stype.loc[ch]} not in {bsc_stypes}")
        

    print('-- Molecular calculations succesfully complete!')
    print('-----------------------------------------')
    print('')

        
    return(molec, molec_info, meteo)
        
def from_us_std(pressure, temperature, elevation, heights):
    
    channels = heights.channel
    
    bins = heights.bins
    
    properties = ['Pressure', 'Temperature', 
                  'Relative_Humidity', 'Number_Density']
    
    meteo = xr.DataArray(coords = [properties, channels, bins],
                         dims = ['properties', 'channel', 'bins'])
    
    molec_info = pd.Series()
    
    molec_info['Molecular_atmosphere_method'] = 'US_Standard_Atmosphere'
    
    atm = us_std.Atmosphere(t_r = temperature, p_r = pressure, alt = elevation)

    for ch in channels:
    
        ch_d = dict(channel = ch)
        
        P = np.array([atm.pressure(h) for h in heights.loc[ch_d].values + elevation]) 
        
        T = np.array([atm.temperature(h) for h in heights.loc[ch_d].values + elevation]) 
        
        RH = np.zeros(P.size)
        
        N = number_density_at_pt(pressure = P, 
                                 temperature = T, 
                                 relative_humidity = RH, 
                                 ideal = True)
        
        meteo.loc[ch_d].loc[dict(properties = 'Pressure')] = P
        
        meteo.loc[ch_d].loc[dict(properties = 'Temperature')] = T
        
        meteo.loc[ch_d].loc[dict(properties = 'Relative_Humidity')] = RH

        meteo.loc[ch_d].loc[dict(properties = 'Number_Density')] = N
    
    return(meteo, molec_info)

def from_rsonde(path, heights, elevation):

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
    
    molec_info['Molecular_atmosphere_method'] = 'Radiosonde'

    attrs = file.attrs
    
    keys_alt = ['Altitude_meter_asl', 'Latitude_degrees_north', 
                'Longitude_degrees_east']
    
    keys_org = ['Sounding_Station_Name', 
                'Sounding_Start_Date', 'Sounding_Start_Time_UT',
                'WMO_Station_Number', 'WBAN_Station_Number']

    for key in keys_alt:
        if key in attrs.keys():
            molec_info[f'Sounding_{key}'] = attrs[key]

    for key in keys_org:
        if key in attrs.keys():
            molec_info[key] = attrs[key]
    
    properties = ['Pressure', 'Temperature', 
                  'Relative_Humidity', 'Number_Density']
        
    meteo = xr.DataArray(coords = [properties, channels, bins],
                         dims = ['properties', 'channel', 'bins'])
    
    for ch in channels:
        
        height_arr = heights.loc[ch].values + elevation
        
        # Extrapolate near ground if the first radiosonde altitude is above the lidar ground altitude
        if H[0] > height_arr[0]:
            H = np.insert(H, 0, height_arr[0])
            P = np.insert(P, 0, P[0])
            T = np.insert(T, 0, T[0])
            RH = np.insert(RH, 0, RH[0])

        P_f = interp1d(H, P, bounds_error = False, fill_value = np.nan)
        
        T_f = interp1d(H, T, bounds_error = False, fill_value = np.nan)
        
        RH_f = interp1d(H, RH, bounds_error = False, fill_value = np.nan)

        
        P_i = P_f(height_arr)
        
        T_i = T_f(height_arr)
        
        RH_i = RH_f(height_arr)
        
        N_i = number_density_at_pt(pressure = P_i, 
                                   temperature = T_i, 
                                   relative_humidity = RH_i, 
                                   ideal = True)
        
        ch_d = dict(channel = ch)
                
        meteo.loc[ch_d].loc[dict(properties = 'Pressure')] = P_i
        
        meteo.loc[ch_d].loc[dict(properties = 'Temperature')] = T_i
        
        meteo.loc[ch_d].loc[dict(properties = 'Relative_Humidity')] = RH_i
        
        meteo.loc[ch_d].loc[dict(properties = 'Number_Density')] = N_i
    
    return(meteo, molec_info)

def saturation_vapour_pressure(P, T):
    
    ewp = (1.0007 + 3.46 * 1E-6 * P) * 6.1121 * \
        np.exp(17.502 * (T - 273.15)/(240.97 + T - 273.15))
        
    return(ewp)

def transmittance(x_range, extinction):

    #Calculate transmittance
    Transmittance = np.exp(-cumtrapz(extinction, x_range))
    
    # Copy the first transmission values to the 0 range bin
    Transmittance = np.insert(Transmittance, 0, 1)

    return(Transmittance)
    
    
