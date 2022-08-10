"""
@author: N. Siomos & P. Paschou
"""
import numpy as np
from scipy.integrate import cumtrapz
import xarray as xr
from readers import read_atm
import os
from molecular import height_scales, us_std
import pandas as pd

def elastic(x_range, P, T, info):
    
    ids = info.molec_list.values

    cs_e = info.cs_e.values.astype(float)
    bs = info.bs.values.astype(float)
    
    a_mol = profiles(x_range, P, T, cs_e, ids)
    b_mol = profiles(x_range, P, T, bs, ids)
    
    T_atm = trans(x_range, a_mol)

    #Calculate attenuated backscatter    
    attn = b_mol*np.power(T_atm, 2)

    attn = xr.DataArray(attn,
                        coords=[ids, x_range],
                        dims=['id', 'range'])
    
    info_el = info.copy().set_index('molec_list')
    info_el.loc[:,'cs_r'] = cs_e
    info_el.loc[:,'wave_r'] = info_el.loc[:,'wave_e'].values
    
    return(a_mol, b_mol, attn, info_el)

def depol(x_range, P, T, info):
    
    ids = list(info.molec_list.values)

    ids_p = [id_join + '_P' for id_join in ids] 
    ids_s = [id_join + '_S' for id_join in ids] 
    
    d_mol = info.d_mol.values.astype(float)
    
    cs_e = info.cs_e.values.astype(float)
    bs_p = info.bs.values.astype(float)/(1. + d_mol)
    bs_s = d_mol*info.bs.values.astype(float)/(1. + d_mol)
    
    a_mol_p = profiles(x_range, P, T, cs_e, ids_p)
    a_mol_s = profiles(x_range, P, T, cs_e, ids_s)
    
    b_mol_p = profiles(x_range, P, T, bs_p, ids_p)
    b_mol_s = profiles(x_range, P, T, bs_s, ids_s)
    
    T_atm_p = trans(x_range, a_mol_p)
    T_atm_s = trans(x_range, a_mol_s)

    #Calculate attenuated backscatter    
    attn_p = b_mol_p*np.power(T_atm_p, 2)
    attn_s = b_mol_s*np.power(T_atm_s, 2)

    attn_p = xr.DataArray(attn_p,
                        coords=[ids_p, x_range],
                        dims=['id', 'range'])

    attn_s = xr.DataArray(attn_s,
                        coords=[ids_s, x_range],
                        dims=['id', 'range'])
    
    a_mol = xr.concat([a_mol_p, a_mol_s], dim = 'id')
    b_mol = xr.concat([b_mol_p, b_mol_s], dim = 'id')
    attn = xr.concat([attn_p, attn_s], dim = 'id')
    
    info_dp = info.copy().set_index('molec_list')
    info_dp.loc[:,'cs_r'] = cs_e
    info_dp.loc[:,'wave_r'] = info_dp.loc[:,'wave_e'].values
    
    info_dp_p = pd.DataFrame(data = info_dp.copy().values,
                             columns = info_dp.columns,
                             index = ids_p)
    info_dp_s = pd.DataFrame(data = info_dp.copy().values,
                             columns = info_dp.columns,
                             index = ids_s)    
    
    info_dp_p.loc[:,'bs'] = bs_p
    info_dp_s.loc[:,'bs'] = bs_s
    
    info_dp = pd.concat([info_dp, info_dp_p, info_dp_s])
    
    return(a_mol, b_mol, attn, info_dp)

def raman(x_range, P, T, info):
    
    ids = info.molec_list.values

    cs_e = info.cs_e.values.astype(float)
    cs_r = info.cs_r.values.astype(float)
    bs = info.bs.values.astype(float)
    
    a_mol_e = profiles(x_range, P, T, cs_e, ids)
    a_mol_r = profiles(x_range, P, T, cs_r, ids)
    b_mol = profiles(x_range, P, T, bs, ids)
    
    T_atm_e = trans(x_range, a_mol_e)
    T_atm_r = trans(x_range, a_mol_r)

    #Calculate attenuated backscatter    
    attn = b_mol*T_atm_e*T_atm_r

    attn = xr.DataArray(attn,
                        coords=[ids, x_range],
                        dims=['id', 'range'])
    
    info_rm = info.copy().set_index('molec_list')
    
    return(a_mol_e, a_mol_r, b_mol, attn, info_rm)

def profiles(x_range, P, T, xs, ids):
    
    prof = np.nan*np.zeros((len(ids), x_range.size))

    for j in range(len(ids)):
        prof[j,:] = xs[j]*P[:]/T[:]
        
    prof = xr.DataArray(prof,
                        coords=[ids, x_range],
                        dims=['id', 'range'])
    
    return(prof)

def trans(x_range, a_mol):
    
    #Calculate transmittance
    Transmittance = np.exp(-cumtrapz(a_mol, x_range))
    
    # Copy the first transmission values to the 0 range bin
    Transmittance = np.insert(Transmittance, 0, 1, axis = 1)

    return(Transmittance)

def atmospheric_profiles(dir_inp, imol, alt_lr, gnd_us, gnd_alt, SZA):

    if imol == 3:
        alt_rs, P, T, RH, imol = read_atm.rsonde(os.path.join(dir_inp, 'rsonde'))

    if imol == 2:
        alt_rs, P, T, RH, imol = read_atm.model(os.path.join(dir_inp, 'model'))

    if imol == 1:
        #check if T and P are given from the user
        gnd_us.temp, gnd_us.pres = us_NTP(gnd_us.temp, gnd_us.pres)
        
        #Create the temperature and pressure profiles from US standard atmosphere
        atmosphere=us_std.Atmosphere(gnd_us.temp, gnd_us.pres, gnd_alt)
        alt_bins = alt_lr.shape[0]
        P = np.zeros(alt_bins) # Pressure in hPa
        T = np.zeros(alt_bins) # Temperature in K
        RH = np.zeros(alt_bins) # Relative Humidity, set to 0 for us std 
        alt_rs = alt_lr    
        for i in range(0, alt_bins):
            P[i]=atmosphere.pressure(alt_lr[i])
            T[i]=atmosphere.temperature(alt_lr[i])

    # alt_rs, P, T = height_scales.profiles_interp(alt_rs, P, T, alt_lr)

    alt_rs, P, T, RH = height_scales.profiles_interp(alt_rs, P, T, RH, alt_lr)

    #number_density
    N = number_density_at_pt(P, T, RH)

    x_H2O = molar_fraction_water_vapour(P, T, RH)
    
    N_H2O = N*x_H2O
    
    return(alt_rs, P, T, RH, N, N_H2O)


def us_NTP(temp, pres):
    '''
    If the temperature and pressure are not given from the user, the Normal
    Temperature and Pressure state (NTP) are used instead in the U.S. Standard Atmosphere 
    model

    Parameters
    ----------
    temp : FLOAT
        The ground temperature value (in Kelvin) for the U.S. Standard Atmosphere model.
    pres : FLOAT
        The ground pressure value (in hPa) for the U.S. Standard Atmosphere model.

    Returns
    -------
    temp : FLOAT
        The ground temperature value (in Kelvin) for the U.S. Standard Atmosphere model.
    pres : FLOAT
        The ground pressure value (in hPa) for the U.S. Standard Atmosphere model.
    '''
    if temp == '' or pres == '':
        
        if temp == '':
            temp = 293.15
            
        if pres == '':    
            pres = 1013.25

        print('---- Warning! No temperature and/or pressure values at ground found for U.S. '+\
              'standard atmosphere --> Using 293.15 K and/or 1013.25 hPa instead..')
            
    return(temp, pres)

    """
    The following functions are obtained from:
    https://gitlab.com/ioannis_binietoglou/lidar_molecular
    
    ################################################################################ 
    # The MIT License (MIT)
    #
    # Copyright (c) 2015, Ioannis Binietoglou
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in all
    # copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    # SOFTWARE.
    ################################################################################    
    """


def number_density_at_pt(pressure, temperature, relative_humidity, ideal=False):
    """ Calculate the number density for a given temperature and pressure,
    taking into account the compressibility of air.
    
    Parameters
    ----------
    pressure: float or array
       Pressure in hPa
    temperature: float or array
       Temperature in K
    relative_humidity: float or array (?)
       ? The relative humidity of air (Check)
    ideal: boolean
       If False, the compressibility of air is considered. If True, the 
       compressibility is set to 1.
    
    Returns
    -------
    n: array or array
       Number density of the atmosphere [m-3]   
    """
    Xw = molar_fraction_water_vapour(pressure, temperature, relative_humidity)
        
    if ideal:
        Z = 1
    else:    
        Z = compressibility_of_moist_air(pressure, temperature, Xw)

    p_pa = pressure * 100.  # Pressure in pascal

    k_b = 1.38064852 * 10**-23  # Boltzmann constant in J/K

    n = p_pa / (Z * temperature * k_b)
    
    return(n)

def compressibility_of_moist_air(pressure, temperature, molar_fraction):
    """ Compressibility of moist air.
    
    Parameters
    ----------
    pressure: float
       Atmospheric pressure [hPa]
    temperature: float
       Atmospehric temperature [K]   
    molar_fraction: float
       Molar fraction.
       
    Note
    ----
    Eg. 16 of Tomasi et al. is missing a bracket. The formula of Ciddor 1996
    was used instead.
    
    """
    a0 = 1.58123e-6  # K Pa-1
    a1 = -2.9331e-8  # Pa-1
    a2 = 1.1043e-10  # K Pa-1
    b0 = 5.707e-6  # K Pa-1
    b1 = -2.051e-8  # Pa-1
    c0 = 1.9898e-4  # Pa-1
    c1 = -2.376e-6  # Pa-1
    d0 = 1.83e-11  # K2 Pa-2
    d1 = -7.65e-9  # K2 Pa-2

    p = pressure * 100.  # in Pa
    T = np.array(temperature, dtype=float)
    Tc = temperature - 273.15  # in C

    Xw = molar_fraction

    Z = 1 - (p / T) * (a0 + a1 * Tc + a2 * Tc ** 2 + (b0 + b1 * Tc) * Xw + \
                       (c0 + c1 * Tc) * Xw ** 2) + (p / T) ** 2 * (d0 + d1 * Xw ** 2)
    return(Z)

def molar_fraction_water_vapour(pressure, temperature, relative_humidity):
    """ Molar fraction of water vapor. 
    
    Parameters
    ----------
    pressure: float
       Total pressure [hPa]
    temperature: float
       Atmospehric temperature [K] 
    relative_humidity:
       Relative humidity from 0 to 100 [%]
   
    """
    # Convert units
    p = pressure  # In hPa
    h = relative_humidity / 100.  # From 0 to 1

    # Calculate water vapor pressure
    f = enhancement_factor_f(pressure, temperature)
    svp = saturation_vapor_pressure(temperature)

    p_wv = h * f * svp  # Water vapor pressure

    Xw = p_wv / p

    return(Xw)


def enhancement_factor_f(pressure, temperature):
    """ Enhancement factor.
    
    Parameters
    ----------
    pressure: float
       Atmospheric pressure [hPa]
    temperature: float
       Atmospehric temperature [K]

    """
    T = temperature
    p = pressure * 100.  # In Pa

    f = 1.00062 + 3.14e-8 * p + 5.6e-7 * (T - 273.15) ** 2

    return(f)


def saturation_vapor_pressure(temperature):
    """ Saturation vapor pressure of water of moist air.
    
    Note: In original documentation, this was specified as the saturation pressure of 
    pure water vapour. This seems wrong. 
    
    
    Parameters
    ----------
    temperature: float
       Atmospheric temperature [K] 
    
    Returns
    -------
    E: float
       Saturation vapor pressure [hPa]
             
    References
    ----------
    Ciddor, P. E.: Refractive index of air: new equations for the visible and near 
    infrared, Appl. Opt., 35(9), 1566-1573, doi:10.1364/AO.35.001566, 1996.
    
    Davis, R. S.: Equation for the Determination of the Density of 
    Moist Air (1981/91), Metrologia, 29(1), 67, doi:10.1088/0026-1394/29/1/008, 1992.
    
    """
    T = temperature
    E = np.exp(1.2378847e-5 * T ** 2 - 1.9121316e-2 * T +
               33.93711047 - 6343.1645 / T)
    return(E / 100.)  # In hPa

