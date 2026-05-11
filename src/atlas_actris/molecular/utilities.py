""" Utility functions used in all modules. """

from molecular.constants import k_b
import numpy as np

def number_density_at_pt(pressure, temperature, relative_humidity, ideal=False):
    """ Calculate the number density for a given temperature and pressure,
    taking into account the compressibility of air.
    
    Parameters
    ----------
    pressure: float or array
       Pressure in Pa
    temperature: float or array
       Temperature in K
    relative_humidity: float or array 
       The relative humidity of air (between 0 and 1)
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

    p_pa = pressure  # Pressure in pascal

    n = p_pa / (Z * temperature * k_b)
    
    return n
    
def saturation_vapour_pressure(pressure, temperature):
    
    """
Saturation Vapour Pressure (Buck Equation, 1981)

This formula calculates the saturation vapour pressure over liquid water.
It is a refined Magnus-type equation with an enhancement factor that corrects
for the effect of ambient air pressure on water vapour (non-ideal behaviour).

Formula:
    e_s = (1.0007 + 3.46e-6 * p) * 6.1121 * exp(17.502 * t / (240.97 + t))

Where:
    t = air temperature in °C  (T[K] - 273.15)
    p = ambient pressure in hPa (or mbar)
    e_s = saturation vapour pressure in hPa

Source:
    Buck, A. L. (1981). New equations for computing vapor pressure and enhancement factor.
    Journal of Applied Meteorology and Climatology, 20(12), 1527–1532.

Wikipedia reference (general Magnus formula background):
    https://en.wikipedia.org/wiki/Magnus_formula
    
    Parameters: 
    -------
    pressure: scalar float or 1D array
        Pressure in Pa
    temperature: scalar float or 1D array
        Temperature in K

    Returns
    -------
    e_s: Saturation Vapour pressure in Pa
    
"""
    
    P = pressure / 100. #in hPa
    
    T = temperature - 273.15 #in K
    
    e_s = (1.0007 + 3.46 * 1E-6 * P) * 6.1121 * \
        np.exp(17.502 * T /(240.97 + T)) #in hPa
       
    e_s = e_s / 100.
    
    return(e_s)

def rh_to_pressure(rh, temperature):
    """ Convert relative humidity to water vapour partial pressure.
    
    Parameters
    ----------
    rh: float
       Relative humidity from 0 to 1
    temperature: float
       Temperature [K]
       
    Returns
    -------
    p_wv: float
       Water vapour pressure [hPa].
    """
    svp = saturation_vapor_pressure(temperature)
    h = rh
    
    p_wv = h * svp
    return p_wv


def pressure_to_rh(partial_pressure, temperature):
    """ Convert water vapour partial pressure to relative humidity.

    Parameters
    ----------
    partial_pressure: float
       Water vapour partial pressure [hPa] 
    temperature: float
       Temperature [K]

    Returns
    -------
    rh: float
       Relative humidity from 0 to 1 [%].
    """
    svp = saturation_vapor_pressure(temperature)

    rh = partial_pressure / svp

    return rh

def molar_fraction_water_vapour(pressure, temperature, relative_humidity):
    """ Molar fraction of water vapor. 
    
    Parameters
    ----------
    pressure: float
       Total pressure [hPa]
    temperature: float
       Atmospehric temperature [K] 
    relative_humidity:
       Relative humidity from 0 to 1 [%]
    """
    # Convert units
    p = pressure  # In hPa
    h = relative_humidity  # From 0 to 1

    # Calculate water vapor pressure
    f = enhancement_factor_f(pressure, temperature)
    svp = saturation_vapor_pressure(temperature)

    p_wv = h * f * svp  # Water vapor pressure

    Xw = p_wv / p

    return Xw


def enhancement_factor_f(pressure, temperature):
    """ Enhancement factor.
    
    Parameters
    ----------
    pressure: float
       Atmospheric pressure [Pa]
    temperature: float
       Atmospehric temperature [K]    
    """
    T = temperature
    p = pressure  # In Pa

    f = 1.00062 + 3.14e-8 * p + 5.6e-7 * (T - 273.15) ** 2

    return f


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
       Saturation vapor pressure [Pa]
             
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
    return E  # In Pa


def compressibility_of_moist_air(pressure, temperature, molar_fraction):
    """ Compressibility of moist air.
    
    Parameters
    ----------
    pressure: float
       Atmospheric pressure [Pa]
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

    p = pressure  # in Pa
    T = np.array(temperature, dtype=float)
    Tc = temperature - 273.15  # in C

    Xw = molar_fraction

    Z = 1 - (p / T) * (a0 + a1 * Tc + a2 * Tc ** 2 + (b0 + b1 * Tc) * Xw + \
                       (c0 + c1 * Tc) * Xw ** 2) + (p / T) ** 2 * (d0 + d1 * Xw ** 2)
    return Z
