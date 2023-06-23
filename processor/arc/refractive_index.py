""" This file includes functions to calculate the refractive index of air
according to Ciddor (1996, 2002), summarized by Tomasi et al. (2005).
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from .constants import R, eps_o, k_b, R, NA

def n_air(wavelength, pressure, temperature, C, relative_humidity):
    """ Calculate the refractive index of air. 
    
    Parameters
    ----------
    wavelength : float
       Light wavelegnth [nm]
    pressure : float
       Atmospheric pressure [hPa]
    temperature : float
       Atmospehric temperature [K]
    C : float
       Concentration of CO2 [ppmv]
    relative_humidity : float
       Relative humidity from 0 to 100 [%]
    
    Returns
    -------
    n_air : float
       Refractive index of air.
    """

    Xw = molar_fraction_water_vapour(pressure, temperature, relative_humidity)

    rho_axs, _, _ = moist_air_density(1013.25, 288.15, C, 0)
    rho_ws, _, _ = moist_air_density(13.33, 293.15, 0, 1)  # C not relevant

    _, rho_a, rho_w = moist_air_density(pressure, temperature, C, Xw)

    n_axs = n_standard_air_with_CO2(wavelength, C)
    n_ws = n_water_vapor(wavelength)

    n = 1 + (rho_a / rho_axs) * (n_axs - 1) + (rho_w / rho_ws) * (n_ws - 1)

    return n


def moist_air_density(pressure, temperature, C, Xw):
    """ Calculate the moist air density using the BIPM (Bureau International des
    Poids et Mesures) 1981/91 equation. See Tomasi et al. (2005), eq. 12.
    
    Parameters
    ----------
    pressure: float
       Total pressure [hPa]
    temperature: float
       Atmospehric temperature [K]
    C: float
       CO2 concentration [ppmv]
    Xw: float
       Molar fraction of water vapor
    """
    Ma = molar_mass_dry_air(C)  # in kg/mol--  Molar mass  dry air.
    Mw = 0.018015  # in kg/mol -- Molar mass of water vapour. 

    Z = compressibility_of_moist_air(pressure, temperature, Xw)

    P = pressure * 100.  # In Pa
    T = temperature

    rho = P * Ma / (Z * R * T) * (1 - Xw * (1 - Mw / Ma))

    rho_air = (1 - Xw) * P * Ma / (Z * R * T)
    rho_wv = Xw * P * Mw / (Z * R * T)

    return rho, rho_air, rho_wv


def molar_mass_dry_air(C):
    """ Molar mass of dry air, as a function of CO2 concentration.
    
    Parameters
    ----------
    C: float
       CO2 concentration [ppmv]
    
    Returns
    -------
    Ma: float
       Molar mass of dry air [km/mol]
    """
    C1 = 400.

    Ma = 10 ** -3 * (28.9635 + 12.011e-6 * (C - C1))

    return Ma


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

    return Xw


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
    return E / 100.  # In hPa


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
    return Z


def n_standard_air(wavelength):
    """The refractive index of air at a specific wavelength with CO2 concentration 450 ppmv. 
    
    Calculated for standard air at T = 15C, P=1013.25hPa, e = 0, C=450ppmv. 
    (see Tomasi, 2005, eg. 17).
      
    Parameters
    ----------
    wavelength : float
       Wavelength [nm]
    
    Returns
    -------
    ns : float
       Refractivity of standard air with C = 450ppmv
    """

    wl_micrometers = wavelength / 1000.0  # Convert nm to um

    s = 1 / wl_micrometers  # the reciprocal of wavelength
    c1 = 5792105.
    c2 = 238.0185
    c3 = 167917.
    c4 = 57.362
    ns = 1 + (c1 / (c2 - s ** 2) + c3 / (c4 - s ** 2)) * 1e-8
    return ns


def n_standard_air_with_CO2(wavelength, C):
    """ The refractive index of air at a specific wavelength including random CO2. 
    
    Calculated for standard air at T = 15C, P=1013.25hPa, e = 0. 
    (see Tomasi, 2005, eq. 18)
      
    Parameters
    ----------
    wavelength : float
       Wavelength [nm]
    C : float
       CO2 concentration [ppmv]
       
    Returns
    -------
    n_axs : float
       Refractive index of air for the specified CO2 concentration.
    """
    C2 = 450.  # ppmv

    n_as = n_standard_air(wavelength)

    n_axs = 1 + (n_as - 1) * (1 + 0.534e-6 * (C - C2))

    return n_axs


def n_water_vapor(wavelength):
    """ Refractive index of water vapour. 

    Calculated for T = 20C, e=1333Pa  (see Tomasi, 2005, eq. 19)
    
    Parameters
    ----------
    wavelength: float
       Wavelength [nm]
    
    Returns
    -------
    n_wv : float
       Refractive index of water vapour.
    """
    wl_micrometers = wavelength / 1000.0  # Convert nm to um

    s = 1 / wl_micrometers  # the reciprocal of wavelength

    c1 = 1.022
    c2 = 295.235
    c3 = 2.6422
    c4 = 0.032380  # Defined positive
    c5 = 0.004028

    n_ws = 1 + c1 * (c2 + c3 * s ** 2 - c4 * s ** 4 + c5 * s ** 6) * 1e-8

    return n_ws

def n_N2(wavelength, temperature=293.15, pressure=1013.25, method='combined', extrapolate=False):
    """ Refractive index of nitrogen (N2). 

    Calculated for the given number density (temperature  & pressure) and wavelength
    
    Parameters
    ----------
    wavelength : float
       Wavelength in vacuum [nm]
      
    temperature :  float
       Temperature of the gas in Kelvin
       
       Defaults to 273.15K (0C) (STP conditions)
      
    pressure : float
       Pressure of the gas in hecto Pascal (hPa)
       
       Defaults to 1000 hPa (STP conditions)
       
       It is recommended to use the partial pressure of the gas in a mixture
       
    method : string
       One of the following:  griesmann_burnett, boerzsoenyi,
       peck_and_khanna, combined
       
       In all cases the calculations are performed for the given pressure
       and temperature conditions and not the ones where the formulas have
       been defined!
       
       griesmann_burnett:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 145-270nm
           
           U. Griesmann and J. H. Burnett. Refractivity of nitrogen gas in the 
           vacuum ultraviolet, Opt. lett. 24, 1699-1701 (1999)

       boerzsoenyi:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 400-1000nm
                     
           A. Börzsönyi, Z. Heiner, M. P. Kalashnikov, A. P. Kovács, 
           and K. Osvay, Dispersion measurement of inert gases and gas 
           mixtures at 800 nm, Appl. Opt. 47, 4856-4863 (2008)

       peck_khanna:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 468-2058nm
           
           E. R. Peck and B. N. Khanna. Dispersion of nitrogen, 
           J. Opt. Soc. Am. 56, 1059-1063 (1966)
    
       combined:
           
           From 145 to 270nm the formula of griesmann_and_burnett is applied.
           
           From 270 to 400nm the formula of boerzsoenyi is extrapollated
           
           From 400 to 1000nm the formula of Boerzsoenyi is applied. 

           From 1000 to 2058nm the formula of Peck_Khanna is applied.

    extrapolate : bool
       If False nan values are returned for wavelengths out of the valid limits
       of the methods

    Returns
    -------
    n : float
       Refractive index of water vapour.
    
    a : float
       Isotropic polarizability of N2 in farad * m3 (SI)
       
    N : float
       Gas number density of N2 in m-3

    """

    gas = "N2"
    if method == 'griesmann_burnett':
        wv_llim = 145.
        wv_ulim = 270.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = griesmann_and_burnett_N2_isotropic_polarizability(wavelength)

    elif method == 'boerzsoenyi':
        wv_llim = 400.
        wv_ulim = 1000.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = boerzsoenyi_N2_isotropic_polarizability(wavelength)

    elif method == 'peck_khanna':
        wv_llim = 468.
        wv_ulim = 2058.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = peck_and_khanna_N2_isotropic_polarizability(wavelength)

    elif method == 'combined':
        wv_llim = 145.
        wv_change_1 = 288.
        wv_change_2 = 1000.
        wv_ulim = 2058.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a_griesmann_and_burnett = griesmann_and_burnett_N2_isotropic_polarizability(wavelength)
        a_boerzsoenyi = boerzsoenyi_N2_isotropic_polarizability(wavelength)
        a_peck_and_khanna = peck_and_khanna_N2_isotropic_polarizability(wavelength)

        a = np.hstack((a_griesmann_and_burnett[wavelength < wv_change_1],
                       a_boerzsoenyi[(wavelength >= wv_change_1) & (wavelength <= wv_change_2)],
                       a_peck_and_khanna[wavelength > wv_change_2]))
                
        if len(a) == 1:
            a = a[0]
            
    else:
        raise ValueError('Method {} not implemented for n_N2. Should be one of "griesmann_and_burnett", "boerzsoenyi", "peck_and_khanna" or "combined"'.format(method))

    N = 100. * pressure / (k_b * temperature)

    n_f = np.sqrt((3. * eps_o + 2. * a * N) / (3. * eps_o - a * N))

    if not extrapolate:
        mask = (wavelength < wv_llim) | (wavelength > wv_ulim)
        n_f[mask] = np.nan  # nan values in wavelengths out of the valid ranges

    return n_f, a, N


def n_O2(wavelength, temperature=293.15, pressure=1013.25, method='combined', extrapolate=False):
    """ Refractive index of oxygen (O2). 

    Calculated for the given number density (temperature  & pressure) and wavelength
    
    Parameters
    ----------
    wavelength : float  or np.ndarray
       Wavelength in vacuum [nm]
      
    temperature : float
       Temperature of the gas in Kelvin
       
       Defaults to 293.15K (20C) (STP conditions according to NIST)
      
    pressure : float
       Pressure of the gas in hecto Pascal (hPa)
       
       Defaults to 1013.25 hPa  (STP conditions according to NIST)
       
       It is recommended to use the partial pressure of the gas in a mixture
       
    method : string
       One of the following: zhang, smith, combined
       
       In all cases the calculations are performed for the given pressure
       and temperature conditions and not the ones where the formulas have
       been defined!
       
       smith:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 185-288nm
           
           P. L. Smith, M. C. E. Huber, W. H. Parkinson. Refractivities of 
           H2, He, O2, CO, and Kr for 168≤λ≤288 nm,
           Phys Rev. A 13, 199-203 (1976)
           
       zhang:
           
           Original formula corresponds to T = 293.15K (20C), e = 1013.25 hPa
           
           Valid in the range 400-1800nm
           
           The correction of Kren is applied to the original formula of Zhang
           
           J. Zhang, Z. H. Lu, and L. J. Wang. Precision refractive index 
           measurements of air, N2, O2, Ar, and CO2 with a frequency comb, 
           Appl. Opt. 47, 3143-3151 (2008)
           
           P. Křen. Comment on "Precision refractive index measurements of air,
           N2, O2, Ar, and CO2 with a frequency comb", 
           Appl. Opt. 50, 6484-6485 (2011)
           
       combined:
           
           From 185 to  288nm the formula of Smith is applied. 
           
           From 288 to  400nm the formula of Smith is extrapollated
           
           From 400 to 1800nm the formula of Smith is applied.

    extrapolate : bool
       If False nan values are returned for wavelengths out of the valid limits
       of the methods
       
       
    Returns
    -------
    n : float
       Refractive index of water vapour.
    
    a : float
       Isotropic polarizability of O2  in farad * m3 (SI)
       
    N : float
       Gas number density of O2 in m-3
    """

    gas = 'O2'
    
    if method == 'smith':
        wv_llim = 185.
        wv_ulim = 288.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = smith_O2_isotropic_polarizability(wavelength)

    elif method == 'zhang':
        wv_llim = 400.
        wv_ulim = 1800.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = zhang_O2_isotropic_polarizability(wavelength)

    elif method == 'combined':
        wv_llim = 185.
        wv_change = 288.
        wv_ulim = 1800.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a_smith = smith_O2_isotropic_polarizability(wavelength)
        a_zhang = zhang_O2_isotropic_polarizability(wavelength)

        a = np.hstack((a_smith[wavelength < wv_change], a_zhang[wavelength >= wv_change]))
        
        if len(a) == 1:
            a = a[0]
            
    else:
        raise ValueError('Method {} not implemented for n_O2. Should be one of "zhang", "smith" or "combined"'.format(method))

    N = 100. * pressure / (k_b * temperature)

    n_f = np.sqrt((3. * eps_o + 2. * a * N) / (3. * eps_o - a * N))

    if not extrapolate:
        mask = (wavelength < wv_llim) | (wavelength > wv_ulim)
        n_f[mask] = np.nan  # nan values in wavelenghts out of the valid ranges

    return n_f, a, N


def n_Ar(wavelength, temperature=293.15, pressure=1013.25, method='combined', extrapolate=False):
    """ Refractive index of Argon (Ar). 

    Calculated for the given number density (temperature  & pressure) and wavelength
    
    Parameters
    ----------
    wavelength : float
       Wavelength in vacuum [nm]
      
    temperature : float
       Temperature of the gas in Kelvin
       
       Defaults to 293.15K (20C) (STP conditions accodring to NIST)
      
    pressure : float
       Pressure of the gas in hecto Pascal (hPa)
       
       Defaults to 1013.25 hPa (STP conditions according to NIST)
       
       It is recommended to use the partial pressure of the gas in a mixture
       
    method: string
       One of the following: bideau_mehu_larsen, boerzsoenyi, peck_fisher, combined
       
       In all cases the calculations are performed for the given pressure
       and temperature conditions and not the ones where the formulas have
       been defined!
       
       bideau_mehu_larsen:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 141-567nm
           
           A. Bideau-Mehu, Y. Guern, R. Abjean, A. Johannin-Gilles. 
           Measurement of refractive indices of neon, argon, krypton and xenon 
           in the 253.7-140.4 nm wavelength range. Dispersion relations and 
           estimated oscillator strengths of the resonance lines. J. Quant. 
           Spectrosc. Rad. Transfer 25, 395-402 (1981)
           
           T. Larsén. Beitrag zur Dispersion der Edelgase. Z. 
           Physik 88, 389-394 (1934)
           
           *Sellmeier formula is derived by the authors of ref. 1 
           using their own data in the 0.1404-0.2537 μm range combined with 
           data from ref. 2 at longer wavelengths.

       boerzsoenyi:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 400-1000nm
                     
           A. Börzsönyi, Z. Heiner, M. P. Kalashnikov, A. P. Kovács, 
           and K. Osvay, Dispersion measurement of inert gases and gas 
           mixtures at 800 nm, Appl. Opt. 47, 4856-4863 (2008)
           
       peck_fisher:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 468-2058nm
                     
           E. R. Peck and D. J. Fisher. Dispersion of argon, J. 
           Opt. Soc. Am. 54, 1362-1364 (1964)

       combined:
           
           From 181 to 1694nm the formula of Bideau-Mehu_Larsen is applied.

           From 400 to 1000nm the formula of Boerzsoenyi is applied.           

           From 468 to 2058nm the formula of Peck_and_Fischer is applied.           

    extrapolate: bool
       If False nan values are returned for wavelengths out of the valid limits
       of the methods
       
       
    Returns
    -------
    n : float
       Refractive index of water vapour.
    
    a : float
       Isotropic polarizability of Ar in farad * m3 (SI)
       
    N : float
       Gas number density of Ar in m-3

    """

    gas = 'Ar'
    
    if method == 'bideau_mehu_larsen':
        wv_llim = 141.
        wv_ulim = 567.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = bideau_mehu_larsen_Ar_isotropic_polarizability(wavelength)

    elif method == 'boerzsoenyi':
        wv_llim = 400.
        wv_ulim = 1000.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = boerzsoenyi_Ar_isotropic_polarizability(wavelength)

    elif method == 'peck_fisher':
        wv_llim = 468.
        wv_ulim = 2058.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = peck_and_fisher_Ar_isotropic_polarizability(wavelength)

    elif method == 'combined':
        wv_llim = 141.
        wv_change_1 = 400.
        wv_change_2 = 1000.
        wv_ulim = 2058.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a_bideau_mehu_larsen = bideau_mehu_larsen_Ar_isotropic_polarizability(wavelength)
        a_boerzsoenyi = boerzsoenyi_Ar_isotropic_polarizability(wavelength)
        a_peck_and_fisher = peck_and_fisher_Ar_isotropic_polarizability(wavelength)

        a = np.hstack((a_bideau_mehu_larsen[wavelength < wv_change_1],
                       a_boerzsoenyi[(wavelength >= wv_change_1) & (wavelength <= wv_change_2)],
                       a_peck_and_fisher[wavelength > wv_change_2]))
        
        if len(a) == 1:
            a = a[0]

    else:
        raise ValueError(
            'Method {} not implemented for n_Ar. Should be one of "bideau_mehu_larsen", "boerzsoenyi", "peck_fisher" or "combined"'.format(
                method))

    N = 100. * pressure / (k_b * temperature)

    n_f = np.sqrt((3. * eps_o + 2. * a * N) / (3. * eps_o - a * N))

    if not extrapolate:
        mask = (wavelength < wv_llim) | (wavelength > wv_ulim)
        n_f[mask] = np.nan  # nan values in wavelenghts out of the valid ranges

    return n_f, a, N

def n_CO2(wavelength, temperature=293.15, pressure=1013.25, method='combined', extrapolate=False):
    """ Refractive index of carbon nitrogen (CO2). 

    Calculated for the given number density (temperature  & pressure) and wavelength
    
    Parameters
    ----------
    wavelength : float
       Wavelength in vacuum [nm]
      
    temperature : float
       Temperature of the gas in Kelvin
       
       Defaults to 293.15K (20C) (STP conditions according to NIST)
      
    pressure : float
       Pressure of the gas in hecto Pascal (hPa)
       
       Defaults to 1013.25 hPa  (STP conditions according to NIST)
       
       It is recommended to use the partial pressure of the gas in a mixture
       
    method : string
       One of the following: bideau_mehu_larsen, old, combined
       
       In all cases the calculations are performed for the given pressure
       and temperature conditions and not the ones where the formulas have
       been defined!
       
       bideau_mehu:
           
           Original formula corresponds to T = 273.15K (0C), e = 1013.25 hPa
           
           Valid in the range 141-566nm
           
           A. Bideau-Mehu, Y. Guern, R. Abjean and A. Johannin-Gilles. 
           Interferometric determination of the refractive index of carbon 
           dioxide in the ultraviolet region, Opt. Commun. 9, 432-434 (1973)
           

       old:
           
           Original formula corresponds to T = 273.15K (0C), p = 1013.25 hPa
           
           Valid in the range 481-1817nm
           
           J. G. Old, K. L. Gentili, and E. R. Peck. Dispersion of Carbon 
           Dioxide, J. Opt. Soc. Am. 61, 89-90 (1971)

       combined:
           
           From 181 to 1694nm the formula of Bideau-Mehu is applied.

           From 1694 to 1817nm the formula of Old is applied.           
 
    extrapolate: bool
       If False nan values are returned for wavelengths out of the valid limits
       of the methods
       
       
    Returns
    -------
        
    n : float
       Refractive index of water vapour.
    
    a : float
       Isotropic polarizability of CO2 in farad * m3 (SI)
       
    N : float
       Gas number density of CO2 in m-3

    """

    gas = 'CO2'
    
    if method == 'bideau_mehu':
        wv_llim = 181.
        wv_ulim = 1694.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = bideau_mehu_CO2_isotropic_polarizability(wavelength)

    elif method == 'old':
        wv_llim = 481.
        wv_ulim = 1817.
       
        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = old_CO2_isotropic_polarizability(wavelength)

    elif method == 'combined':
        wv_llim = 181.
        wv_change_1 = 1694.
        wv_ulim = 1817.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a_bideau_mehu = bideau_mehu_CO2_isotropic_polarizability(wavelength)
        a_old = old_CO2_isotropic_polarizability(wavelength)

        a = np.hstack((a_bideau_mehu[wavelength <= wv_change_1],
                       a_old[wavelength > wv_change_1]))
                
        if len(a) == 1:
            a = a[0]

    else:
        raise ValueError(
            'Method {} not implemented for n_CO2. Should be one of "bideau_mehu_larsen", "old", "peck_and_khanna" or "combined"'.format(
                method))

    N = 100. * pressure / (k_b * temperature)

    n_f = np.sqrt((3. * eps_o + 2. * a * N) / (3. * eps_o - a * N))

    if not extrapolate:
        mask = (wavelength < wv_llim) | (wavelength > wv_ulim)
        n_f[mask] = np.nan  # nan values in wavelengths out of the valid ranges

    return n_f, a, N


def n_H2O(wavelength, temperature=293.15, pressure=1013.25, method='combined', extrapolate=False):
    """ Refractive index of water vapor (H2O). 

    Calculated for the given number density (temperature  & pressure) and wavelength
    
    Parameters
    ----------
    wavelength : float
       Wavelength in vacuum [nm]
      
    temperature : float
       Temperature of the gas in Kelvin
       
       Defaults to 293.15K (20C) (STP conditions according to NIST)
      
    pressure : float
       Pressure of the gas in hecto Pascal (hPa)
       
       Defaults to 1013.25 hPa  (STP conditions according to NIST)
       
       It is recommended to use the partial pressure of the gas in a mixture
       
    method: string
       One of the following: cidor, combined
       
       In all cases the calculations are performed for the given pressure
       and temperature conditions and not the ones where the formulas have
       been defined!
       
       cidor:
           
           Original formula corresponds to T = 293.15K (20C), e = 13.33 hPa
           
           Valid in the range 350-1200nm (claimed by Cidor to be able to work 
                                          out of the original measurement range)
           
           P. E. Ciddor, “Refractive index of air: new equations for the
           visible and near infrared,” Appl. Opt. 35, 1566 –1573 (1996).
    
       combined:
           
           From 350 to 1200nm the formula of Cidor is applied.
           
    extrapolate: bool
       If False nan values are returned for wavelengths out of the valid limits
       of the methods
       
       
    Returns
    -------
    n : float
       Refractive index of water vapour.
    
    a : float
       Isotropic polarizability of H2O in farad * m3 (SI) 
       
    N : float
       Gas number density of H2O in m-3

    """

    gas = 'H2O'
    
    if method == 'cidor':
        wv_llim = 350.
        wv_ulim = 1200.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = cidor_wv_isotropic_polarizability(wavelength)

    elif method == 'combined':
        wv_llim = 350.
        # wv_change_1 = 1200.
        wv_ulim = 1200.

        if wavelength < wv_llim or wavelength > wv_ulim:
            raise Exception(f'The selected wavelength ({wavelength}nm) is out of the {method} method limits for {gas}: {wv_llim} to {wv_ulim}nm')
       
        a = cidor_wv_isotropic_polarizability(wavelength)

        # a = np.hstack((a_bideau_mehu[wavelength <= wv_change_1],
        # a_old[wavelength > wv_change_1]))
    else:
        raise ValueError(
            'Method {} not implemented for n_H2O. Should be one of "cidor" or "combined"'.format(
                method))

    N = 100. * pressure / (k_b * temperature)

    n_f = np.sqrt((3. * eps_o + 2. * a * N) / (3. * eps_o - a * N))

    if not extrapolate:
        mask = (wavelength < wv_llim) | (wavelength > wv_ulim)
        n_f[mask] = np.nan  # nan values in wavelenghts out of the valid ranges

    return n_f, a, N



# Functions for the calculations of the isotropic polarizabilities for each 
# atmospheric gas species

def zhang_O2_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of O2 

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 101325
    T_o = 293.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 1.181494E-4
    c2 = 9.708931E-3
    c3 = 75.4
    n_o = 1. + c1 + c2 / (c3 - np.power(x, -2))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def smith_O2_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of O2  

    """

    e_o = 101325
    T_o = 273.15

    N_o = e_o / (k_b * T_o)  # Number density in m^-3 in the measurement conditions of the RI

    wv_o = np.array([288.24, 263.21, 252.93, 252.49, 252.00, 251.69, 251.51,
                     250.77, 243.59, 221.74, 221.16, 220.87, 212.48, 205.88,
                     198.90, 198.64, 198.32, 198.06, 197.92, 197.76, 191.85,
                     184.30])

    n_o = 1. + np.array([291.6, 298.6, 302.4, 302.7, 302.7, 302.9, 303.0, 303.3,
                         306.5, 320.1, 320.6, 320.8, 328.5, 336.1, 346.2, 346.7,
                         347.1, 347.6, 347.9, 348.2, 361.2, 379.2]) * 1E-6

    f_n = interp1d(wv_o, n_o, bounds_error=False, fill_value='extrapolate')

    a = (3. * eps_o) * (1. / N_o) * (np.power(f_n(x), 2) - 1.) / (np.power(f_n(x), 2) + 2.)

    return a


def bideau_mehu_CO2_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of CO2 

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 101325
    T_o = 273.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 6.99100E-2
    c2 = 166.175
    c3 = 1.44720E-3
    c4 = 79.609
    c5 = 6.42941E-5
    c6 = 56.3064
    c7 = 5.21306E-5
    c8 = 46.0196
    c9 = 1.46847E-6
    c10 = 0.0584738

    n_o = 1. + c1 / (c2 - np.power(x, -2)) + c3 / (c4 - np.power(x, -2)) + \
          c5 / (c6 - np.power(x, -2)) + c7 / (c8 - np.power(x, -2)) + \
          c9 / (c10 - np.power(x, -2))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def old_CO2_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of CO2 

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 101325
    T_o = 273.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 1.54489E-6
    c2 = 5.84738E-2
    c3 = 8.3091927E-2
    c4 = 210.9241
    c5 = 2.8764190E-3
    c6 = 60.122959

    n_o = 1. + c1 / (c2 - np.power(x, -2)) + c3 / (c4 - np.power(x, -2)) + \
          c5 / (c6 - np.power(x, -2))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def griesmann_and_burnett_N2_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of N2 

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 101325
    T_o = 273.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 1.9662731
    c2 = 22086.66
    c3 = 2.7450825E-2
    c4 = 133.85688

    n_o = 1. + c1 / (c2 - np.power(x, -2)) + c3 / (c4 - np.power(x, -2))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def boerzsoenyi_N2_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of N2 

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 100000
    T_o = 273.00

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 39209.95E-8
    c2 = 1146.24E-6
    c3 = 18806.48E-8
    c4 = 13.476E-3

    n_o = np.sqrt(1. + c1 * np.power(x, 2) / (np.power(x, 2) - c2) + c3 * np.power(x, 2) / (np.power(x, 2) - c4))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def peck_and_khanna_N2_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of N2 

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 101325
    T_o = 273.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 6.8552E-5
    c2 = 3.243157E-2
    c3 = 144.
    n_o = 1. + c1 + c2 / (c3 - np.power(x, -2))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def cidor_wv_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of water vapor

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 1333
    T_o = 293.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 1.022
    c2 = 295.235
    c3 = 2.6422
    c4 = 0.032380
    c5 = 0.004028

    n_o = 1. + c1 * (c2 + c3 * np.power(x, -2) - c4 * np.power(x, -4) + \
                     c5 * np.power(x, -6)) * 1e-8

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def bideau_mehu_larsen_Ar_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of Argon

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 101325
    T_o = 273.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 2.50141E-3
    c2 = 91.012
    c3 = 5.00283E-4
    c4 = 87.892
    c5 = 5.22343E-2
    c6 = 214.02

    n_o = 1. + c1 / (c2 - np.power(x, -2)) + c3 / (c4 - np.power(x, -2)) + \
          c5 / (c6 - np.power(x, -2))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def boerzsoenyi_Ar_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of Argon

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 100000
    T_o = 273.00

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 20332.29E-8
    c2 = 206.12E-6
    c3 = 34458.31E-8
    c4 = 8.066E-3

    n_o = np.sqrt(1. + c1 * np.power(x, 2) / (np.power(x, 2) - c2) + \
                  c3 * np.power(x, 2) / (np.power(x, 2) - c4))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a


def peck_and_fisher_Ar_isotropic_polarizability(x):
    """
    
    Parameters
    ----------
    x: float
       Wavelength in vacuum [nm]
       
       
    Returns
    -------
    
    a : float
       Isotropic polarizability of Argon

    """

    x = x * 1E-3  # the original formula is expressed in μm

    e_o = 101325
    T_o = 273.15

    N_o = e_o / (k_b * T_o)  # Number density in the measurement conditions of the RI

    c1 = 6.7867E-5
    c2 = 3.0182943E-2
    c3 = 144

    n_o = 1. + c1 + c2 / (c3 - np.power(x, -2))

    a = (3. * eps_o) * (1. / N_o) * (np.power(n_o, 2) - 1.) / (np.power(n_o, 2) + 2.)

    return a

def nonideal_gas(gas, pressure = 1013.25, temperature = 293.15): # TO DO not ready yet
    
    pressure = pressure * 1E2
    # Van der Waals constant from https://en.wikipedia.org/wiki/Van_der_Waals_constants_(data_page)
    if gas == 'N2':
        a = 1.370
        b = 0.03870
    elif gas == 'O2':
        a = 1.382
        b = 0.03186
    elif gas == 'Ar':
        a = 1.355
        b = 0.03201
    elif gas == 'CO2':
        a = 3.640
        b = 0.04267
    elif gas == 'H2O':
        a = 5.536
        b = 0.03049
        
    a = a * 1E-1 # Convert to SI
    b = b * 1E-3 # Convert to SI
    
    number_density_ideal = ideal_gas(pressure=pressure, temperature=temperature)
    
    number_density = \
        fsolve(van_der_waals, x0 = number_density_ideal, args = (pressure,temperature,a,b),xtol = 1E-8)[0] #in m-3
            
    return number_density

def ideal_gas(pressure = 1013.25, temperature = 293.15):
    
    pressure = pressure * 1E2
    
    number_density = pressure / (k_b * temperature) #in m-3
    
    return number_density
    

def van_der_waals(N, pressure, temperature, a, b):
    
    f = (pressure + a*(N/NA)**2)*(1. - b * N/NA) - (N/NA) * R * temperature 
    
    return f