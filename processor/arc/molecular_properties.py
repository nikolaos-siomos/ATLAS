""" This module calculates the atmospheric gas properties, including
King's factor, epsilon, gamma squared and depolarization. The calculations
are based on the following publications:

Tomasi, C., Vitale, V., Petkov, B., Lupi, A. and Cacciari, A.: Improved 
algorithm for calculations of Rayleigh-scattering optical depth in standard 
atmospheres, Applied Optics, 44(16), 3320, doi:10.1364/AO.44.003320, 2005.

She, C.-Y. Spectral structure of laser light scattering revisited: 
bandwidths of nonresonant scattering lidars. Appl. Opt. 40, 4875-4884 (2001)

Chance, K. V. & Spurr, R. J. D. Ring effect studies: Rayleigh scattering, 
including molecular parameters for rotational Raman scattering, and the 
Fraunhofer spectrum. Applied Optics 36, 5224 (1997).
"""

import numpy as np


def kings_factor_N2(wavenumber):
    """ Approximates the King's correction factor for a specific wavenumber.

    According to Bates, the agreement with experimental values is
    "rather better than 1 per cent."

    Parameters
    ----------
    wavenumber : float
       Wavenumber (defined as 1/lamda) in m-1
    
    Returns
    -------
    Fk : float
       Kings factor for N2

    Notes
    -----
    The King's factor is estimated as:

    .. math::

       F_{N_2} = 1.034 + 3.17 \cdot 10^{-4} \cdot \lambda^{-2}

    where :math:`\lambda` is the wavelength in micrometers.

    References
    ----------
    Tomasi, C., Vitale, V., Petkov, B., Lupi, A. & Cacciari, A. Improved
    algorithm for calculations of Rayleigh-scattering optical depth in standard
    atmospheres. Applied Optics 44, 3320 (2005).

    Bates, D. R.: Rayleigh scattering by air, Planetary and Space Science, 32(6),
    785-790, doi:10.1016/0032-0633(84)90102-8, 1984.
    """
    lamda_m = 1 / wavenumber
    lamda_um = lamda_m * 10 ** 6  # Convert to micrometers, as in the paper

    Fk = 1.034 + 3.17e-4 * lamda_um ** -2

    return Fk


def kings_factor_O2(wavenumber):
    """ Approximates the King's correction factor for a specific wavelength.

    Parameters
    ----------
    wavenumber : float
       wavenumber (defined as 1/lamda) in m-1 
    
    Returns
    -------
    Fk : float
       King's factor for O2

    Notes
    -----
    The King's factor is estimated as:

    .. math::

       F_{N_2} = 1.096 + 1.385 \cdot 10^{-3} \cdot \lambda^{-2} + 1.448 \cdot 10^{-4} \cdot \lambda^{-4}

    where :math:`\lambda` is the wavelength in micrometers.

    References
    ----------
    Tomasi, C., Vitale, V., Petkov, B., Lupi, A. & Cacciari, A. Improved
    algorithm for calculations of Rayleigh-scattering optical depth in standard
    atmospheres. Applied Optics 44, 3320 (2005).
    """
    lamda_m = 1 / wavenumber
    lamda_um = lamda_m * 10 ** 6  # Convert to micrometers, as in the paper

    Fk = 1.096 + 1.385e-3 * lamda_um ** -2 + 1.448e-4 * lamda_um ** -4

    return Fk


def kings_factor_Ar():
    """ Returns the King's correction factor of Ar.

    According to Tomasi et al., 2005 it's wavelength independent and equal to 1.
    The value is taken from Alms et al., 1975
    Returns
    -------
    Fk : float
       King's factor for Ar

    References
    ----------
    Tomasi, C., Vitale, V., Petkov, B., Lupi, A. & Cacciari, A. Improved
    algorithm for calculations of Rayleigh-scattering optical depth in standard 
    atmospheres. Applied Optics 44, 3320 (2005).

    """
    Fk = 1.0
    return Fk


def kings_factor_CO2():
    """ Returns the King's correction factor of CO2.

    According to Tomasi et al., 2005 it's wavelength independent and equal to 1.15.
    The value is taken from  Alms et al., 1975.
    
    Returns
    -------
    Fk : float
       King's factor for CO2

    References
    ----------
    Tomasi, C., Vitale, V., Petkov, B., Lupi, A. & Cacciari, A. Improved
    algorithm for calculations of Rayleigh-scattering optical depth in standard
    atmospheres. Applied Optics 44, 3320 (2005).

    Alms et al. Measurement of the discpersion in polarizability anisotropies. Journal of Chemical Physics (1975)
    TODO: Fix citation.
    """
    Fk = 1.15
    return Fk


def kings_factor_H2O():
    """ Returns the King's correction factor of H2O.

    According to Tomasi et al., 2005 it's wavelength independent and equal to 1.001.
    The value is taken from Sioris et al., 2002.
    Returns
    -------
    Fk : float
       King's factor for H2O

    References
    ----------
    Tomasi, C., Vitale, V., Petkov, B., Lupi, A. & Cacciari, A. Improved
    algorithm for calculations of Rayleigh-scattering optical depth in standard
    atmospheres. Applied Optics 44, 3320 (2005).

    C. E. Sioris, W. F. J. Evans, R. L. Gattinger, I. C. McDade, D. A. Degenstein,
    and E. J. Llewellyn, "Ground-based Ring-effect measurements with the OSIRIS
    development model," Can. J. Phys. 80, 483-491 (2002).
    """
    Fk = 1.001
    return Fk


def kings_factor_atmosphere(wavelength, C=300., p_e=0., p_t=1013.25):
    """ Calculates the King's factor for the atmosphere.

    The calculations assume fixed concentrations for :math:`N_2`, :math:`C_{O_2}`, and :math:`Ar`, equal to:

    .. math::

       c_{N_{2}} = 78.084\%

       c_{O_{2}} = 20.946\%

       c_{Ar} = 0.934\%

    The concentration of :math:`CO_2` and :math:`H_2O` can be defined as input.

    Parameters
    ----------
    wavelength: float or array of floats
       Wavelength in nm
    C: float
       CO2 concentration in ppmv
    p_e: float
       water-vapor pressure in hPa
    p_t: float
       total air pressure in hPa

    Returns
    -------
    Fk: float or array of floats
       Total atmospheric King's factor
    """
    wavelength = np.array(wavelength)

    if not np.all((wavelength >= 200) & (wavelength <= 4000)):
        raise ValueError("King's factor formula is only valid from 0.2 to 4um.")

    # Calculate wavenumber
    lamda_m = wavelength * 10 ** -9
    wavenumber = 1 / lamda_m

    # Individual kings factors
    F_N2 = kings_factor_N2(wavenumber)
    F_O2 = kings_factor_O2(wavenumber)
    F_ar = kings_factor_Ar()
    F_CO2 = kings_factor_CO2()
    F_H2O = kings_factor_H2O()

    # Individual concentrations
    c_n2 = 0.78084
    c_o2 = 0.20946
    c_ar = 0.00934
    c_co2 = 1e-6 * C
    c_h2o = p_e / p_t

    # Total concentration
    c_tot = c_n2 + c_o2 + c_ar + c_co2 + c_h2o

    F_k = (c_n2 * F_N2 + c_o2 * F_O2 + c_ar * F_ar + c_co2 * F_CO2 + c_h2o * F_H2O) / c_tot
    return F_k


def epsilon_N2(wavenumber):
    """ Returns the epsilon parameters of N2 for a given wavelength.
    
    epsilon = (gamma / alpha )**2 is the relative anisotropy.
    
    Ra according to  She (2000).
    
    She, C.-Y. Spectral structure of laser light scattering revisited: 
    bandwidths of nonresonant scattering lidars. 
    Appl. Opt. 40, 4875-4884 (2001)
    
    Parameters
    ----------
    wavenumber : float
       wavenumber (defined as 1/lamda) in m-1 

    Returns
    -------
    e : float
       ...
    """

    Fk = kings_factor_N2(wavenumber)

    e = (Fk - 1.) * 9 / 2.

    return e


def epsilon_O2(wavenumber):
    """ Returns the epsilon parameters of O2 for a given wavelength.
    
    epsilon = (gamma / alpha )**2 is the relative ansiotropy.
    
    Ra according to  She (2000).
    
    She, C.-Y. Spectral structure of laser light scattering revisited: 
    bandwidths of nonresonant scattering lidars. 
    Appl. Opt. 40, 4875-4884 (2001).
    
    Parameters
    ----------
    wavenumber : float
       wavenumber (defined as 1/lamda) in m-1 
    """

    Fk = kings_factor_O2(wavenumber)

    e = (Fk - 1) * 9 / 2.

    return e

def epsilon_Ar(wavenumber):
    """ Returns the epsilon parameters of Ar for a given wavelength.
    
    epsilon = (gamma / alpha )**2 is the relative ansiotropy.
    
    Ra according to  She (2000).
    
    She, C.-Y. Spectral structure of laser light scattering revisited: 
    bandwidths of nonresonant scattering lidars. 
    Appl. Opt. 40, 4875-4884 (2001).
    
    Parameters
    ----------
    wavenumber : float
       wavenumber (defined as 1/lamda) in m-1 
    """

    Fk = kings_factor_Ar()

    e = (Fk - 1) * 9 / 2.

    return e

def epsilon_CO2(wavenumber):
    """ Returns the epsilon parameters of CO2 for a given wavelength.
    
    epsilon = (gamma / alpha )**2 is the relative ansiotropy.
    
    Ra according to  She (2000).
    
    She, C.-Y. Spectral structure of laser light scattering revisited: 
    bandwidths of nonresonant scattering lidars. 
    Appl. Opt. 40, 4875-4884 (2001).
    
    Parameters
    ----------
    wavenumber : float
       wavenumber (defined as 1/lamda) in m-1 
    """

    Fk = kings_factor_CO2()

    e = (Fk - 1) * 9 / 2.

    return e

def epsilon_H2O(wavenumber):
    """ Returns the epsilon parameters of H2O for a given wavelength.
    
    epsilon = (gamma / alpha )**2 is the relative ansiotropy.
    
    Ra according to  She (2000).
    
    She, C.-Y. Spectral structure of laser light scattering revisited: 
    bandwidths of nonresonant scattering lidars. 
    Appl. Opt. 40, 4875-4884 (2001).
    
    Parameters
    ----------
    wavenumber : float
       wavenumber (defined as 1/lamda) in m-1 
    """

    Fk = kings_factor_H2O()

    e = (Fk - 1) * 9 / 2.

    return e


def epsilon_atmosphere(wavelength, C=385., p_e=0., p_t=1013.25):
    """ Calculate the combined epsilon for the atmosphere.
    
    epsilon = (gamma / alpha )**2 is the relative anisotropy. 
    
    Parameters
    ----------
    wavelength : float
       Light wavelength in nm
    C : float
       CO2 concentration [ppmv].
    p_e : float
       water-vapor pressure [hPa] 
    p_t : float
       Total atmospheric pressure pressure [hPa] 
       
    Returns
    -------
    epsilon:
       Relative anisotropy of the atmosphere
    """
    Fk = kings_factor_atmosphere(wavelength, C=C, p_e=p_e, p_t=p_t)
    epsilon = (Fk - 1.) * 9. / 2.
    return epsilon


# def gamma_square_N2(wavenumber, ignore_range=False):
#     """ Returns the gamma squared parameter for N2 for a given wavelength.
    
#     The empirical fit is take from:
    
#     Chance, K. V. & Spurr, R. J. D. Ring effect studies: Rayleigh scattering, 
#     including molecular parameters for rotational Raman scattering, and the 
#     Fraunhofer spectrum. Applied Optics 36, 5224 (1997).

#     The fit is valid for the 200nm - 1000nm range.
    
#     Parameters
#     ----------
#     wavenumber : float
#        wavenumber (defined as 1/lamda) in m-1
#     ignore_range: bool
#        If set to True, it will ignore the valid range of the formula (200nm - 1000nm)
    
#     Returns
#     -------
#     g : float
#        gamma squared in (F * m2)^2, F --> farad
#     """
#     wavelength = 1. / wavenumber * 10 ** 9 # in nm
#     if not ignore_range:
#         if (wavelength < 200.) or (wavelength > 1000):
#             raise ValueError("The empirical formula for gamma-squared is valid only between 200nm and 1000nm. Set ignore_range=False if you need to use it anyway.")

#     wavenumber_um = wavenumber * 10 ** -6

#     g = -6.01466 + 2385.57 / (186.099 - wavenumber_um ** 2)

#     g *= 10 ** -25  # Correct scale

#     return g ** 2  # Return gamma squared


# def gamma_square_O2(wavenumber, ignore_range=False):
#     """ Returns the gamma squared parameter for O2 for a given wavelength.
    
#     The empirical fit is take from:
    
#     Chance, K. V. & Spurr, R. J. D. Ring effect studies: Rayleigh scattering, 
#     including molecular parameters for rotational Raman scattering, and the 
#     Fraunhofer spectrum. Applied Optics 36, 5224 (1997).

#     The fit is valid for the 200nm - 1000nm range.

#     Parameters
#     ----------
#     wavenumber : float
#        wavenumber (defined as 1/lamda) in m-1 
#     ignore_range: bool
#        If set to True, it will ignore the valid range of the formula (200nm - 1000nm)

#     Returns
#     -------
#     g : float
#        gamma squared in (F * m2)^2, F --> farad
#     """
#     wavelength = 1. / wavenumber * 10 ** 9 # in nm
#     if not ignore_range:
#         if (wavelength < 200.) or (wavelength > 1000):
#             raise ValueError("The empirical formula for gamma-squared is valid only between 200nm and 1000nm. Set ignore_range=False if you need to use it anyway.")

#     wavenumber_um = wavenumber * 10 ** -6

#     g = 0.07149 + 45.9364 / (48.2716 - wavenumber_um ** 2)

#     g *= 10 ** -24  # Correct scale

#     return g ** 2  # Return gamma squared


def rho_atmosphere(wavelength, C=300., p_e=0., p_t=1013.25):
    """ Calculate the depolarization factor of the atmosphere. 
    
    Parameters
    ----------
    wavelength : float or array of floats
       Wavelength in nm
    C : float
       CO2 concentration in ppmv
    p_e : float
       water-vapor pressure [hPa]
    p_t : float
       total air pressure [hPa]

    Returnssa
    -------
    rho: float or array of floats
       Depolarization factor
    """
    F_k = kings_factor_atmosphere(wavelength, C=C, p_e=p_e, p_t=p_t)
    rho = (6 * F_k - 6) / (7 * F_k + 3)
    return rho

def delta_cabannes(epsilon):

    """ Calculate the depolarization ratio of the Cabannes line for one gas.
    Behrend and Nakamura 2002, eq 16
    
    Parameters
    ----------
    epsilon (scalar or numpy array):
       Relative anisotropy of the atmosphere: epsilon = (gamma / alpha )**2
       
    Returns
    -------
    delta_cabannes: 
        The linear depolarization ratio of the Cabannes line.
    """
    
    delta_cabannes = 3. * epsilon / (180. + 4. * epsilon)
    
    return(delta_cabannes)
    
def delta_cabannes_air(epsilon, gamma_sq, wv, conc, full = False):

    """ Calculate the depolarization ratio of the Cabannes line for one gas.
    Behrend and Nakamura 2002, combining eq. 6 and 7 and summing for all the gases
    
    Parameters
    ----------
    epsilon: list of scalars or 1D numpy arrays
       Relative anisotropy of the atmosphere: epsilon = (gamma / alpha )**2
    
    conc: list of scalars
       The relative concentration of each gas in the mixture (same size as epsilon)
       
    gamma_sq : list of scalars or 1Dnumpy arrays (same shape as epsilon)
       gamma squared in (F m2)^2, F --> farad
    
    wv : scalar or 1D numpy arrays (same size as epsilon)
       The incident wavelength in nm?
    
    Returns
    -------
    delta_air: 
        The linear depolarization ratio of the Cabannes line.
    """
    
    wv = np.array(wv) * 1E-9 # convert in m
    
    conc = np.array(conc)
    
    epsilon = np.array(epsilon)
    
    gamma_sq = np.array(gamma_sq)
    
    k = np.pi ** 2 / wv **4

    for i in range(len(conc)):
        
        alpha_sq = gamma_sq[i] / epsilon[i]
        
        cross_section_p = k * (alpha_sq + gamma_sq / 45.)

        cross_section_s = k * gamma_sq / 60.
        
    cross_section_p_cab = np.sum(cross_section_p * conc[:, np.newaxis], axis = 0)

    cross_section_s_cab = np.sum(cross_section_s * conc[:, np.newaxis], axis = 0)

    delta_cab = cross_section_s_cab / cross_section_p_cab

    
    if full:
        return(delta_cab, cross_section_p, cross_section_s)
    
    else:
        return delta_cab
    