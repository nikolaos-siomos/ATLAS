"""
Retrieval of aerosol optical properties from elastic lidar signals.

*The following functions are obtained from:
https://gitlab.com/ioannis_binietoglou/lidar-processing/

The MIT License (MIT)

Copyright (c) 2015, Ioannis Binietoglou

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
from scipy.integrate import cumtrapz


def klett_backscatter_aerosol(range_corrected_signal, lidar_ratio_aerosol, beta_molecular, index_reference,
                              reference_range, sc_ratio_ref, bin_length, lidar_ratio_molecular=8.37758041):
    r"""Calculation of aerosol backscatter coefficient using Klett algorithm.

    The method also calculates aerosol backscatter above the reference altitude using forward integration approach.

    Parameters
    ----------
    range_corrected_signal : float.
       The range corrected signal.
    lidar_ratio_aerosol : float.
       The aerosol lidar ratio.
    beta_molecular : array_like
       The molecular backscatter coefficient. (m^-1 * sr^-1)
    index_reference : integer
       The index of the reference height. (bins)
    reference_range : integer
       The reference height range. (bins) * To calculate a mean value around the reference height
    sc_ratio_ref : float
       The scattering ratio (aer. bsc + mol. bsc / mol. bsc) value on the reference height.
    bin_length : float
       The vertical bin length. (m)
    lidar_ratio_molecular : float
       The molecular lidar ratio. Default value is :math:`8 \pi/3` which is a typical approximation.

    Returns
    -------
    beta_aerosol: float
       The aerosol backscatter coefficient. (m^-1 * sr^-1)

    Notes
    -----
    We estimate aerosol backscatter using the equation.

    .. math::
       \beta_{aer}(R) = \frac{A}{B-C} - \beta_{mol}(R)

    where

    .. math::
       A &= S(R) \cdot exp(-2\int_{R_{0}}^{R} [L_{aer}(r)-L_{mol}] \cdot \beta_{mol}(r) dr)

       B &= \frac{S(R_0)}{\beta_{aer}(R_{0})+\beta_{mol}(R_0)}

       C &= -2 \int_{R_0}^{R} L_{aer}(r) \cdot S(r) \cdot T(r, R_0) dr

    with

    .. math::
        T(r,R_0) = exp(-2\int_{R_0}^{r}[L_{aer}(r')-L_{mol}] \cdot \beta_{mol}(r') \cdot dr')

    and

    * :math:`R` the distance from the source,
    * :math:`R_0` the distance between the source and the reference region,
    * :math:`\beta_{aer}` the aerosol backscatter coefficient,
    * :math:`\beta_{mol}` the molecular backscatter coefficient,
    * :math:`S(R)` the range corrected signal,
    * :math:`P` the signal due to particle and molecular scattering,
    * :math:`L_{aer}` the aerosol lidar ratio (extinction-to-backscatter coefficient),
    * :math:`L_{mol}` the molecular lidar ratio.

    Note that `lidar_ratio_molecular` should correspond to the `beta_molecular` i.e. they should both correspond
    to total or Cabannes signal.

    References
    ----------
    Ansmann, A. and Muller, D.: Lidar and Atmospheric Aerosol Particles,
    in Lidar:  Range-Resolved Optical Remote Sensing of the Atmosphere, vol. 102,
    edited by C. Weitkamp, Springer, New York., 2005. p. 111.    
    """
    
    beta_aerosol_reference  = (sc_ratio_ref - 1.) * beta_molecular[index_reference]
    # Get molecular reference values
    beta_molecular_reference, range_corrected_signal_reference = \
        _get_reference_values(beta_molecular, index_reference,
                              range_corrected_signal, reference_range)
  
    range_corrected_signal[index_reference] = range_corrected_signal_reference    
        
    # Calculate the Tau-integral and Tau for each bin. Eq. 4.11 of Weitkamp
    numerator_integral_argument = (lidar_ratio_aerosol - lidar_ratio_molecular) * beta_molecular
    numerator_integral = _integrate_from_reference(numerator_integral_argument, index_reference, bin_length)
    tau = np.exp(-2. * numerator_integral)

    # Calculate the integral of the denominator
    denominator_integral_argument = lidar_ratio_aerosol * range_corrected_signal * tau
    denominator_integral = _integrate_from_reference(denominator_integral_argument, index_reference, bin_length)

    # Calculate the numerator and denominator
    numerator = range_corrected_signal * tau
    denominator = range_corrected_signal_reference / (
        beta_aerosol_reference + beta_molecular_reference) - 2. * denominator_integral

    # Sum of aerosol and molecular backscatter coefficients.
    beta_sum = numerator / denominator

    # Aerosol backscatter coefficient.
    beta_aerosol = beta_sum - beta_molecular

    return (beta_aerosol)


def _get_reference_values(beta_molecular, index_reference, range_corrected_signal, reference_range):
    """
    Determine the reference value for Klett retrieval.

    Parameters
    ----------
    beta_molecular : array_like
       The molecular backscatter coefficient. (m^-1 * sr^-1)
    index_reference : integer
       The index of the reference height. (bins)
    range_corrected_signal : float.
       The range corrected signal.
    reference_range : integer
       The reference height range. (bins)

    Returns
    -------
    beta_molecular_reference : float
       The reference molecular value
    range_corrected_signal_reference : float
       The reference value for the range corrected signal       
    """

    range_corrected_signal_reference = np.mean(
        range_corrected_signal[(index_reference - reference_range):(index_reference + reference_range)])

    beta_molecular_reference = beta_molecular[index_reference]
    
    return (beta_molecular_reference, range_corrected_signal_reference)


def _integrate_from_reference(integral_argument, index_reference, bin_length):
    """
    Calculate the cumulative integral the `integral_argument` from and below the reference point.

    Parameters
    ----------
    integral_argument : array_like
       The argument to integrate
    index_reference : integer
       The index of the reference height. (bins)
    bin_length : float
       The vertical bin length. (m)

    Returns
    -------
    tau_integral : array_like
       The cumulative integral from the reference point.
    """
    # Integrate from the reference point towards the beginning
    tau_integral_below = cumtrapz(integral_argument[:index_reference + 1][::-1], dx=-bin_length)[::-1]

    # Integrate from the reference point towards the end
    tau_integral_above = cumtrapz(integral_argument[index_reference:], dx=bin_length)

    # Join the arrays and set a 0 value for the reference point.
    tau_integral = np.concatenate([tau_integral_below, np.zeros(1), tau_integral_above])

    return (tau_integral)
