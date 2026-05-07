#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 15:49:30 2022

@author: nick
"""

from .molecular_properties import epsilon_N2, epsilon_O2, epsilon_Ar, epsilon_CO2, epsilon_H2O

from .refractive_index import n_N2, n_O2, n_Ar, n_CO2, n_H2O

import numpy as np

# def N2(wavelength, relative_concentration = 0.78084, alpha_method = 'combined'):
def N2(wavelength, pressure = 1013.25, temperature = 293.15, relative_concentration = 0.79, alpha_method = 'combined'):
    
    wavenumber = 10 ** 9 / wavelength
    
    RI_N2, alpha, N_N2 = n_N2(wavelength = wavelength, pressure = pressure, temperature = temperature, method = alpha_method, extrapolate=True)
    alpha_square = alpha ** 2
    
    epsilon = epsilon_N2(wavenumber)
    

    gamma_square = alpha_square * epsilon

    N2_parameters = {'name': "N_{2}",
                     'molecule_type': 'linear',
                     'B0': 1.989570 * 1E2, # from Weitkamp 2005, in m-1
                     'D0': 5.76 * 1E-4, # from Behrend 2002, in m-1
                     'I': 1, # from Behrend 2002
                     'alpha_square': alpha_square,  
                     'gamma_square': gamma_square,  
                     'epsilon': epsilon,  
                     'g': [6, 3], # from Behrend 2002 for even/odd J - alowed molecular spin values depending on J
                     'relative_concentration': relative_concentration,
                     'refractive_index': RI_N2,
                     'number_concentration': N_N2,
                     'pressure': pressure,
                     'temperature': temperature}   
    
    return(N2_parameters)

# def O2(wavelength, relative_concentration = 0.20946, alpha_method = 'combined'):
def O2(wavelength, pressure = 1013.25, temperature = 293.15, relative_concentration = 0.21, alpha_method = 'combined'):
    
    wavenumber = 10 ** 9 / wavelength
    
    RI_O2, alpha, N_O2 = n_O2(wavelength = wavelength, pressure = pressure, temperature = temperature, method = alpha_method, extrapolate=True)
    alpha_square = alpha ** 2
    
    epsilon = epsilon_O2(wavenumber)
    
    gamma_square = alpha_square * epsilon

    O2_parameters = {'name': "O_{2}",
                     'molecule_type': 'linear',
                     'B0': 1.43768 * 1E2, # from Behrend 2002, in m-1
                     'D0': 4.85 * 1E-4, # from Behrend 2002, in m-1
                     'I': 0, # from Behrend 2002
                     'alpha_square': alpha_square, 
                     'gamma_square': gamma_square,  
                     'epsilon': epsilon, 
                     'g': [0, 1], # from Behrend 2002 for even/odd J
                     'relative_concentration': relative_concentration,
                     'refractive_index': RI_O2,
                     'number_concentration': N_O2,
                     'pressure': pressure,
                     'temperature': temperature}   
    
    return(O2_parameters)

# def Ar(wavelength, relative_concentration = 0.009340, alpha_method = 'combined'):
def Ar(wavelength, pressure = 1013.25, temperature = 293.15, relative_concentration = 0., alpha_method = 'combined'):
    
    wavenumber = 10 ** 9 / wavelength

    RI_Ar, alpha, N_Ar = n_Ar(wavelength = wavelength, pressure = pressure, temperature = temperature, method = alpha_method, extrapolate=True)
    alpha_square = alpha ** 2
    
    epsilon = epsilon_Ar(wavenumber)
    
    gamma_square = alpha_square * epsilon

    Ar_parameters = {'name': "Ar",
                     'molecule_type': 'spherical',
                     'I': 0, # because the proton and neutron number is even and the spins cancel each other out
                     'alpha_square': alpha_square,  
                     'gamma_square': gamma_square,  
                     'epsilon': epsilon,  
                     'g': [0, 1], # similar to O2 and CO2 based on the nuclear spin
                     'relative_concentration': relative_concentration,
                     'refractive_index': RI_Ar,
                     'number_concentration': RI_Ar,
                     'pressure': pressure,
                     'temperature': temperature}     
    
    return(Ar_parameters)

# def CO2(wavelength, relative_concentration = 0.00036, alpha_method = 'combined'):
def CO2(wavelength, pressure = 1013.25, temperature = 293.15, relative_concentration = 0., alpha_method = 'combined'):
    
    wavenumber = 10 ** 9 / wavelength
    
    RI_CO2, alpha, N_CO2 = n_CO2(wavelength = wavelength, pressure = pressure, temperature = temperature, method = alpha_method, extrapolate=True)
    alpha_square = alpha ** 2
    
    epsilon = epsilon_CO2(wavenumber)
    
    gamma_square = alpha_square * epsilon

    CO2_parameters = {'name': "CO_{2}",
                      'molecule_type': 'linear',
                      'B0': 0.39027 * 1E2, #from Sioris 2001 PhD thesis, in m-1
                      'D0': 1.27E-7 * 1E2, #from Sioris 2001 PhD thesis, in m-1
                      'I': 0, #from Sioris 2001 PhD thesis
                      'alpha_square': alpha_square,
                      'gamma_square': gamma_square,
                      'epsilon': epsilon, 
                      'g': [1, 0], #from Sioris 2001 PhD thesis for even/odd J
                      'relative_concentration': relative_concentration,
                      'refractive_index': RI_CO2,
                      'number_concentration': RI_CO2,
                      'pressure': pressure,
                      'temperature': temperature} 
    
    return(CO2_parameters)

def H2O(wavelength, pressure = 13.33, temperature = 293.15, relative_concentration = 0., alpha_method = 'combined'):
    
    wavenumber = 10 ** 9 / wavelength
    
    RI_H2O, alpha, N_H2O = n_H2O(wavelength = wavelength, pressure = pressure, temperature = temperature, method = alpha_method, extrapolate=True)
    alpha_square = alpha ** 2
    
    epsilon = epsilon_H2O(wavenumber)
    
    gamma_square = alpha_square * epsilon

    H2O_parameters = {'name': "H_{2}O",
                      'molecule_type': 'asymetric',
                      'alpha_square': alpha_square,  
                      'gamma_square': gamma_square,  
                      'epsilon': epsilon, 
                      'relative_concentration': relative_concentration,
                      'refractive_index': RI_H2O,
                      'number_concentration': RI_H2O,
                      'pressure': pressure,
                      'temperature': temperature}  
    
    return(H2O_parameters)

def air(wavelength, gases = ['N2', 'O2', 'Ar', 'CO2', 'H2O'], concentrations = [0.78084, 0.20946, 0.009340, 0.00036, 0.]):
    
    relative_concentrations = np.array(concentrations) / np.sum(concentrations)
    
    all_parameters = []
    for i in range(len(gases)):
        if gases[i] == 'N2':
            all_parameters.append(N2(wavelength, relative_concentration = relative_concentrations[i]))
        elif gases[i] == 'O2':
            all_parameters.append(O2(wavelength, relative_concentration = relative_concentrations[i]))
        elif gases[i] == 'Ar':
            all_parameters.append(Ar(wavelength, relative_concentration = relative_concentrations[i]))
        elif gases[i] == 'CO2':
            all_parameters.append(CO2(wavelength, relative_concentration = relative_concentrations[i]))
        elif gases[i] == 'H2O':
            all_parameters.append(H2O(wavelength, relative_concentration = relative_concentrations[i]))
        else:
            raise ValueError('Gas {} not implemented. The currently supported gases are "N2", "O2", "Ar", "CO2", "H2O"'.format(gases[i]))
    
    return(all_parameters)
