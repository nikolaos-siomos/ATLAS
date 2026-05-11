#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:42:17 2025

@author: nikos
"""

# Unit conversion functions
def km_asl_to_m_asl(x):
    return(1E3 * x)

def m_agl_to_m_asl(x, ground = 0.):
    return(x + ground)

def km_agl_to_m_asl(x, ground = 0.):
    return(1E3 * x + ground)

def geo_to_asl(x):
    Re = 6.371E6
    return(x * Re / (Re - x))

def Pa_to_hPa(x):
    return(1E-2 * x)

def hPa_to_Pa(x):
    return(100 * x)

def atm_to_hPa(x):
    return(x * 1013.25)

def atm_to_Pa(x):
    return(x * 101325)

def C_to_K(x):
    return(x + 273.16)

def Cx10_to_K(x):
    return(x/10. + 273.16)

def fraction_to_percent(x):
    return(100. * x)

def percent_to_fraction(x):
    return(x / 100.)
