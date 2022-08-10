# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 16:31:29 2021

@author: Peristera
"""
import matplotlib as mpl

# Return the matplotlib params to its default values
def mpl_default():
    
    # keep the info about interactive option 
    interactivefig = mpl.rcParams['interactive']
    mpl.rcParams.update(mpl.rcParamsDefault)
    
    mpl.rcParams['interactive'] = interactivefig
    return()