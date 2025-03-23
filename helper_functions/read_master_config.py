"""
@author: P. Paschou & N. Siomos
"""
import configparser
import re
import numpy as np
import pandas as pd
import sys

class config():
    
    def __init__(self, path):
        """Reads the config file at the given path"""
        
        parser = configparser.ConfigParser()
        parser.optionxform = str
        parser.read(path, encoding="utf-8")

# # General
#         if parser.has_section('general'):

#             self.gen = read_section(parser['general'], dtype = object, squeeze = True)

#         else:
            
#             raise Exception("-- Error: No 'general' section is provided in the configuration files. Please include a section with at least the mandatory fields!")

# Converter
        if parser.has_section('converter'):

            self.cnv = read_section(parser['converter'], dtype = object, squeeze = True)

        else:
            
            self.cnv = dict()

# Processor
        if parser.has_section('preprocessor'):

            self.prs = read_section(parser['preprocessor'], dtype = object, squeeze = True)

        else:
            
            self.prs = dict()

# Quicklooks
        if parser.has_section('quicklooks'):

            self.qck = read_section(parser['quicklooks'], dtype = object, squeeze = True)

        else:
            
            self.qck = dict()

# Rayleigh Fit
        if parser.has_section('rayleigh_fit'):

            self.ray = read_section(parser['rayleigh_fit'], dtype = object, squeeze = True)

        else:
            self.ray = dict()

# Telecover
        if parser.has_section('telecover'):

            self.tlc = read_section(parser['telecover'], dtype = object, squeeze = True)

        else:
            
            self.tlc = dict()

# Polarization Calibration
        if parser.has_section('polarization_calibration'):

            self.pcb = read_section(parser['polarization_calibration'], dtype = object, squeeze = True, exception = ['ch_r','ch_t','K','G_R','G_T','H_R','H_T','R_to_T_transmission_ratio'])

        else:
            
            self.pcb = dict()
            
# Dark
        if parser.has_section('dark'):

            self.drk = read_section(parser['dark'], dtype = object, squeeze = True)

        else:
            
            self.drk = dict()


# -------- END OF CLASS

def read_section(section, dtype=object, skip_vars=[], squeeze = False, exception = []):
    # Reads the whole or part of the section and returns a Pandas Series
    map_info = dict()
    
    for key in section:
        if key not in skip_vars:
            arr = [i.strip() for i in re.split(',', section[key]) if i !='']
            if len(arr) > 0:
                arr = assume_type(arr)
                if squeeze and len(arr) == 1 and key not in exception:
                    map_info[key] = arr[0]
                else:
                    map_info[key] = arr
    
    return(map_info)

def assume_type(arr):
    
    for i in range(len(arr)):
        if arr[i] == 'True':
            arr[i] = True
        elif arr[i] == 'False':
            arr[i] = False
        elif '.' in arr[i] and arr[i].replace('.', '', 1).replace('-', '', 1).isdigit() == True:
            arr[i] = float(arr[i])
        elif arr[i].isdigit():
            arr[i] = int(arr[i])        
        
    return(arr)

def comma_split(var, dtype):
    
    if var != '':
        var = re.split(',', var)
    
        var = np.array([item.strip() for item in var], 
                       dtype = dtype) #trimming the spaces
    else:
        var=[]
    return(var)