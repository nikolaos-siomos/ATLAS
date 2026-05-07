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
        parser.read(path, encoding="utf-8")

# Intercomparison
        if parser.has_section('intercomparison'):

            self.cmp = read_section(parser['intercomparison'], dtype = object, squeeze = True)

        else:
            
            raise Exception("-- Error: No 'intercomparison' section is provided in the configuration files. Please include a section with at least the mandatory fields!")

# -------- END OF CLASS

def read_section(section, dtype=object, skip_vars=[], squeeze = False):
    # Reads the whole or part of the section and returns a Pandas Series
    map_info = dict()
    
    type_map = {}

    for key in section:
        if key not in skip_vars:
            arr = [i.strip() for i in re.split(',', section[key]) if i !='']
            if len(arr) > 0:
                arr = assume_type(arr)
                if squeeze and len(arr) == 1:
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
        elif '.' in arr[i] and arr[i].replace('.', '', 1).isdigit() == True:
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