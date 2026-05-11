#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:23:10 2026

@author: nikos
"""

def shallow_copy(input_data):
    
    output_data = input_data.copy()
    
    for k, v in input_data.items():
        if isinstance(v, dict):
            output_data[k] = v.copy()
            
    return output_data
