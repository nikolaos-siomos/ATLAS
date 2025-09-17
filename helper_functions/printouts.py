#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 18:07:53 2025

@author: nikos
"""

def print_header(text: str) -> None:
    
    count = getattr(print_header, "_count", 0) + 1
    print_header._count = count
    
    print(' ')    
    print('-----------------------------------------------')
    print(f'{count}) {text}')
    print('-----------------------------------------------')
    print(' ')    

def print_subsection(text: str) -> None:
    
    print('-----------------------------------------------')
    print(f'{text}')   
        
