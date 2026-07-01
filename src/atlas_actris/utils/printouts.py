#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 18:07:53 2025

@author: nikos
"""

import sys

termination_ids = {
    1: "No files to process! The QA test folders either don't exist or are empty.",
    2: "No channels to process! Check the provided channel_recorder_id values.",
    3: "Slicing and excluding measurement parts removed all measurements. No signals to process.",
    4: "Screening low shots removed all measurements. No signals to process.",
    5: "Handing overflows removed all measurements. No signals to process.",
    6: "Filtering options removed all channels. No signals to process.",
    }

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
    
def print_entry(text: str) -> None:
    
    print(f'-- {text}')
    
def counter():
    
    if not hasattr(counter, "count"):
        counter.count = 0
    counter.count += 1
    
    return counter.count
        
def endpoint(ID: str):
    
    text = termination_ids[ID]
    
    raise SystemExit(f"Termination ID {ID}: {text} ATLAS terminates here!")

    