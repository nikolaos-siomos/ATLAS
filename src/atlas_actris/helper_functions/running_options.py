#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 18:56:34 2023

@author: nikos
"""

def auto_set_process(process):
    
    processing = {'ray' : True,
                  'tlc' : True,
                  'pcb' : True,
                  'drk' : True}

    if 'ray' not in process: processing['ray'] = False
    if 'tlc' not in process: processing['tlc'] = False
    if 'pcb' not in process: processing['pcb'] = False
    if 'drk' not in process: processing['drk'] = False
    if 'off' in process: 
        processing['ray'] = False
        processing['tlc'] = False
        processing['pcb'] = False
        processing['drk'] = False

        
    return(processing)

def auto_set_reprocess(processing, quick_run):
    
    if quick_run == True:
        reprocess = {'converter' : False,
                     'preprocessor' : False,
                     'visualizer' : True}
    else:
        reprocess = {'converter' : True,
                     'preprocessor' : True,
                     'visualizer' : True}
    
    return(reprocess)

def auto_set_quicklook(process_qck):

    quicklook = {'ray' : True,
                 'tlc' : True,
                 'pcb' : True,
                 'drk' : True}
    
    
    if 'ray' not in process_qck: quicklook['ray'] = False
    if 'tlc' not in process_qck: quicklook['tlc'] = False
    if 'pcb' not in process_qck: quicklook['pcb'] = False
    if 'drk' not in process_qck: quicklook['drk'] = False
    if 'off' in process_qck: 
        quicklook['ray'] = False
        quicklook['tlc'] = False
        quicklook['pcb'] = False
        quicklook['drk'] = False
        
    return(quicklook)