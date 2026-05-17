#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:26:14 2026

@author: nikos
"""

def validate_selector(processor, selector):
    required = {'db', 'item_id'}
    missing = required - set(selector.keys())
    if missing:
        raise KeyError(f"Selector is missing keys: {missing}")

    db_name = selector['db']

    if not hasattr(processor, db_name):
        setattr(processor, db_name, {})  # create empty dict

def get_from_selector(processor, selector):
    validate_selector(processor, selector)
    container = getattr(processor, selector['db'])

    if selector['item_id'] not in container:
        raise KeyError(
            f"Item_id '{selector['item_id']}' not found in db '{selector['db']}'. "
            f"Available item_ids: {list(container.keys())}"
        )

    return container[selector['item_id']]

def initialize_selector(processor, selector):
    validate_selector(processor, selector)
    container = getattr(processor, selector['db'])
    container.setdefault(selector['item_id'], {})

def set_from_selector(processor, selector, value):
    validate_selector(processor, selector)
    container = getattr(processor, selector['db'])
    container.setdefault(selector['item_id'], value)

def get_source_map(io_id):
    source_map = {
        'profile':{'db':'profile_db', 'item_id':io_id},
        'profile_error':{'db':'profile_error_db', 'item_id':io_id},
        'profile_mask':{'db':'profile_mask_db', 'item_id':io_id},
        'meteo':{'db':'meteo_db', 'item_id':io_id},
        'molecular':{'db':'molecular_db', 'item_id':io_id},
        'background':{'db':'background_db', 'item_id':io_id},
        'background_error':{'db':'background_error_db', 'item_id':io_id},
        'background_mask':{'db':'background_mask_db', 'item_id':io_id},
        'time_mask':{'db':'time_mask_db', 'item_id':io_id},
        'channel_mask':{'db':'channel_mask_db', 'item_id':io_id},
        'range':{'db':'range_db', 'item_id':io_id},
        'height_agl':{'db':'height_agl_db', 'item_id':io_id},
        'height_asl':{'db':'height_asl_db', 'item_id':io_id},
        'system_info':{'db':'system_info_db', 'item_id':io_id},
        'time_info':{'db':'time_info_db', 'item_id':io_id},
        'channel_info':{'db':'channel_info_db', 'item_id':io_id},
        'pol_cal_info':{'db':'pol_cal_info_db', 'item_id':io_id},
        'radiosonde_info':{'db':'radiosonde_info_db', 'item_id':io_id},
        'gluing_info':{'db':'gluing_info_db', 'item_id':io_id},
        'shots':{'db':'shots_db', 'item_id':io_id},
        }

    return source_map
    
def get_io_map(io_id:str, sources:list = None):

    source_map = get_source_map(io_id)
    
    if sources is not None:
        invalid_keys = [key for key in sources if key not in source_map]
        
        if invalid_keys:
            raise KeyError(f"Provided keys are invalid: {invalid_keys}")
            
        io_map = {s: source_map[s] for s in sources if s in source_map}
    else:
        io_map = source_map
    
    return io_map