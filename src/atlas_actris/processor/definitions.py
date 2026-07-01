#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 13:27:36 2026

@author: nikos
"""

profile_instances = [
    'profile', 
    'profile_mean', 
    'profile_low_res', 
    'profile_high_res'
    ]

profile_error_instances = [
    'profile_error', 
    'profile_error_mean', 
    'profile_error_low_res', 
    'profile_error_high_res'
    ]

background_instances = [
    'background', 
    'background_mean', 
    'background_low_res', 
    'background_high_res'
    ]

background_error_instances = [
    'background_error', 
    'background_error_mean', 
    'background_error_low_res', 
    'background_error_high_res'
    ]

profile_error_map = dict(zip(profile_instances,profile_error_instances))
background_map = dict(zip(profile_instances,background_instances))
background_error_map = dict(zip(profile_instances,background_error_instances))

assign_drk = {
    "ray": "drk_ray",
    "ray_pcb": "drk_ray_pcb",

    "pcb_p45": "drk_pcb",
    "pcb_m45": "drk_pcb",
    "pcb_aux_p45": "drk_pcb_aux",
    "pcb_aux_m45": "drk_pcb_aux",

    "tlc_north": "drk_tlc",
    "tlc_east": "drk_tlc",
    "tlc_south": "drk_tlc",
    "tlc_west": "drk_tlc",

    "tlc_inner": "drk_tlc_rin",
    "tlc_outer": "drk_tlc_rin",

    "trg": "drk_trg",
    "dtm_fi": "drk_dtm",
    "dtm_fo": "drk_dtm",
    "dtm": "drk_dtm",
    "cam": "drk_cam",
    }