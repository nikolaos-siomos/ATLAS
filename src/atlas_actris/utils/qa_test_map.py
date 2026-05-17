#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 00:26:44 2026

@author: nikos
"""

qa_test_map = {
    'quicklook':['ray','tlc','pcb','drk','ray_pcb','pcb_aux','dtm','nsf'],
    'rayleigh_fit':['ray','ray_pcb'],
    'sector_telecover':['tlc_north','tlc_east','tlc_south','tlc_west'],
    'ring_telecover':['tlc_outer','tlc_inner'],
    'polarization_calibration':['pcb_p45','pcb_m45','ray_pcb'],
    'dark':['drk'],
    'zero_bin':['trg'],
    'filter_transmission':['pcb_p45','pcb_m45','pcb_aux_p45','pcb_aux_m45'],
    'dead_time':['dtm_fi','dtm_fo']
    }