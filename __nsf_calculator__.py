#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:23:31 2024

@author: nikos
"""

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from processor.readers.read_files import short_reader
from processor.lidar_processing.signal import smoothing

system_info, channel_info, time_info, time_info_d,\
    signal, signal_d, shots, shots_d = \
           short_reader(fpath = '/home/nikos/Big_Data/WALLE/data/22092022/netcdf/converter/walle_20220922_141336_ray_ATLAS_0.4.6_prepro.nc', 
                        exclude_telescope_type = '', 
                        exclude_channel_type = '', 
                        exclude_acquisition_mode = '', 
                        exclude_channel_subtype = '', 
                        use_channels = None)
           
# Calculate the Noise Scale Factor (NSF)
zone = [12000, 15000]

shots_d = shots_d.median().values
shots = shots.median().values

signal_d_m = smoothing(sig = signal_d.mean(dim='time'), 
                       smoothing_window = 500,
                       smoothing_sbin = 750,
                       smoothing_ebin = -1)

signal_b = signal - signal_d_m

signal_d = signal_d - signal_d_m

ΔV_b = signal[:,:,zone[0]:zone[1]].std(dim = 'bins') * np.sqrt(shots)
Vm_b = signal[:,:,zone[0]:zone[1]].mean(dim = 'bins') 

ΔV_d = signal_d[0,:,zone[0]:zone[1]].std(dim = 'bins') * np.sqrt(shots_d)

NSF = np.sqrt(np.power(ΔV_b, 2) - np.power(ΔV_d, 2)) / np.sqrt(Vm_b)

# SNF_m = np.mean()