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
           short_reader(fpath = '/home/nikos/Big_Data/Intercomparison/Analysis/TONI/203_246_953_20251003_2/netcdf/converter/gar_203_246_953_20251003_110416_ray_ATLAS_0.5.0_prepro.nc', 
                        exclude_telescope_type = '', 
                        exclude_channel_type = '', 
                        exclude_acquisition_mode = '', 
                        exclude_channel_subtype = '', 
                        use_channels = None)
           
# Calculate the Noise Scale Factor (NSF)
zone = [8000, 15000]

shots_d = shots_d.median().values
shots = shots.median().values

frames = signal.time.size
frames_d = signal_d.time.size

signal_d_m = signal_d.mean(dim='time')

signal_d_sm,_ = smoothing(sig = signal_d_m, 
                         smoothing_window = 500,
                         smoothing_sbin = 750,
                         smoothing_ebin = -1)

signal_m = signal.mean(dim='time')

signal_dc_m = signal_m - signal_d_sm

signal_d_m = signal_d_m - signal_d_sm

ΔV_dc = signal_dc_m[:,zone[0]:zone[1]].std(dim = 'bins') * np.sqrt(shots * frames)
Vm_dc = signal_dc_m[:,zone[0]:zone[1]].mean(dim = 'bins') 

ΔV_d = signal_d_m[:,zone[0]:zone[1]].std(dim = 'bins') * np.sqrt(shots_d * frames_d)

NSF = np.sqrt(np.power(ΔV_dc, 2) - np.power(ΔV_d, 2)) / np.sqrt(Vm_dc)

# SNF_m = np.mean()