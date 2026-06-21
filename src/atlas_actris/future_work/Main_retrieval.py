# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 11:36:27 2025

@author: Peristera
"""
import pathlib
import os
import glob
import xarray as xr
from products import klett_bsc

# parent_dir = os.getcwd()
# input_dir = os.path.join(parent_dir,'input')
# input_dir = r'C:\Users\Peristera\Documents\#NOA\Algorithms\DIVA_lidar_algorithm\input'
# input_file = glob.glob(os.path.join(input_dir,'*.nc'))[0]
input_file = r"C:\Users\Peristera\Documents\#NOA\Algorithms\DIVA_lidar_algorithm\input\mun_9_9_9_20251113_204513_ray_ATLAS_0.5.1_unofficial_prepro.nc"
filename = pathlib.Path(input_file).stem
filename_dir = os.path.split(input_file)[0]


sig_ds = xr.open_dataset(input_file, decode_times=False)


prod, height_levels, info_prod = [],[], [] # initialize prod and height_levels (xarray) and info_prod (dataframe)

# Klett Retrieval
# -- user INPUTS -- 
lr_const = 40. # sr
sr_ref = 1.03 # 1.02
ref_height = 5.5e3 # m
ref_hwindow = 500. # m


# Call the Klett function
prod, height_levels, info_prod = klett_bsc(prod, height_levels, info_prod, sig_ds, lr_const, sr_ref, ref_height, ref_hwindow)


# Export to netcdf 
export_ds = xr.Dataset() #prod.to_dataset(dim='product')
export_ds = export_ds.assign_coords({'product':prod.product.values,
                                     'bins':prod.bins.values,
                                     'time':sig_ds.time.values})

export_ds['level2_profiles'] = xr.Variable(['product','bins'], prod)

export_ds['height_levels'] = xr.Variable(['product','bins'], height_levels)

for col in info_prod.columns:
    export_ds[col] = xr.Variable(['product'], info_prod[col])

export_ds['time'] = sig_ds.time.copy()
export_ds['time'].values = sig_ds.time.values

# dataframe to xarray dataset
# info_ds = info_prod.to_xarray()
nc_name = f'{filename}_L2_profiles.nc'

parent_dir = os.path.split(filename_dir)[0]
output_dir = os.path.join(parent_dir,'output')
os.makedirs(output_dir, exist_ok=True)

# #Saves the netCDF file
export_ds.to_netcdf(os.path.join(output_dir,nc_name)) # format=nc_format with nc_format = 'NETCDF3_CLASSIC'
export_ds.close()
