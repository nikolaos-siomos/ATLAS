#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings
import xarray as xr
import numpy as np
import netCDF4 as nc
from .readers.parse_ray_args import call_parser, check_parser
from .readers.check import check_channels
from .plotting import make_axis, make_title, make_plot
from .writters import make_header, export_ascii
from .tools import normalize, average, differentiate, curve_fit

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
    # Check the command line argument information
    args = check_parser(args)

    print('-----------------------------------------')
    print('Initializing the Rayleigh Fit...')
    print('-----------------------------------------')
    
    # Read the quicklook file
    data = xr.open_dataset(args['input_file'])

    # Extract signal
    sig = data.Range_Corrected_Signals
    sig = sig.copy().where(sig != nc.default_fillvals['f8'])
    
    # Extract Attenuated Backscatter
    atb = data.Attenuated_Backscatter
    
    # Extract IFF info
    dwl = data.Detected_Wavelength
    ewl = data.Emitted_Wavelength
    bdw = data.Channel_Bandwidth
            
    mol_method = data.Molecular_atmosphere_method
    if 'Sounding_Station_Name' in data.attrs:
        st_name = data.Sounding_Station_Name
    else: st_name = ''
    if 'Sounding_Start_Date' in data.attrs:
        rs_start_date = data.Sounding_Start_Date
    else: rs_start_date = ''
    if 'Sounding_Start_Time_UT' in data.attrs:
        rs_start_time = data.Sounding_Start_Time_UT
    else: rs_start_time = ''
    if 'WMO_Station_Number' in data.attrs:
        wmo_id = data.WMO_Station_Number
    else: wmo_id = ''
    if 'WBAN_Station_Number' in data.attrs:
        wban_id = data.WBAN_Station_Number
    else: wban_id = ''
        
    
    # Check if the parsed channels exist
    channels = \
        check_channels(sel_channels = args['channels'], 
                       all_channels = data.channel.values,
                       exclude_telescope_type = args['exclude_telescope_type'], 
                       exclude_channel_type = args['exclude_channel_type'], 
                       exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                       exclude_channel_subtype = args['exclude_channel_subtype'])
    
    # iterate over the channels
    for ch in channels:
        print(f"-- channel: {ch}")

        ch_d = dict(channel = ch)
        sig_ch = sig.loc[ch_d].values
        atb_ch = atb.loc[ch_d].values
        dwl_ch = dwl.copy().loc[ch_d].values
        ewl_ch = ewl.copy().loc[ch_d].values
        bdw_ch = bdw.copy().loc[ch_d].values
        
        # Create the y axis (height/range)
        x_lbin, x_ubin, x_llim, x_ulim, x_vals, x_label = \
            make_axis.rayleigh_x(heights = data.Height_levels.loc[ch_d].values, 
                                 ranges = data.Range_levels.loc[ch_d].values,
                                 x_lims = args['x_lims'], 
                                 use_dis = args['use_range'])
    
        # Smoothing
        if args['smooth']:
            if not isinstance(args['smoothing_window'],list):
                from .tools.smoothing import sliding_average_1D_fast as smooth_1D
            else:
                from .tools.smoothing import sliding_average_1D as smooth_1D

            y_vals_sm, y_errs = \
                smooth_1D(y_vals = sig_ch.copy(), 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')
    
        else:
            y_vals_sm = sig_ch.copy()
            y_errs = np.nan * y_vals_sm
        
        
        nder, mder, cder, rerr, rsem, mshp, coef, mneg = \
            curve_fit.stats(y1 = sig_ch.copy(),
                            y2 = atb_ch.copy(), 
                            x  = x_vals)   
            
        norm_region, mmol, idx, ismol = \
            curve_fit.scan(mder = mder, cder = cder,
                           rsem = rsem, mshp = mshp, 
                           mneg = mneg,
                           dflt_region = args['normalization_region'],
                           auto_fit = args['auto_fit'])
                
        nder_c = float(nder[idx].values)
        mder_c = float(mder[idx].values)
        rsem_c = float(rsem[idx].values)
        coef_c = float(coef[idx].values)
        
        # rsem_c = coef_c[0] * sem_c / atb_c
        
        # Create the y axis (signal)
        y_llim, y_ulim, y_label = \
            make_axis.rayleigh_y(sig = coef_c * y_vals_sm[slice(x_lbin,x_ubin+1)], 
                                 atb = atb_ch.copy()[slice(x_lbin,x_ubin+1)], 
                                 y_lims = args['y_lims'],
                                 wave = dwl_ch,
                                 use_lin = args['use_lin_scale'])
        
        # Make title
        title = make_title.rayleigh(start_date = data.RawData_Start_Date,
                                    start_time = data.RawData_Start_Time_UT, 
                                    end_time = data.RawData_Stop_Time_UT, 
                                    lidar = data.Lidar_Name, 
                                    channel = ch, 
                                    zan = data.Laser_Pointing_Angle,
                                    loc = data.Station_Name,
                                    dwl = dwl_ch,
                                    ewl = ewl_ch,
                                    bdw = bdw_ch,
                                    smooth = args['smooth'],
                                    sm_lims = args['smoothing_range'],
                                    sm_win = args['smoothing_window'],
                                    sm_expo = args['smooth_exponential'],
                                    mol_method = mol_method,
                                    st_name = st_name,
                                    rs_start_date = rs_start_date,
                                    rs_start_time = rs_start_time, 
                                    wmo_id = wmo_id,
                                    wban_id = wban_id)
        
    
        # Make plot filename
        fname = f'{data.Measurement_ID}_{data.Lidar_Name}_ray_{ch}_ATLAS_{__version__}.png'
    
        # Make the png file
        fpath = make_plot.rayleigh(dir_out = args['output_folder'], 
                                   fname = fname, title = title,
                                   dpi_val = args['dpi'],
                                   color_reduction = args['color_reduction'],
                                   use_lin = args['use_lin_scale'],
                                   norm_region = norm_region,
                                   x_vals = x_vals, 
                                   y1_vals = y_vals_sm,
                                   y2_vals = atb_ch,
                                   y1_errs = y_errs,
                                   coef = coef_c,
                                   rsem = rsem_c,
                                   rslope = nder_c,
                                   pval = mder_c,
                                   x_lbin = x_lbin, x_ubin = x_ubin,
                                   x_llim = x_llim, x_ulim = x_ulim, 
                                   y_llim = y_llim, y_ulim = y_ulim, 
                                   x_label = x_label, y_label = y_label,
                                   x_tick = args['x_tick']) 
        
        # Make ascii file header
        header = \
            make_header.rayleigh(start_date = data.RawData_Start_Date, 
                                 start_time = data.RawData_Start_Time_UT, 
                                 start_time_sec = data.Raw_Data_Start_Time.values,
                                 stop_time_sec = data.Raw_Data_Stop_Time.values, 
                                 wave = dwl_ch, 
                                 lidar = data.Lidar_Name, 
                                 loc = data.Station_Name, 
                                 meas_id = data.Measurement_ID, 
                                 channel = ch, 
                                 norm_region = norm_region, 
                                 st_name = st_name, 
                                 rs_start_date = rs_start_date,
                                 rs_start_time = rs_start_time, 
                                 wmo_id = wmo_id, 
                                 wban_id = wban_id)
        
        # Make the ascii filename
        ascii_name = f'{data.Measurement_ID}_{data.Lidar_Name}_ray_{ch}_ATLAS_{__version__}.txt'

        # Export to ascii (Volker's format)        
        export_ascii.rayleigh(dir_out = args['output_folder'], 
                              fname = ascii_name, 
                              alt = x_vals, 
                              atb = atb_ch, 
                              rcs = sig_ch, 
                              header = header)
        
        fname_mask = f'{data.Measurement_ID}_{data.Lidar_Name}_ray_{ch}_mask_ATLAS_{__version__}.png'

        make_plot.rayleigh_mask(dir_out = args['output_folder'], 
                                fname = fname_mask, title = title,
                                dpi_val = args['dpi'],
                                color_reduction = args['color_reduction'],
                                mmol = mmol, mder = mder, cder = cder, 
                                mshp = mshp, mneg = mneg, rsem = rsem)
        
        # Add metadata to the quicklook plot
        from PIL import Image
        from PIL import PngImagePlugin
       
        METADATA = {"processing_software" : f"ATLAS_{data.version}",
                    "measurement_id" : f"{data.Measurement_ID}",
                    "channel" : f"{ch}",
                    "smooth" : f"{args['smooth']}",
                    "smoothing_exponential" : f"{args['smooth_exponential']}",
                    "smoothing_range" : f"{args['smoothing_range']}",
                    "smoothing_window": f"{args['smoothing_window']}",
                    "dpi" : f"{args['dpi']}",
                    "color_reduction" : f"{args['color_reduction']}",
                    "normalization_region" : f"{args['normalization_region']}",
                    "use_lin_scale" : f"{args['use_lin_scale']}",
                    "use_range" : f"{args['use_range']}",
                    "x_lims" : f"{args['x_lims']}",
                    "y_lims" : f"{args['y_lims']}",
                    "x_tick" : f"{args['x_tick']}"}

                
        im = Image.open(fpath)
        meta = PngImagePlugin.PngInfo()
    
        for x in METADATA.keys():
            meta.add_text(x, METADATA[x])
            
        im.save(fpath, "png", pnginfo = meta)
        
    print('-----------------------------------------')
    print(' ')
    return()

if __name__ == '__main__':
    # Get the command line argument information
    args = call_parser()

    # Call main
    main(args)
