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
import os
from PIL import Image
from PIL import PngImagePlugin
       
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
    
    dead_time = data.Dead_Time
    daq_trigger_offset = data.DAQ_Trigger_Offset
    background_low_bin = data.Background_Low_Bin
    background_high_bin = data.Background_High_Bin
    raw_data_range_resolution = data.Raw_Data_Range_Resolution
    
    # Extract SCC info
    station_id = data.Station_ID.lower()
        
    try: lidar_id = data.Lidar_ID
    except: lidar_id = ''
    
    try: version_id = data.Version_ID
    except: version_id = ''
    
    try: config_id = data.Configuration_ID
    except: config_id = ''
    
    try: config_name = data.Configuration_Name
    except: config_name = ''
    
    try: scc_id = data.channel_ID
    except: scc_id = ''

    # Extract date info
    start_date = data.RawData_Start_Date
    start_time = data.RawData_Start_Time_UT
    stop_time = data.RawData_Stop_Time_UT
    
    # Relative Stnadard error of the mean upper limit for the auto_fit
    rsem_lim = 0.02
            
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
        scc_id_ch = scc_id.copy().loc[ch_d].values
        
        dead_time_ch = dead_time.copy().loc[ch_d].values
        daq_trigger_offset_ch = daq_trigger_offset.loc[ch_d].values
        background_low_bin_ch = background_low_bin.loc[ch_d].values
        background_high_bin_ch = background_high_bin.loc[ch_d].values
        raw_data_range_resolution_ch = raw_data_range_resolution.loc[ch_d].values
        
        background_mode_ch, background_low_ch, background_high_ch, \
            first_signal_rangebin_ch, trigger_delay_ch =\
                atlas_to_scc_triggering(background_low_bin = background_low_bin_ch, 
                                        background_high_bin = background_high_bin_ch, 
                                        daq_trigger_offset = daq_trigger_offset_ch, 
                                        raw_data_range_resolution = raw_data_range_resolution_ch)
        
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

            atb_ch, y_errs_mol = \
                smooth_1D(y_vals = atb_ch.copy(), 
                          x_vals = x_vals,
                          x_sm_lims = args['smoothing_range'],
                          x_sm_win = args['smoothing_window'],
                          expo = args['smooth_exponential'],
                          err_type = 'std')
        
        else:
            y_vals_sm = sig_ch.copy()
            y_errs = np.nan * y_vals_sm   
        
        ulim = 16.
        min_win = 1.
        max_win = 4.
        
        if args['normalization_region'][0] < args['cross_check_lim']:
            raise Exception(f"The provided cross_check_lim {args['cross_check_lim']} is lower than the lower edge of the normalization region {args['normalization_region'][0]}. Please revise the rayleight_fit section of the settings file")
        if args['normalization_region'][0] > ulim:
            ulim = args['normalization_region'][1]
        if args['normalization_region'][1] - args['normalization_region'][0] < min_win:
            min_win = args['normalization_region'][1] - args['normalization_region'][0]
        if args['normalization_region'][1] - args['normalization_region'][0] > max_win:
            max_win = args['normalization_region'][1] - args['normalization_region'][0]
        
        # Check for a fit range
        rsem, nder, mfit, msem, mder, msec, mshp, mcrc, coef = \
            curve_fit.stats(y1 = sig_ch.copy(),
                            y2 = atb_ch.copy(), 
                            x  = x_vals,
                            min_win = 1.,
                            max_win = 4.,
                            step = 0.1,
                            llim = args['cross_check_lim'],
                            ulim = 16.,
                            rsem_lim = rsem_lim)   
            
        norm_region, idx, fit = \
            curve_fit.scan(mfit = mfit,
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
        title = make_title.rayleigh(start_date = start_date,
                                    start_time = start_time, 
                                    end_time = stop_time, 
                                    lidar = data.Lidar_Name, 
                                    channel = ch, 
                                    station_id = station_id,
                                    lidar_id = lidar_id,
                                    version_id = version_id,
                                    config_id = config_id,
                                    config_name = config_name,
                                    scc_id = scc_id_ch,
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
        parts = [str(station_id), str(lidar_id), str(version_id), str(config_id), str(start_date), str(start_time), 'ray', str(ch), str(scc_id_ch), 'ATLAS', str(__version__)]
        fname = "_".join([part for part in parts if len(part) > 0]) + ".png"
        # fname = f'{station_id}_{lidar_id}_{version_id}_{config_id}_{start_date}_{start_time}_ray_{ch}_{scc_id_ch}_ATLAS_{__version__}.png'
    
        # Make the png file
        fpath = make_plot.rayleigh(dir_out = os.path.join(args['output_folder'],'plots'), 
                                   fname = fname, title = title,
                                   dpi_val = args['dpi'],
                                   color_reduction = args['color_reduction'],
                                   use_lin = args['use_lin_scale'],
                                   norm_region = norm_region,
                                   x_vals = x_vals, 
                                   y1_vals = y_vals_sm,
                                   y2_vals = atb_ch,
                                   y1_errs = y_errs,
                                   y2_errs = np.nan * np.zeros(y_errs.shape),
                                   coef = coef_c,
                                   rsem = rsem_c,
                                   rslope = nder_c,
                                   pval = mder_c,
                                   rsem_lim = rsem_lim,
                                   fit = fit,
                                   auto_fit = args['auto_fit'],
                                   x_lbin = x_lbin, x_ubin = x_ubin,
                                   x_llim = x_llim, x_ulim = x_ulim, 
                                   y_llim = y_llim, y_ulim = y_ulim, 
                                   x_label = x_label, y_label = y_label,
                                   x_tick = args['x_tick'],
                                   label_1 = 'measured',
                                   label_2 = 'molecular') 
        
        # Make ascii file header
        header = \
            make_header.rayleigh(start_date = start_date,
                                 start_time = start_time, 
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
        parts = [str(station_id), str(lidar_id), str(version_id), str(config_id), str(start_date), str(start_time), 'ray', str(ch), str(scc_id_ch), 'ATLAS', str(__version__)]
        ascii_name = "_".join([part for part in parts if len(part) > 0]) + ".txt"
        # ascii_name = f'{station_id}_{lidar_id}_{version_id}_{config_id}_{start_date}_{start_time}_ray_{ch}_{scc_id_ch}_ATLAS_{__version__}.txt'

        # Export to ascii (Volker's format)        
        export_ascii.rayleigh(dir_out = args['output_folder'], 
                              fname = ascii_name, 
                              alt = x_vals, 
                              atb = atb_ch, 
                              rcs = sig_ch, 
                              header = header)
        
        if args['auto_fit']:
            parts = [str(station_id), str(lidar_id), str(version_id), str(config_id), str(start_date), str(start_time), 'ray', str(ch), str(scc_id_ch), 'mask', 'ATLAS', str(__version__)]
            fname_mask = "_".join([part for part in parts if len(part) > 0]) + ".png"
            # fname_mask = f'{station_id}_{lidar_id}_{version_id}_{config_id}_{start_date}_{start_time}_ray_{ch}_{scc_id_ch}_mask_ATLAS_{__version__}.png'
    
            # Make title
            title_mask = \
                make_title.rayleigh_mask(start_date = start_date,
                                         start_time = start_time, 
                                         end_time = stop_time, 
                                         lidar = data.Lidar_Name, 
                                         channel = ch, 
                                         station_id = station_id,
                                         lidar_id = lidar_id,
                                         version_id = version_id,
                                         config_id = config_id,
                                         config_name = config_name,
                                         scc_id = scc_id_ch,
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
            
            make_plot.rayleigh_mask(dir_out = os.path.join(args['output_folder'],'plots'), 
                                    fname = fname_mask, title = title_mask,
                                    dpi_val = args['dpi'],
                                    color_reduction = args['color_reduction'],
                                    mfit = mfit, mder = mder, msec = msec, 
                                    mshp = mshp, mcrc = mcrc, rsem = rsem,
                                    rsem_lim = rsem_lim)
        
        # Add metadata to the quicklook plot
        METADATA = {"processing_software" : f"ATLAS_{data.version}",
                    "measurement_type" : "ray",
                    "station_id" : f"{station_id}",
                    "lidar_id" : f"{lidar_id}",
                    "version_id" : f"{version_id}",
                    "config_id" : f"{config_id}",
                    "smooth" : f"{args['smooth']}",
                    "smoothing_exponential" : f"{args['smooth_exponential']}",
                    "smoothing_range" : f"{args['smoothing_range']}",
                    "smoothing_window": f"{args['smoothing_window']}",
                    "dpi" : f"{args['dpi']}",
                    "color_reduction" : f"{args['color_reduction']}",
                    "fit" : f"{fit}",
                    "auto_fit" : f"{args['auto_fit']}",
                    "rsem" : f"{rsem_c}",
                    "rslope" : f"{nder_c}",
                    "normalization_region" : f"{args['normalization_region']}",
                    "use_lin_scale" : f"{args['use_lin_scale']}",
                    "use_range" : f"{args['use_range']}",
                    "x_lims" : f"{args['x_lims']}",
                    "y_lims" : f"{args['y_lims']}",
                    "x_tick" : f"{args['x_tick']}",
                    "emission_wavelength" : f"{ewl_ch}",
                    "interference_filter_center" : f"{dwl_ch}",
                    "interference_filter_fwhm" : f"{bdw_ch}",
                    "atlas_channel_id" : f"{ch}",
                    "scc_channel_id" : f"{scc_id_ch}",
                    "dead_time" : f"{dead_time_ch}",
                    "daq_trigger_offset" : f"{daq_trigger_offset_ch}",
                    "first_signal_rangebin" : f"{first_signal_rangebin_ch}",
                    "trigger_delay" : f"{trigger_delay_ch}",
                    "background_low_bin" : f"{background_low_bin_ch}",
                    "background_high_bin" : f"{background_high_bin_ch}",
                    "background_mode" : f"{background_mode_ch}",
                    "background_low" : f"{background_low_ch}",
                    "background_high" : f"{background_high_ch}",
                    "maximum_channel_height" : f"{norm_region[-1]}"}

        im = Image.open(fpath)
        meta = PngImagePlugin.PngInfo()
    
        for x in METADATA.keys():
            meta.add_text(x, METADATA[x])
            
        im.save(fpath, "png", pnginfo = meta)
        
    print('-----------------------------------------')
    print(' ')
    return()

def atlas_to_scc_triggering(background_low_bin, background_high_bin,
                            daq_trigger_offset, raw_data_range_resolution):

    if daq_trigger_offset < -50:
        background_mode = 'Pre-Trigger'
        background_low = background_low_bin
        background_high = background_high_bin
        first_signal_rangebin = -daq_trigger_offset
        trigger_delay = -999.
    else:
        background_mode = 'Far Field'
        background_mode = raw_data_range_resolution * background_low_bin
        background_high = raw_data_range_resolution * background_high_bin
        if daq_trigger_offset <= 0:
            first_signal_rangebin = -daq_trigger_offset
            trigger_delay = -999.
        else:
            first_signal_rangebin = -999.
            trigger_delay = daq_trigger_offset * raw_data_range_resolution * 20. / 3.

    return(background_mode, background_low, background_high, first_signal_rangebin, trigger_delay)

if __name__ == '__main__':
    # Get the command line argument information
    args = call_parser()

    # Call main
    main(args)
