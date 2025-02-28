#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:32:07 2022

@author: nick
"""

import os, sys
import netCDF4 as nc
import numpy as np

def rayleigh_file(system_info, channel_info, time_info, time_info_d, nc_path,
                  meas_ID, sig, sig_d, shots, shots_d, 
                  P = None, T = None, radiosonde_file = None):

    print('-----------------------------------------')
    print('Start exporting to a Rayleigh QA file...')
    print('-----------------------------------------')
    
    """Creates the rayleigh netcdf file according to the SCC format 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html
    and exports it to nc_path"""
            
    n_time = sig.time.size
    n_channels = sig.channel.size
    n_points = sig.bins.size

    n_nb_of_time_scales = 1
    n_scan_angles = 1
    
    n_char_fname = len(time_info.filename.values[0])
    
    atlas_channel_label = channel_info.telescope_type.values + \
        channel_info.channel_type.values + \
            channel_info.acquisition_type.values + \
                channel_info.channel_subtype.values + \
                    np.round(channel_info.detected_wavelength.values.astype(float),
                             decimals = 0).astype('str')
                    
    start_time = [np.datetime64(t,'us').item() for t in time_info['start_time']]
    end_time = [np.datetime64(t,'us').item() for t in time_info['end_time']]
    
    start_t = [(dt - start_time[0]).seconds for dt in start_time]
    end_t = [(dt - start_time[0]).seconds for dt in end_time]
        
    Raw_Start_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    Raw_Stop_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    
    Raw_Start_Time[:,0] = start_t
    Raw_Stop_Time[:,0] = end_t
    
    mask_delay = channel_info.daq_trigger_offset.values >= 0.
    Trigger_Delay = channel_info.daq_trigger_offset.values.astype(float) * 150. / channel_info.range_resolution.values.astype(float)
    First_Signal_Bin = - channel_info.daq_trigger_offset.values.astype(int)
    Trigger_Delay[~mask_delay] = 0.
    First_Signal_Bin[mask_delay] = 0
    
    if not isinstance(sig_d,list):
        n_time_bck = sig_d.time.size
        start_time_d = [np.datetime64(t,'us').item() for t in time_info_d['start_time']] 
        end_time_d = [np.datetime64(t,'us').item() for t in time_info_d['end_time']]
        start_t_d = [(dt - start_time_d[0]).seconds for dt in start_time_d]
        end_t_d = [(dt - start_time_d[0]).seconds for dt in end_time_d]
        Bck_Start_Time = np.nan * np.zeros([n_time_bck,n_nb_of_time_scales])
        Bck_Stop_Time = np.nan * np.zeros([n_time_bck,n_nb_of_time_scales])
        Bck_Start_Time[:,0] = start_t_d
        Bck_Stop_Time[:,0] = end_t_d
        n_char_fname_d = len(time_info_d.filename.values[0])
    
    ds = nc.Dataset(nc_path,mode='w')

# Adding Dimensions    
    ds.createDimension('time', n_time)
    ds.createDimension('channels', n_channels)
    ds.createDimension('points', n_points)
    ds.createDimension('nb_of_time_scales', n_nb_of_time_scales)
    ds.createDimension('scan_angles', n_scan_angles)
    ds.createDimension('nchar_channel',4)
    ds.createDimension('nchar_filename',n_char_fname)
    if not isinstance(sig_d,list):  
        ds.createDimension('time_bck', n_time_bck)
        ds.createDimension('nchar_filename_d',n_char_fname_d)
    
# Adding Global Parameters
    ds.Altitude_meter_asl = float(system_info.altitude);

    ds.Latitude_degrees_north = float(system_info.latitude);

    ds.Longitude_degrees_east = float(system_info.longitude);
  
    ds.Lidar_Name = system_info.lidar_name;

    ds.Station_Name = system_info.station_name;

    ds.Station_ID = system_info.station_id.lower();

    if 'lidar_id' in system_info.index:
        ds.Lidar_ID = system_info.lidar_id;

    if 'version_name' in system_info.index:
        ds.Version_Name = system_info.version_name;

    if 'version_id' in system_info.index:
        ds.Version_ID = system_info.version_id;  
    
    if 'configuration_name' in system_info.index:
        ds.Configuration_Name = system_info.configuration_name;

    if 'configuration_id' in system_info.index:
        ds.Configuration_ID = system_info.configuration_id;

    ds.Measurement_ID = meas_ID;

    ds.Measurement_type = 'ray';

    if radiosonde_file:
       
       ds.Sounding_File_Name = radiosonde_file;
    
    if not isinstance(sig_d,list):

        ds.RawBck_Start_Date = start_time_d[0].strftime('%Y%m%d');

        ds.RawBck_Start_Time_UT = start_time_d[0].strftime('%H%M%S');

        ds.RawBck_Stop_Time_UT = end_time_d[-1].strftime('%H%M%S');

    ds.RawData_Start_Date = start_time[0].strftime('%Y%m%d');
    
    ds.RawData_Start_Time_UT = start_time[0].strftime('%H%M%S');
    
    ds.RawData_Stop_Time_UT = end_time[-1].strftime('%H%M%S');

# Adding Variables
    make_nc_var(ds, name = 'Acquisition_Mode', value = channel_info.acquisition_mode.values, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'ADC_resolution', value = channel_info.analog_to_digital_resolution.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_Low_Bin', value = channel_info.background_low_bin.values.astype(float), dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High_Bin', value = channel_info.background_high_bin.values.astype(float), dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Background_Low', value = channel_info.background_low.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High', value = channel_info.background_high.values, dtype = 'float', dims = ('channels',))
        
    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Background_Profile', value = sig_d.values, dtype = 'float', dims = ('time_bck', 'channels', 'points',))

    make_nc_var(ds, name = 'channel_ID', value = channel_info.scc_channel_id.values, dtype = 'int', dims = ('channels',))

    make_nc_str(ds, name = 'atlas_channel_label', value = atlas_channel_label, dims = ('channels','nchar_channel'), length = 4)    
    
    make_nc_var(ds, name = 'DAQ_Range', value = channel_info.data_acquisition_range.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Dead_Time', value = channel_info.dead_time.values, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Dead_Time_Correction_Type', value = channel_info.dead_time_correction_type.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Detected_Wavelength', value = channel_info.detected_wavelength.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Emitted_Wavelength', value = channel_info.emitted_wavelength.values, dtype = 'float', dims = ('channels',))

    make_nc_str(ds, name = 'Filename', value = time_info.filename.values, dims = ('time','nchar_filename'), length = n_char_fname)    

    if not isinstance(sig_d,list):
        make_nc_str(ds, name = 'Filename_Bck', value = time_info_d.filename.values, dims = ('time_bck','nchar_filename_d'), length = n_char_fname_d)    
            
    make_nc_var(ds, name = 'id_timescale', value = np.zeros(n_channels), dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'Channel_Bandwidth', value = channel_info.channel_bandwidth.values, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Laser_Pointing_Angle', value = np.array(n_scan_angles * [system_info.zenith_angle]), dtype = 'float', dims = ('scan_angles',))

    make_nc_var(ds, name = 'Laser_Pointing_Azimuth_Angle', value = np.array(n_scan_angles * [system_info.azimuth_angle]), dtype = 'float', dims = ('scan_angles',))
    
    make_nc_var(ds, name = 'Laser_Pointing_Angle_of_Profiles', value = np.zeros([n_time, n_nb_of_time_scales]), dtype = 'int', dims = ('time', 'nb_of_time_scales',))

    make_nc_var(ds, name = 'Laser_Repetition_Rate', value = channel_info.laser_repetition_rate.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Laser_Shots', value = shots.values, dtype = 'int', dims = ('time', 'channels',))
    
    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Background_Shots', value = shots_d.values, dtype = 'int', dims = ('time_bck', 'channels',))
    
    make_nc_var(ds, name = 'Raw_Lidar_Data', value = sig.values, dtype = 'float', dims = ('time', 'channels', 'points',))
    
    make_nc_var(ds, name = 'Raw_Data_Range_Resolution', value = channel_info.range_resolution.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Raw_Data_Start_Time', value = Raw_Start_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))
    
    make_nc_var(ds, name = 'Raw_Data_Stop_Time', value = Raw_Stop_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))

    make_nc_var(ds, name = 'Trigger_Delay', value = Trigger_Delay, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'First_Signal_Bin', value = First_Signal_Bin, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'DAQ_Trigger_Offset', value = channel_info.daq_trigger_offset.values, dtype = 'int', dims = ('channels',))
    
    if not isinstance(sig_d,list):
        
        make_nc_var(ds, name = 'Bck_Data_Start_Time', value = Bck_Start_Time, dtype = 'int', dims = ('time_bck', 'nb_of_time_scales',))
    
        make_nc_var(ds, name = 'Bck_Data_Stop_Time', value = Bck_Stop_Time, dtype = 'int', dims = ('time_bck', 'nb_of_time_scales',))
    
    if radiosonde_file:
    
        make_nc_var(ds, name = 'Molecular_Calc', value = 1, dtype = 'int')
    
    else:
        
        make_nc_var(ds, name = 'Molecular_Calc', value = 0, dtype = 'int')
        
        make_nc_var(ds, name = 'Pressure_at_Lidar_Station', value = P, dtype = 'float')
        
        make_nc_var(ds, name = 'Temperature_at_Lidar_Station', value = T, dtype = 'float')

    ds.close()
    
    return()

def telecover_file(system_info, channel_info, time_info, time_info_d, nc_path, 
                   meas_ID, sig, sig_d, shots, shots_d):

    print('-----------------------------------------')
    print('Start exporting to a Telecover QA file...')
    print('-----------------------------------------')

    """Creates the telecover netcdf file according to the SCC format 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html
    and exports it to nc_path"""
            
    n_time = sig.time.size
    n_channels = sig.channel.size
    n_points = sig.bins.size
    
    n_nb_of_time_scales = 1
    n_scan_angles = 1
    
    n_char_fname = len(time_info.filename.values[0])

    atlas_channel_label = channel_info.telescope_type.values + \
        channel_info.channel_type.values + \
            channel_info.acquisition_type.values + \
                channel_info.channel_subtype.values + \
                    np.round(channel_info.detected_wavelength.values.astype(float),
                             decimals = 0).astype('str')
    
    start_time = [np.datetime64(t,'us').item() for t in time_info['start_time']]
    end_time = [np.datetime64(t,'us').item() for t in time_info['end_time']]
    
    start_t = [(dt - start_time[0]).seconds for dt in start_time]
    end_t = [(dt - start_time[0]).seconds for dt in end_time]
        
    Raw_Start_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    Raw_Stop_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    
    Raw_Start_Time[:,0] = start_t
    Raw_Stop_Time[:,0] = end_t
    
    mask_delay = channel_info.daq_trigger_offset.values >= 0.
    Trigger_Delay = channel_info.daq_trigger_offset.values.astype(float) * 150. / channel_info.range_resolution.values.astype(float)
    First_Signal_Bin = - channel_info.daq_trigger_offset.values.astype(int)
    Trigger_Delay[~mask_delay] = 0.
    First_Signal_Bin[mask_delay] = 0
    
    if not isinstance(sig_d,list):
        n_time_bck = sig_d.time.size
        start_time_d = [np.datetime64(t,'us').item() for t in time_info_d['start_time']] 
        end_time_d = [np.datetime64(t,'us').item() for t in time_info_d['end_time']]
        start_t_d = [(dt - start_time_d[0]).seconds for dt in start_time_d]
        end_t_d = [(dt - start_time_d[0]).seconds for dt in end_time_d]
        Bck_Start_Time = np.nan * np.zeros([n_time_bck,n_nb_of_time_scales])
        Bck_Stop_Time = np.nan * np.zeros([n_time_bck,n_nb_of_time_scales])
        Bck_Start_Time[:,0] = start_t_d
        Bck_Stop_Time[:,0] = end_t_d
        n_char_fname_d = len(time_info_d.filename.values[0])

    ds = nc.Dataset(nc_path,mode='w')

# Adding Dimensions    
    ds.createDimension('time', n_time)
    ds.createDimension('channels', n_channels)
    ds.createDimension('points', n_points)
    ds.createDimension('nb_of_time_scales', n_nb_of_time_scales)
    ds.createDimension('scan_angles', n_scan_angles)
    ds.createDimension('nchar_channel',4)
    ds.createDimension('nchar_filename',n_char_fname)
    if not isinstance(sig_d,list):  
        ds.createDimension('time_bck', n_time_bck)
        ds.createDimension('nchar_filename_d',n_char_fname_d)
    
# Adding Global Parameters
    ds.Altitude_meter_asl = system_info.altitude;

    ds.Latitude_degrees_north = system_info.latitude;

    ds.Longitude_degrees_east = system_info.longitude;

    ds.Lidar_Name = system_info.lidar_name;

    ds.Station_Name = system_info.station_name;

    ds.Station_ID = system_info.station_id.lower();

    if 'lidar_id' in system_info.index:
        ds.Lidar_ID = system_info.lidar_id;

    if 'version_name' in system_info.index:
        ds.Version_Name = system_info.version_name;

    if 'version_id' in system_info.index:
        ds.Version_ID = system_info.version_id;  
    
    if 'configuration_name' in system_info.index:
        ds.Configuration_Name = system_info.configuration_name;

    if 'configuration_id' in system_info.index:
        ds.Configuration_ID = system_info.configuration_id;
  
    ds.Measurement_ID = meas_ID;

    ds.Measurement_type = 'tlc';
    
    if not isinstance(sig_d,list):

        ds.RawBck_Start_Date = start_time_d[0].strftime('%Y%m%d');

        ds.RawBck_Start_Time_UT = start_time_d[0].strftime('%H%M%S');

        ds.RawBck_Stop_Time_UT = end_time_d[-1].strftime('%H%M%S');

    ds.RawData_Start_Date = start_time[0].strftime('%Y%m%d');
    
    ds.RawData_Start_Time_UT = start_time[0].strftime('%H%M%S');
    
    ds.RawData_Stop_Time_UT = end_time[-1].strftime('%H%M%S');

# Adding Variables
    make_nc_var(ds, name = 'Acquisition_Mode', value = channel_info.acquisition_mode.values, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'ADC_resolution', value = channel_info.analog_to_digital_resolution.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_Low_Bin', value = channel_info.background_low_bin.values.astype(float), dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High_Bin', value = channel_info.background_high_bin.values.astype(float), dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Background_Low', value = channel_info.background_low.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High', value = channel_info.background_high.values, dtype = 'float', dims = ('channels',))
        
    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Background_Profile', value = sig_d.values, dtype = 'float', dims = ('time_bck', 'channels', 'points',))

    make_nc_var(ds, name = 'channel_ID', value = channel_info.scc_channel_id.values, dtype = 'int', dims = ('channels',))

    make_nc_str(ds, name = 'atlas_channel_label', value = atlas_channel_label, dims = ('channels','nchar_channel'), length = 4)    
    
    make_nc_var(ds, name = 'DAQ_Range', value = channel_info.data_acquisition_range.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Dead_Time', value = channel_info.dead_time.values, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Dead_Time_Correction_Type', value = channel_info.dead_time_correction_type.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Detected_Wavelength', value = channel_info.detected_wavelength.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Emitted_Wavelength', value = channel_info.emitted_wavelength.values, dtype = 'float', dims = ('channels',))

    make_nc_str(ds, name = 'Filename', value = time_info.filename.values, dims = ('time','nchar_filename'), length = n_char_fname)    

    if not isinstance(sig_d,list):
        make_nc_str(ds, name = 'Filename_Bck', value = time_info_d.filename.values, dims = ('time_bck','nchar_filename_d'), length = n_char_fname_d)    
            
    make_nc_var(ds, name = 'id_timescale', value = np.zeros(n_channels), dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'Channel_Bandwidth', value = channel_info.channel_bandwidth.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Laser_Pointing_Angle', value = np.array(n_scan_angles * [system_info.zenith_angle]), dtype = 'float', dims = ('scan_angles',))

    make_nc_var(ds, name = 'Laser_Pointing_Azimuth_Angle', value = np.array(n_scan_angles * [system_info.azimuth_angle]), dtype = 'float', dims = ('scan_angles',))

    make_nc_var(ds, name = 'Laser_Pointing_Angle_of_Profiles', value = np.zeros([n_time, n_nb_of_time_scales]), dtype = 'int', dims = ('time', 'nb_of_time_scales',))

    make_nc_var(ds, name = 'Laser_Repetition_Rate', value = channel_info.laser_repetition_rate.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Laser_Shots', value = shots.values, dtype = 'int', dims = ('time', 'channels',))

    make_nc_var(ds, name = 'Trigger_Delay', value = Trigger_Delay, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'First_Signal_Bin', value = First_Signal_Bin, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'DAQ_Trigger_Offset', value = channel_info.daq_trigger_offset.values, dtype = 'int', dims = ('channels',))
    
    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Background_Shots', value = shots_d.values, dtype = 'int', dims = ('time_bck', 'channels',))
    
    make_nc_var(ds, name = 'Raw_Lidar_Data', value = sig.values, dtype = 'float', dims = ('time', 'channels', 'points',))
    
    make_nc_var(ds, name = 'Raw_Data_Range_Resolution', value = channel_info.range_resolution.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Raw_Data_Start_Time', value = Raw_Start_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))
    
    make_nc_var(ds, name = 'Raw_Data_Stop_Time', value = Raw_Stop_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))
    
    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Bck_Data_Start_Time', value = Bck_Start_Time, dtype = 'int', dims = ('time_bck', 'nb_of_time_scales',))
    
        make_nc_var(ds, name = 'Bck_Data_Stop_Time', value = Bck_Stop_Time, dtype = 'int', dims = ('time_bck', 'nb_of_time_scales',))
    
    make_nc_var(ds, name = 'sector', value = time_info.sector.values, dtype = 'int', dims = ('time',))
    
    ds.close()
    
    return()

def polarization_calibration_file(
        system_info, channel_info, time_info, time_info_d, nc_path,
        meas_ID, sig, sig_d, shots, shots_d, molecular_calc = [], 
        P = [], T = [], radiosonde_file = None, rayleigh = []):

    print('-----------------------------------------')
    print('Start exporting to a Calibration QA file...')
    print('-----------------------------------------')
    
    """Creates the polarization calibration netcdf file according to the SCC format 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html
    and exports it to nc_path"""
        
    
    n_time = sig.time.size
    n_channels = sig.channel.size
    n_points = sig.bins.size
    
    n_nb_of_time_scales = 1
    n_scan_angles = 1

    n_char_fname = len(time_info.filename.values[0])

    atlas_channel_label = channel_info.telescope_type.values + channel_info.channel_type.values + channel_info.acquisition_type.values + channel_info.channel_subtype.values
    
    start_time = [np.datetime64(t,'us').item() for t in time_info['start_time']]
    end_time = [np.datetime64(t,'us').item() for t in time_info['end_time']]
    
    start_t = [(dt - start_time[0]).seconds for dt in start_time]
    end_t = [(dt - start_time[0]).seconds for dt in end_time]
        
    Raw_Start_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    Raw_Stop_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    
    Raw_Start_Time[:,0] = start_t
    Raw_Stop_Time[:,0] = end_t
    
    mask_delay = channel_info.daq_trigger_offset.values >= 0.
    Trigger_Delay = channel_info.daq_trigger_offset.values.astype(float) * 150. / channel_info.range_resolution.values.astype(float)
    First_Signal_Bin = - channel_info.daq_trigger_offset.values.astype(int)
    Trigger_Delay[~mask_delay] = 0.
    First_Signal_Bin[mask_delay] = 0
    
    if not isinstance(sig_d,list):
        n_time_bck = sig_d.time.size
        start_time_d = [np.datetime64(t,'us').item() for t in time_info_d['start_time']] 
        end_time_d = [np.datetime64(t,'us').item() for t in time_info_d['end_time']]
        start_t_d = [(dt - start_time_d[0]).seconds for dt in start_time_d]
        end_t_d = [(dt - start_time_d[0]).seconds for dt in end_time_d]
        Bck_Start_Time = np.nan * np.zeros([n_time_bck,n_nb_of_time_scales])
        Bck_Stop_Time = np.nan * np.zeros([n_time_bck,n_nb_of_time_scales])
        Bck_Start_Time[:,0] = start_t_d
        Bck_Stop_Time[:,0] = end_t_d
        n_char_fname_d = len(time_info_d.filename.values[0])

    ds = nc.Dataset(nc_path,mode='w')

# Adding Dimensions
    ds.createDimension('time', n_time)
    ds.createDimension('channels', n_channels)
    ds.createDimension('points', n_points)
    ds.createDimension('nb_of_time_scales', n_nb_of_time_scales)
    ds.createDimension('scan_angles', n_scan_angles)
    ds.createDimension('nchar_channel',4)
    ds.createDimension('nchar_filename',n_char_fname)
    if not isinstance(sig_d,list):  
        ds.createDimension('time_bck', n_time_bck)
        ds.createDimension('nchar_filename_d',n_char_fname_d)

# Adding Global Parameters
    ds.Altitude_meter_asl = system_info.altitude;

    ds.Latitude_degrees_north = system_info.latitude;

    ds.Longitude_degrees_east = system_info.longitude;

    ds.Lidar_Name = system_info.lidar_name;

    ds.Station_Name = system_info.station_name;

    ds.Station_ID = system_info.station_id.lower();

    if 'lidar_id' in system_info.index:
        ds.Lidar_ID = system_info.lidar_id;

    if 'version_name' in system_info.index:
        ds.Version_Name = system_info.version_name;

    if 'version_id' in system_info.index:
        ds.Version_ID = system_info.version_id;  
    
    if 'configuration_name' in system_info.index:
        ds.Configuration_Name = system_info.configuration_name;

    if 'configuration_id' in system_info.index:
        ds.Configuration_ID = system_info.configuration_id;
  
    ds.Measurement_ID = meas_ID;
    
    ds.Measurement_type = 'pcb';

    if radiosonde_file:
       
       ds.Sounding_File_Name = radiosonde_file;
    
    if not isinstance(sig_d,list):

        ds.RawBck_Start_Date = start_time_d[0].strftime('%Y%m%d');

        ds.RawBck_Start_Time_UT = start_time_d[0].strftime('%H%M%S');

        ds.RawBck_Stop_Time_UT = end_time_d[-1].strftime('%H%M%S');

    ds.RawData_Start_Date = start_time[0].strftime('%Y%m%d');
    
    ds.RawData_Start_Time_UT = start_time[0].strftime('%H%M%S');
    
    ds.RawData_Stop_Time_UT = end_time[-1].strftime('%H%M%S');

# Adding Variables
    make_nc_var(ds, name = 'Acquisition_Mode', value = channel_info.acquisition_mode.values, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'ADC_resolution', value = channel_info.analog_to_digital_resolution.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_Low_Bin', value = channel_info.background_low_bin.values.astype(float), dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High_Bin', value = channel_info.background_high_bin.values.astype(float), dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Background_Low', value = channel_info.background_low.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High', value = channel_info.background_high.values, dtype = 'float', dims = ('channels',))
        
    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Background_Profile', value = sig_d.values, dtype = 'float', dims = ('time_bck', 'channels', 'points',))

    make_nc_var(ds, name = 'channel_ID', value = channel_info.scc_channel_id.values, dtype = 'int', dims = ('channels',))

    make_nc_str(ds, name = 'atlas_channel_label', value = atlas_channel_label, dims = ('channels','nchar_channel'), length = 4)    
    
    make_nc_var(ds, name = 'DAQ_Range', value = channel_info.data_acquisition_range.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Dead_Time', value = channel_info.dead_time.values, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Dead_Time_Correction_Type', value = channel_info.dead_time_correction_type.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Detected_Wavelength', value = channel_info.detected_wavelength.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Emitted_Wavelength', value = channel_info.emitted_wavelength.values, dtype = 'float', dims = ('channels',))

    make_nc_str(ds, name = 'Filename', value = time_info.filename.values, dims = ('time','nchar_filename'), length = n_char_fname)    

    if not isinstance(sig_d,list):
        make_nc_str(ds, name = 'Filename_Bck', value = time_info_d.filename.values, dims = ('time_bck','nchar_filename_d'), length = n_char_fname_d)    
            
    make_nc_var(ds, name = 'id_timescale', value = np.zeros(n_channels), dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'Channel_Bandwidth', value = channel_info.channel_bandwidth.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Laser_Pointing_Angle', value = np.array(n_scan_angles * [system_info.zenith_angle]), dtype = 'float', dims = ('scan_angles',))

    make_nc_var(ds, name = 'Laser_Pointing_Azimuth_Angle', value = np.array(n_scan_angles * [system_info.azimuth_angle]), dtype = 'float', dims = ('scan_angles',))
  
    make_nc_var(ds, name = 'Laser_Pointing_Angle_of_Profiles', value = np.zeros([n_time, n_nb_of_time_scales]), dtype = 'int', dims = ('time', 'nb_of_time_scales',))

    make_nc_var(ds, name = 'Laser_Repetition_Rate', value = channel_info.laser_repetition_rate.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Laser_Shots', value = shots.values, dtype = 'int', dims = ('time', 'channels',))

    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Background_Shots', value = shots_d.values, dtype = 'int', dims = ('time_bck', 'channels',))
    
    make_nc_var(ds, name = 'Raw_Lidar_Data', value = sig.values, dtype = 'float', dims = ('time', 'channels', 'points',))
    
    make_nc_var(ds, name = 'Raw_Data_Range_Resolution', value = channel_info.range_resolution.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Raw_Data_Start_Time', value = Raw_Start_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))
    
    make_nc_var(ds, name = 'Raw_Data_Stop_Time', value = Raw_Stop_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))

    make_nc_var(ds, name = 'Trigger_Delay', value = Trigger_Delay, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'First_Signal_Bin', value = First_Signal_Bin, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'DAQ_Trigger_Offset', value = channel_info.daq_trigger_offset.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'calibrator_position', value = time_info.position.values, dtype = 'int', dims = ('time',))
    
    if not isinstance(sig_d,list):
        make_nc_var(ds, name = 'Bck_Data_Start_Time', value = Bck_Start_Time, dtype = 'int', dims = ('time_bck', 'nb_of_time_scales',))
    
        make_nc_var(ds, name = 'Bck_Data_Stop_Time', value = Bck_Stop_Time, dtype = 'int', dims = ('time_bck', 'nb_of_time_scales',))
    
    if radiosonde_file:
    
        make_nc_var(ds, name = 'Molecular_Calc', value = 1, dtype = 'int')
    
    else:
        
        make_nc_var(ds, name = 'Molecular_Calc', value = 0, dtype = 'int')
        
        make_nc_var(ds, name = 'Pressure_at_Lidar_Station', value = P, dtype = 'float')
        
        make_nc_var(ds, name = 'Temperature_at_Lidar_Station', value = T, dtype = 'float')

    if isinstance(rayleigh, list) == False:
        
        ds.Rayleigh_File_Name = rayleigh;
    
    ds.close()
    
    return()

def dark_file(system_info, channel_info, time_info, nc_path, 
              meas_ID, sig, shots):
    
    print('-----------------------------------------')
    print('Start exporting to a Standalone Dark QA file...')
    print('-----------------------------------------')

    """Creates the rayleigh netcdf file according to the SCC format 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html
    and exports it to nc_path"""
            
    n_time = sig.time.size
    n_channels = sig.channel.size
    n_points = sig.bins.size
    
    n_nb_of_time_scales = 1
    n_scan_angles = 1
    
    n_char_fname = len(time_info.filename.values[0])

    atlas_channel_label = channel_info.telescope_type.values + \
        channel_info.channel_type.values + \
            channel_info.acquisition_type.values + \
                channel_info.channel_subtype.values + \
                    np.round(channel_info.detected_wavelength.values.astype(float),
                             decimals = 0).astype('str')
                    
    start_time = [np.datetime64(t,'us').item() for t in time_info['start_time']]
    end_time = [np.datetime64(t,'us').item() for t in time_info['end_time']]
    
    start_t = [(dt - start_time[0]).seconds for dt in start_time]
    end_t = [(dt - start_time[0]).seconds for dt in end_time]
        
    Raw_Start_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    Raw_Stop_Time = np.nan * np.zeros([n_time,n_nb_of_time_scales])
    
    Raw_Start_Time[:,0] = start_t
    Raw_Stop_Time[:,0] = end_t
    
    mask_delay = channel_info.daq_trigger_offset.values >= 0.
    Trigger_Delay = channel_info.daq_trigger_offset.values.astype(float) * 150. / channel_info.range_resolution.values.astype(float)
    First_Signal_Bin = - channel_info.daq_trigger_offset.values.astype(int)
    Trigger_Delay[~mask_delay] = 0.
    First_Signal_Bin[mask_delay] = 0
    
    ds = nc.Dataset(nc_path,mode='w')

# Adding Dimensions
    ds.createDimension('time', n_time)
    ds.createDimension('channels', n_channels)
    ds.createDimension('points', n_points)
    ds.createDimension('nb_of_time_scales', n_nb_of_time_scales)
    ds.createDimension('scan_angles', n_scan_angles)
    ds.createDimension('nchar_channel',4)
    ds.createDimension('nchar_filename',n_char_fname)
        
# Adding Global Parameters
    ds.Altitude_meter_asl = system_info.altitude;

    ds.Latitude_degrees_north = system_info.latitude;

    ds.Longitude_degrees_east = system_info.longitude;

    ds.Lidar_Name = system_info.lidar_name;

    ds.Station_Name = system_info.station_name;

    ds.Station_ID = system_info.station_id.lower();

    if 'lidar_id' in system_info.index:
        ds.Lidar_ID = system_info.lidar_id;

    if 'version_name' in system_info.index:
        ds.Version_Name = system_info.version_name;

    if 'version_id' in system_info.index:
        ds.Version_ID = system_info.version_id;  
    
    if 'configuration_name' in system_info.index:
        ds.Configuration_Name = system_info.configuration_name;

    if 'configuration_id' in system_info.index:
        ds.Configuration_ID = system_info.configuration_id;
  
    ds.Measurement_ID = meas_ID;
    
    ds.Measurement_type = 'drk';

    ds.RawData_Start_Date = start_time[0].strftime('%Y%m%d');

    ds.RawData_Start_Time_UT = start_time[0].strftime('%H%M%S');

    ds.RawData_Stop_Time_UT = end_time[-1].strftime('%H%M%S');

# Adding Variables
    make_nc_var(ds, name = 'Acquisition_Mode', value = channel_info.acquisition_mode.values, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'ADC_resolution', value = channel_info.analog_to_digital_resolution.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_Low_Bin', value = channel_info.background_low_bin.values.astype(float), dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High_Bin', value = channel_info.background_high_bin.values.astype(float), dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Background_Low', value = channel_info.background_low.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Background_High', value = channel_info.background_high.values, dtype = 'float', dims = ('channels',))
        
    make_nc_var(ds, name = 'Raw_Lidar_Data', value = sig.values, dtype = 'float', dims = ('time', 'channels', 'points',))

    make_nc_var(ds, name = 'channel_ID', value = channel_info.scc_channel_id.values, dtype = 'int', dims = ('channels',))

    make_nc_str(ds, name = 'atlas_channel_label', value = atlas_channel_label, dims = ('channels','nchar_channel'), length = 4)    
    
    make_nc_var(ds, name = 'DAQ_Range', value = channel_info.data_acquisition_range.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Dead_Time', value = channel_info.dead_time.values, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Dead_Time_Correction_Type', value = channel_info.dead_time_correction_type.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Detected_Wavelength', value = channel_info.detected_wavelength.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Emitted_Wavelength', value = channel_info.emitted_wavelength.values, dtype = 'float', dims = ('channels',))

    make_nc_str(ds, name = 'Filename', value = time_info.filename.values, dims = ('time','nchar_filename'), length = n_char_fname)    
                    
    make_nc_var(ds, name = 'id_timescale', value = np.zeros(n_channels), dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'Channel_Bandwidth', value = channel_info.channel_bandwidth.values, dtype = 'float', dims = ('channels',))
    
    make_nc_var(ds, name = 'Laser_Pointing_Angle', value = np.array(n_scan_angles * [system_info.zenith_angle]), dtype = 'float', dims = ('scan_angles',))

    make_nc_var(ds, name = 'Laser_Pointing_Azimuth_Angle', value = np.array(n_scan_angles * [system_info.azimuth_angle]), dtype = 'float', dims = ('scan_angles',))

    make_nc_var(ds, name = 'Laser_Pointing_Angle_of_Profiles', value = np.zeros([n_time, n_nb_of_time_scales]), dtype = 'int', dims = ('time', 'nb_of_time_scales',))

    make_nc_var(ds, name = 'Laser_Repetition_Rate', value = channel_info.laser_repetition_rate.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Laser_Shots', value = shots.values, dtype = 'int', dims = ('time', 'channels',))
        
    make_nc_var(ds, name = 'Raw_Data_Range_Resolution', value = channel_info.range_resolution.values, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'Trigger_Delay', value = Trigger_Delay, dtype = 'float', dims = ('channels',))

    make_nc_var(ds, name = 'First_Signal_Bin', value = First_Signal_Bin, dtype = 'int', dims = ('channels',))

    make_nc_var(ds, name = 'DAQ_Trigger_Offset', value = channel_info.daq_trigger_offset.values, dtype = 'int', dims = ('channels',))
    
    make_nc_var(ds, name = 'Raw_Data_Start_Time', value = Raw_Start_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))
    
    make_nc_var(ds, name = 'Raw_Data_Stop_Time', value = Raw_Stop_Time, dtype = 'int', dims = ('time', 'nb_of_time_scales',))

    ds.close()
    
    return()

def radiosonde_file(nc_path, st_name, wmo_id, wban_id, 
                    date, time, ground, lat, lon, atmo):

    print('-----------------------------------------')
    print('Start exporting to a radiosonde QA file...')
    print('-----------------------------------------')
    
    """Creates the radiosonde netcdf file according to the SCC format 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html
    and exports it to nc_path"""

    n_points = atmo.height.size

    ds = nc.Dataset(nc_path,mode='w')

# Adding Dimensions
    ds.createDimension('points', n_points)
        
# Adding Global Parameters    
    ds.Altitude_meter_asl = ground;

    if lat != None:
        ds.Latitude_degrees_north = lat;

    if lon != None:
        ds.Longitude_degrees_east = lon;

    ds.Measurement_type = 'rs';

    ds.Sounding_Start_Date = date;

    ds.Sounding_Start_Time_UT = f'{time}00';

    if st_name != None:
        ds.Sounding_Station_Name = st_name;
        
    if wmo_id != None:
        ds.WMO_Station_Number = str(wmo_id);

    if wban_id != None:
        ds.WBAN_Station_Number = str(wban_id);

# Adding Variables
    make_nc_var(ds, name = 'Altitude', value = atmo.height.values, dtype = 'float', dims = ('points',))

    make_nc_var(ds, name = 'Pressure', value = atmo.loc[dict(parameters = 'P')].values, dtype = 'float', dims = ('points',))

    make_nc_var(ds, name = 'Temperature', value = atmo.loc[dict(parameters = 'T')].values, dtype = 'float', dims = ('points',))

    if 'RH' in atmo.parameters.values:
        make_nc_var(ds, name = 'RelativeHumidity', value = atmo.loc[dict(parameters = 'RH')].values, dtype = 'float', dims = ('points',))

    ds.close()
    
    return()

def debug_file(path, data, meas_type, label, meas_ID, show_index = False, header = True):
    
    debug_folder = os.path.join(path, 'debug')
    
    fname = os.path.join(debug_folder,f'{meas_type}_{meas_ID}_{label}.txt')
    
    data.to_csv(fname, index = show_index, header = header)
    
    return()

def meas_id(lr_id, time):
    """Creates a sting variable with the measurement ID variable according to: 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html"""
    
    start_date = time.dt.date.values[0].strftime('%Y%m%d')
    
    start_time = time.dt.time.values[0].strftime('%H%M')
        
    meas_ID = f"{start_date}{lr_id}{start_time}"
    
    return(meas_ID)
        
def path(output_folder, system_info, time, meas_type, version):
    """Creates a sting variable with the full path to the output QA file: 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html"""

    start_date = time.dt.date.values[0].strftime('%Y%m%d')

    start_time = time.dt.time.values[0].strftime('%H%M%S')

    station_id = system_info['station_id'].lower()

    try: lidar_id = system_info['lidar_id']
    except: lidar_id = ''
    
    try: version_id = system_info['version_id']
    except: version_id = ''
    
    try: config_id = system_info['configuration_id']
    except: config_id = ''
    
    parts = [str(station_id), str(lidar_id), str(version_id), str(config_id), str(start_date), str(start_time), meas_type, 'ATLAS', str(version), 'prepro']
    fname = "_".join([part for part in parts if len(part) > 0]) + ".nc"
    
    nc_path = os.path.join(output_folder,fname)
    
    if os.path.exists(output_folder) == False:
        os.makedirs(nc_path)
        
    if os.path.exists(nc_path):
        os.unlink(nc_path)
        
    return(nc_path)

def rs_path(output_folder, meas_ID):
    """Creates a sting variable with the full path to the output QA file: 
    https://docs.scc.imaa.cnr.it/en/latest/file_formats/netcdf_file.html"""
    
    nc_path = os.path.join(output_folder,f'rs_{meas_ID[:-4]}00.nc')
    
    if os.path.exists(output_folder) == False:
        os.makedirs(nc_path)
        
    if os.path.exists(nc_path):
        os.unlink(nc_path)
        
    return(nc_path)

def make_nc_var(ds, name, value, dtype, dims = []):  
    """Function called by the *_file functions in order to fascilitate variable
    creation in the netcdf"""
    
    if dtype == 'int':
        func = np.int32
        default_val = nc.default_fillvals['i4']
        
    if dtype == 'float':
        func = np.double
        default_val = nc.default_fillvals['f8']
      
    if len(dims) == 0:
        value = func(value)
    else:
        value[value != value] = default_val
        value = value.astype(dtype)

    var = ds.createVariable(name, func, dims)
    
    if len(dims) == 0:
        var[:] = value

    elif len(dims) == 1:
        var[:] = value

    elif len(dims) == 2:
        var[:,:] = value
        
    elif len(dims) == 3:
        var[:,:,:] = value
        
    else:
        sys.exit('-- Error: 4 or higher dimensional arrays not supported in function make_nc_var')
        
    return()

def make_nc_str(ds, name, value, dims, length):  
    """Function called by the *_file functions in order to fascilitate variable
    creation in the netcdf"""

    value_char = nc.stringtochar(value.astype(f'S{length}'))

    var = ds.createVariable(name, 'S1', dims)
    
    if len(dims) == 1:
        var[:] = value_char

    elif len(dims) == 2:
        var[:,:] = value_char
        
    elif len(dims) == 3:
        var[:,:,:] = value_char
        
    else:
        sys.exit('-- Error: 4 or higher dimensional arrays and scalars not supported in function make_nc_str')
        
    return()
    
    