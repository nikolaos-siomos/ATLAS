#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 21:21:21 2026

@author: nikos
"""

import os
import pandas as pd
import xarray as xr

from __get_T_P_profiles_from_cloudnet__ import _download_from_cloudnet  
from __get_T_P_profiles_from_wyoming__ import _download_wyoming
from utils.select_radiosonde import select_radiosonde_filename, select_manual_radiosonde_file
from utils.printouts import print_header, print_entry
from utils.error_classes import CustomWarning

def find_radiosonde(caller_info, metadata):
    
    print_header("Selecting radiosonde file")
    
    radiosonde_file = caller_info.get('radiosonde_file')
    radiosonde_folder = caller_info.get('radiosonde_folder')
    rsonde_station_wmo_id = caller_info['rsonde_station_wmo_id']
    rsonde_station_name = caller_info['rsonde_station_name']
    cloudnet_station_name = caller_info['cloudnet_station_name']
    
    qa_tests = metadata['time_info'].keys()
    
    allowed_tests = ['ray','ray_pcb']
    
    metadata['radiosonde_info'] = {}
    
    for key in qa_tests:
        
        if key in allowed_tests:
            
            station_id = metadata['system_info'][key].sel({'parameters':'station_id'})
            
            print_entry(key)
            print()
            
            time = metadata['time_info'][key].time.values
            mid_time_dt64 = time[0] + (time[-1] - time[0]) / 2.
        
            mid_stamp = pd.Timestamp(mid_time_dt64)
            mid_date = mid_stamp.strftime("%d.%m.%Y")
            mid_time = mid_stamp.strftime("%H:%M:%S")
            
            if radiosonde_file:
                if radiosonde_folder:
                    CustomWarning(
                        "radiosonde_file was provided manually. "
                        "radiosonde_folder will be ignored."
                    )
                    print()

                radiosonde_info, status = select_manual_radiosonde_file(
                    target=mid_time_dt64,
                    radiosonde_file=radiosonde_file,
                )

                if status == 0:
                    print(f"Radiosonde file manually provided: {radiosonde_info['radiosonde_file']}")
                    print()

                    caller_info['radiosonde_status'] = status
                    radiosonde_info['measurement_time'] = mid_time_dt64
                    radiosonde_info['radiosonde_status'] = status

                    metadata['radiosonde_info'][key] = xr.DataArray(
                        data=list(radiosonde_info.values()),
                        dims=["parameters"],
                        coords={"parameters": list(radiosonde_info.keys())},
                    )

                else:
                    caller_info['radiosonde_status'] = status
                    CustomWarning(
                        "The manually provided radiosonde file could not be used. "
                        "Computations which need molecular profiles will not be performed"
                    )
                    print()

                continue

            wyoming_file_downloaded = False
            cloudnet_file_downloaded = False
            
            print("Downloading will be attempted")
            print()
            
            if rsonde_station_wmo_id:
                dl_status = None
                
                try:
                    dl_status = _download_wyoming(
                        wmo_id = rsonde_station_wmo_id,
                        date = mid_date,
                        time_utc = mid_time,
                        save_dir = radiosonde_folder,
                        )
                    
                except Exception as e:
                    CustomWarning(f"Downloading radiosonde from Wyoming failed:\n{e}")
                    print()
                
                else:
                    if dl_status.ok:
                        wyoming_file_downloaded = True
                        print(f"Downloading status from Wyoming: Downloaded {os.path.basename(dl_status.path)}")
                        print()
                        
                    else:
                        CustomWarning(f"Downloading radiosonde from Wyoming failed: {dl_status.message}")
                        print()
                
            if cloudnet_station_name:
                rs_path = None

                try:
                    rs_path = _download_from_cloudnet(
                        station_name=cloudnet_station_name,
                        date=mid_date,
                        time_utc=mid_time,
                        save_dir=radiosonde_folder,
                    )
                
                except Exception as e:
                    CustomWarning(f"Downloading radiosonde from Cloudnet failed:\n{e}")
                    print()
                
                else:
                    if rs_path:
                        cloudnet_file_downloaded = True
                        print(
                            f"Downloading status from Cloudnet: "
                            f"Downloaded {os.path.basename(rs_path)}"
                        )
                        print()
                    else:
                        CustomWarning(
                            "Downloading radiosonde from Cloudnet failed: "
                            "no file path was returned."
                        )
                        print()
            
            identifiers = {
                'wyoming' : rsonde_station_wmo_id,
                'ecmwf': cloudnet_station_name,
                'scc': station_id,
                'ascii': rsonde_station_name
                }
            radiosonde_info, status = select_radiosonde_filename(
                mid_time_dt64, 
                folder = radiosonde_folder,
                identifiers = identifiers
                )
            
            if status == 0:
                print(f"Radiosonde file detected: {radiosonde_info['radiosonde_file']}")
                print()

            else:
                if wyoming_file_downloaded or cloudnet_file_downloaded:
                    CustomWarning("Downloaded radiosonde files were not selected. Please check the selection rules and file timestamps.")
                    print()
                
                CustomWarning("Radiosonde not found and could not be downloaded. Computations which need molecular profiles will not be performed")
             
            if status == 0:
                
                caller_info['radiosonde_status'] = status
                radiosonde_info['measurement_time'] = mid_time_dt64
                radiosonde_info['radiosonde_status'] = status
                
                metadata['radiosonde_info'][key] = xr.DataArray(
                    data=list(radiosonde_info.values()),
                    dims=["parameters"],
                    coords={"parameters": list(radiosonde_info.keys())},
                    )
        else:
            radiosonde_info = {'radiosonde_status':-1}
            
            metadata['radiosonde_info'][key] = xr.DataArray(
                data=list(radiosonde_info.values()),
                dims=["parameters"],
                coords={"parameters": list(radiosonde_info.keys())},
                )

                        
    return metadata
                