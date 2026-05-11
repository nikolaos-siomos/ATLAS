#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 16:20:32 2026

@author: nikos
"""

import re
import os
import glob
import numpy as np
from typing import Union
from pathlib import Path
from datetime import datetime
from utils.error_classes import CustomWarning

WYOMING_RE = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{4})_wyoming_(?P<wmo_id>\d{5})\.dat$"
)

ECMWF_RE = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{4})_ecmwf_(?P<site>[A-Za-z0-9][A-Za-z0-9_-]*)\.nc$"
)

SCC_RE = re.compile(
    r"^rs_(?P<date>\d{8})(?P<site>[A-Za-z]+?)(?P<time>\d{2})\.nc$"
)

ASCII_RE = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{4}).*\.txt$"
)


def select_radiosonde_filename(
    target: object,
    folder: Union[str, Path],
    time_limit: int = 18,
    priority_time_limit: int = 3,
) -> str:
    """
    Find the best radiosonde file in `folder` using filename timestamps.

    Only files matching these filename formats are considered:
      - Wyoming: <yyyymmdd>_<hhmm>_wyoming_<wmo_id>
      - ECMWF:   <yyyymmdd>_<hhmm>_ecmwf_<site>.nc
      - SCC:     rs_<yyyymmdd><site><hh>.nc
      - ASCII:   <yyyymmdd>_<hhmm>*.txt

    Selection rule:
      1. Prefer the closest Wyoming file within +/- `priority_time_limit` hours.
      2. Otherwise, select the closest supported file within +/- `time_limit` hours.
      3. Ignore unrelated files and folders.

    Parameters
    ----------
    target:
        Target time, usually numpy.datetime64.
    folder:
        Folder containing radiosonde files.
    time_limit:
        Fallback time window in hours. Default: 18.
    priority_time_limit:
        Wyoming priority window in hours. Default: 3.

    Returns
    -------
    output, status:
        output is a dictionary with keys:
          - "radiosonde_file": selected file path
          - "radiosonde_time": datetime parsed from filename
          - "radiosonde_format": "wyoming", "ecmwf", "scc", or "ascii"

        status codes:
          - 0: file found
          - 1: supported files exist, but none within `time_limit`
          - 2: no supported files found

        If status is 1 or 2, output is an empty dictionary.
    """
    
    filenames = glob.glob(os.path.join(folder,'*'))
    
    valid_name_ind = []
    time_stamps = []
    delta_hours = []
    file_formats = []
    
    output = {}
    
    for i, filename in enumerate(filenames):
        name = os.path.basename(filename)
        
        wyoming_match = WYOMING_RE.fullmatch(name)
        ecmwf_match = ECMWF_RE.fullmatch(name)
        scc_match = SCC_RE.fullmatch(name)
        ascii_match = ASCII_RE.fullmatch(name)
        
        match = [wyoming_match, ecmwf_match, scc_match, ascii_match]
        
        if any(match) and os.path.isfile(filename):
            
            if wyoming_match:
                date = wyoming_match.group("date")
                time = wyoming_match.group("time")
                strip_format = "%Y%m%d%H%M"
                file_format = 'wyoming'
                valid = True
                
            elif ecmwf_match:
                date = ecmwf_match.group("date")
                time = ecmwf_match.group("time")
                strip_format = "%Y%m%d%H%M"
                file_format = 'ecmwf'
                valid = True

            elif scc_match:
                date = scc_match.group("date")
                time = scc_match.group("time")
                strip_format = "%Y%m%d%H"
                file_format = 'scc'
                valid = True

            elif ascii_match:
                date = ascii_match.group("date")
                time = ascii_match.group("time")
                strip_format = "%Y%m%d%H%M"
                file_format = 'ascii'
                valid = True

            else:
                valid = False

                # raise Exception("No valid radiosonde file was detected. Supported formats: wyoming, ecmwf, scc, ascii")
            
            if valid:
                valid_name_ind.append(i)
    
                dt = datetime.strptime(date + time, strip_format)
                
                time_stamps.append(np.datetime64(dt,'ns'))          
                
                delta_hours.append((np.datetime64(dt) - target) / np.timedelta64(1, "h"))
    
                file_formats.append(file_format)
    
    if len(valid_name_ind) > 0:
        file_formats = np.array(file_formats)
        delta_hours = np.array(delta_hours)
        valid_name_ind = np.array(valid_name_ind)
        
        mask_wyoming = (file_formats == 'wyoming') & \
            (np.abs(delta_hours) <= priority_time_limit)
            
        mask_time = (np.abs(delta_hours) <= time_limit)
        
        if mask_wyoming.any():
            delta_hours_wyoming = delta_hours[mask_wyoming]
            valid_name_ind_wyoming = valid_name_ind[mask_wyoming]
            min_time_ind = np.argmin(np.abs(delta_hours_wyoming))
            
            target_filename = filenames[valid_name_ind_wyoming[min_time_ind]]
            target_timestamp = time_stamps[valid_name_ind_wyoming[min_time_ind]]
            target_format = file_formats[valid_name_ind_wyoming[min_time_ind]]
            
            output = {"radiosonde_file" : target_filename,
                      "radiosonde_format" : target_format,
                      "radiosonde_time" : target_timestamp,
                }
            
            status = 0 
            
        elif mask_time.any():
            delta_hours_time = delta_hours[mask_time]
            valid_name_ind_time = valid_name_ind[mask_time]
            min_time_ind = np.argmin(np.abs(delta_hours_time))
            
            target_filename = filenames[valid_name_ind_time[min_time_ind]]
            target_timestamp = time_stamps[valid_name_ind_time[min_time_ind]]
            target_format = file_formats[valid_name_ind_time[min_time_ind]]

            output = {"radiosonde_file" : target_filename,
                      "radiosonde_format" : target_format,
                      "radiosonde_time" : target_timestamp,
                }
            
            status = 0  
            
        else:
            status = 1
            
            CustomWarning(f"No valid radiosonde file detected within {time_limit} hours from the middle of the measurement inside: {folder}")
            print()
            
    else:
        status = 2

        CustomWarning(f"No valid radiosonde file detected inside: {folder}")
        print()
        
    return output, status
    
