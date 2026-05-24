"""
@author: Nikos Siomos
"""
import os, sys, glob, contextlib
import numpy as np
from readers.read_scc import read_dataset as reader_scc
from readers.read_licel import read_dataset as reader_licel
from readers.read_polly_xt import read_dataset as reader_polly_xt
from readers.read_licel_matlab import read_dataset as reader_licel_matlab
from readers.read_polly_xt_first import read_dataset as reader_polly_xt_first
import xarray as xr
import pandas as pd
from datetime import datetime
import re
from typing import Callable, Optional
from typing import Any, Dict, List, Tuple
from utils.printouts import print_header, print_subsection, endpoint
import contextlib, io
from utils.error_classes import FileReaderError

# Example registry of formats -> reading functions
READERS: dict[str, Callable[[str], object]] = {
    "scc": reader_scc,
    "licel": reader_licel,
    "polly_xt": reader_polly_xt,
}

reader_menu: dict[str, Callable[[str], object]] = {
    "scc": reader_scc,
    "licel": reader_licel,
    "polly_xt": reader_polly_xt,
    "licel_matlab": reader_licel_matlab,
    "polly_xt_first": reader_polly_xt_first,
}


def measurement_type(meas_key: str) -> str:
    """
    Infer the raw-reader measurement type from a caller_info measurement key.

    This mirrors the parser-side measurement_type helper and replaces the old
    flat mtype_* entries in caller_info. The returned value is still the value
    expected by the raw file readers.
    """

    if meas_key.startswith("drk"):
        return "drk"

    if meas_key.startswith("ray"):
        return "nrm"

    if meas_key in ["trg", "dtm", "nsf"]:
        return "nrm"

    if meas_key in ["tlc_north", "tlc_east", "tlc_south", "tlc_west"]:
        return "tlc_qua"

    if meas_key in ["tlc_inner", "tlc_outer"]:
        return "tlc_rin"

    if meas_key == "pcb_p45":
        return "pcb_p45"

    if meas_key == "pcb_m45":
        return "pcb_m45"

    if meas_key == "pcb_aux_p45":
        return "pcb_p45"

    if meas_key == "pcb_aux_m45":
        return "pcb_m45"

    if meas_key.startswith("cam"):
        return "cam"

    raise FileReaderError(f"Could not infer measurement type for {meas_key}")


def _iter_physical_paths(d: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Return physical measurement folders from caller_info['paths']."""

    paths = d.get("paths", {})

    if paths is None:
        paths = {}

    if not isinstance(paths, dict):
        raise FileReaderError("caller_info['paths'] must be a dictionary")

    return [
        (meas_key, filepath)
        for meas_key, filepath in paths.items()
        if filepath is not None and not meas_key.startswith("cam")
    ]

# def downcast_float_safe(da: xr.DataArray, tol=1e-6) -> xr.DataArray:
#     """
#     Downcast a float64 DataArray to float32 if conversion
#     error is within tolerance. Otherwise return unchanged.

#     Parameters
#     ----------
#     da : xr.DataArray
#         Input DataArray.
#     tol : float
#         Relative/absolute tolerance for allclose comparison.

#     Returns
#     -------
#     xr.DataArray
#         Downcasted DataArray (float32) if safe, else original.
#     """
#     if np.issubdtype(da.dtype, np.floating) and da.dtype == np.float64:
#         arr = da.values
#         arr32 = arr.astype(np.float32)
#         if np.allclose(arr, arr32.astype(np.float64), rtol=tol, atol=tol, equal_nan=True):
#             return da.astype(np.float32)
        
#     return da

def downcast_float_safe(da: xr.DataArray, tol: float = 1e-6) -> xr.DataArray:
    """
    Downcast float64 DataArray to float32 only if the full-array conversion
    error is within tolerance.

    For Dask-backed arrays, this performs one full lazy computation to check
    safety, but the returned astype operation remains lazy.
    
    Parameters
    ----------
    da : xr.DataArray
        Input DataArray.
    tol : float
        Relative/absolute tolerance for allclose comparison.

    Returns
    -------
    xr.DataArray
        Downcasted DataArray (float32) if safe, else original.
    """

    if not (np.issubdtype(da.dtype, np.floating) and da.dtype == np.float64):
        return da

    da32 = da.astype(np.float32)

    close = xr.apply_ufunc(
        np.isclose,
        da,
        da32.astype(np.float64),
        kwargs={
            "rtol": tol,
            "atol": tol,
            "equal_nan": True,
        },
        dask="allowed",
    )

    is_safe = bool(close.all().compute())

    if is_safe:
        return da32

    return da

def infer_format(d: Dict[str, Any], station_id: str, debug: bool = False) -> str:
    """
    Infer the raw file format by trying the registered readers on the physical
    folders listed in caller_info["paths"].

    The new parser stores only folders that should be read physically under
    caller_info["paths"].
    """

    print_header("Infering raw file format")

    # Exceptional readers are detected based on the station ID.
    if station_id == "evo":
        return "polly_xt_first"

    if station_id in ["brc", "run"]:
        return "licel_matlab"

    physical_paths = _iter_physical_paths(d)

    if len(physical_paths) == 0:
        endpoint(1)

    # For the rest of the stations, infer the reader by trial and error.
    for meas_key, filepath in physical_paths:
        meas_type = measurement_type(meas_key)

        print(f"Trying readers with {meas_key} dataset:")

        for fmt, reader in READERS.items():
            print(f"--{fmt}: ", end="")

            try:
                if debug:
                    with contextlib.redirect_stdout(io.StringIO()):
                        _, _, _, sig_raw, _ = reader(filepath, meas_type=meas_type)
                else:
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        _, _, _, sig_raw, _ = reader(filepath, meas_type=meas_type)

                if isinstance(sig_raw, list) and len(sig_raw) == 0:
                    print("Wrong reader")
                    continue

                print("Correct reader")
                print(f"File format: {fmt}")
                return fmt

            except Exception as e:
                print("Wrong reader")
                if debug:
                    print(e)
                continue

    raise FileReaderError("The raw QA file format is not supported")


def flexible_reader(
    d: Dict[str, Any],
    downscale: bool = True,
    chunk: bool = True,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, str]]:
    """
    Read only physical lidar signal folders from caller_info["paths"].

    The parser is responsible for deciding which logical measurements are
    aliases of an already-loaded physical folder. 
    """

    print_header("Reading lidar signals...")

    file_format = d["raw_file_format"]

    profiles: Dict[str, Any] = {}

    metadata: Dict[str, Any] = {
        "system_info": {},
        "channel_info": {},
        "time_info": {},
        "shots": {},
        "mtype": {},
    }
    
    physical_paths = _iter_physical_paths(d)

    for meas_key, filepath in physical_paths:
        meas_type = measurement_type(meas_key)

        print_subsection(f"{meas_key} dataset")

        system_info, channel_info, time_info, sig_raw, shots = reader_menu[file_format](
            filepath,
            meas_type=meas_type,
        )

        if isinstance(sig_raw, list):
            continue

        if chunk:
            time_chunks = min(50, max(10, sig_raw.sizes["time"]))
            bin_chunks = 4096

            profiles[meas_key] = sig_raw.chunk({
                "time": time_chunks,
                "channel": -1,
                "bins": bin_chunks,
            })
        else:
            profiles[meas_key] = sig_raw

        if downscale:
            profiles[meas_key] = downcast_float_safe(profiles[meas_key])

        metadata["system_info"][meas_key] = system_info
        metadata["channel_info"][meas_key] = channel_info
        metadata["time_info"][meas_key] = time_info
        metadata["shots"][meas_key] = shots
        metadata["mtype"][meas_key] = meas_type

    if profiles == {}:
        endpoint(1)

    return profiles, metadata
            
def radiosonde(finput_rs, delimiter, skip_header, skip_footer, 
               usecols, units, mtime, ground):

    """Extracts the meteorological information out of the 
    raw radiosonde file."""
    
    # Reading
    print('-----------------------------------------')
    print('Start reading radiosonde file...')
    print('-----------------------------------------')
    
    paths = glob.glob(os.path.join(finput_rs,'*_*.txt'))
    
    if len(paths) == 0:
        raise Exception(f"-- Error: No txt file was found in the radiosonde folder: {finput_rs} Please make sure that the radiosonde files are in txt format")
    
    lib_delimiter =  {"S": "",
                      "C": ",",
                      "T": "\t"}
    
    # Unit conversion functions
    def km_asl_to_m_asl(x):
        return(1E3 * x)

    def m_agl_to_m_asl(x, ground = 0.):
        return(x + ground)

    def km_agl_to_m_asl(x, ground = 0.):
        return(1E3 * x + ground)
    
    def geo_to_asl(x):
        Re = 6.371E6
        return(x * Re / (Re - x))
    
    def Pa_to_hPa(x):
        return(1E-2 * x)
    
    def atm_to_hPa(x):
        return(x * 1013.25)
    
    def C_to_K(x):
        return(x + 273.16)

    def Cx10_to_K(x):
        return(x/10. + 273.16)
    
    def fraction_to_percent(x):
        return(100. * x)
    
    if len(paths) == 0 :
        raise Exception("-- Error: No txt file provided in the radiosonde folder! Please provide a single file with the radiosonde data with a filename that starts with 'yyyymmdd_hhmm' and ends with '.txt' ")

    bname = [os.path.basename(path) for path in paths]

    bad_length = [len(name) < 14 for name in bname]
    
    if any(bad_length) :
        raise Exception(f"-- Error: Radiosonde filename with wrong length detected! Please revise the following files: {bname[bad_length]}. They should start with 'yyyymmdd_hhmm' and end with '.txt' ")
    else:
        pattern = "20[0-9]{2}[0-1][0-9][0-3][0-9]_[0-2][0-9][0-5][0-9]"
        bad_format = [not(bool(re.search(pattern,name))) for name in bname]

        if any(bad_format):
            raise Exception(f"-- Error: Radiosonde filename with wrong format detected! Please revise the following files: {bname[bad_format]}. They should start with 'yyyymmdd_hhmm' and end with '.txt' ")
        else:
            dates = [name[:13].split('_')[0] for name in bname]
            times = [name[:13].split('_')[1] for name in bname]

    bad_dates = [int(date[:4]) not in np.arange(1960,9999,1) or \
                 int(date[4:6]) not in np.arange(1,13,1) or \
                 int(date[6:8]) not in np.arange(1,32,1) for date in dates]

    bad_times = [int(time[:2]) not in np.arange(0,24,1) or \
                 int(time[2:4]) not in np.arange(0,60,1) for time in times]
        
    if any(bad_dates):
        raise Exception(f"-- Error: The date provided in at least one radiosonde filename is not correct. Please revise the following files: {np.array(bname)[bad_dates]}. It should start with 'yyyymmdd_hhmm' and end with '.txt' ")

    if any(bad_times):
        raise Exception(f"-- Error: The time provided in at least radiosond filename is not correct. Please revise the following files: {np.array(bname)[bad_times]}. It should start with 'yyyymmdd_hhmm' and end with '.txt' ")
        
    date_dt = np.array([datetime.strptime(date + time,'%Y%m%d%H%M') 
                        for date, time in zip(dates, times)])
    
    delta_t = np.array([(dt - mtime).total_seconds() /3600. for dt in date_dt])

    ind_rs = np.argmin(np.abs(delta_t))
    
    if not any(np.abs(delta_t) < 24):
        raise Exception(f"-- Error: The nearest radiosonde in time {bname[ind_rs]} was launched with a time difference of {np.round(delta_t[ind_rs],decimals=1)} hours with respect to the middle time of the measurement! Please provide a radiosond file with less than 18 hours temporal difference")
    else:
        print(f'-- Selected radiosonde file: {bname[ind_rs]}')
        
    if usecols[3] == None:
        parameters = ['P', 'T']
        usecols = usecols[:3]

    else:
        parameters = ['P', 'T', 'RH']
        
    data = np.genfromtxt(paths[ind_rs],skip_header = skip_header, 
                         skip_footer = skip_footer,
                         delimiter = lib_delimiter[delimiter], 
                         autostrip = True,
                         usecols = np.array(usecols) - 1, dtype = float)

    if units[0] == 'km_asl':
        data[:,0] = km_asl_to_m_asl(data[:,0])

    if units[0] in ['m_agl', 'km_agl']:
        if ground == None:
            raise Exception("-- Error: The altitude parameter of the rsonde_geodata field is mandatory when the radiosond height is in agl (altitude above ground units). Please provide at least 1 float corresponding to the station altitude: --rsonde_geodata 60.0")
        else:
            if units[0] == 'm_agl':
                data[:,0] = m_agl_to_m_asl(data[:,0], ground = ground)
            else:
                data[:,0] = km_agl_to_m_asl(data[:,0], ground = ground)

    # if 'geo' in units[0]:
    #     data[:,0] = geo_to_asl(data[:,0])
        
    if units[1] == 'Pa':
        data[:,1] = Pa_to_hPa(data[:,1])    

    if units[1] == 'atm':
        data[:,1] = atm_to_hPa(data[:,1])   
        
    if units[2] == 'C':
        data[:,2] = C_to_K(data[:,2]) 

    if units[2] == 'Cx10':
        data[:,2] = Cx10_to_K(data[:,2]) 
        
    if units[3] == 'fraction' and len(usecols) == 4:
        data[:,3] = fraction_to_percent(data[:,3])     
        
    alt = data[:,0]
    
    atmo = xr.DataArray(data[:,1:], 
                        coords = [alt, parameters], 
                        dims = ['height', 'parameters'] )
    
    return(dates[ind_rs], times[ind_rs], atmo)

def folder_to_sector(folder):

    fld = ['north','east','south','west','outer','inner']
    sec = [1,2,3,4,5,6]
            
    sector = np.nan * np.zeros(folder.shape)
    
    for i in range(len(fld)):
        sector[folder == fld[i]] = sec[i]
        
    return(sector)

def folder_to_position(folder):

    fld = ['static', '-45', '+45']
    sec = [0, 1, 2]
            
    position = np.nan * np.zeros(folder.shape)
    
    for i in range(len(fld)):
        position[folder == fld[i]] = sec[i]
        
    return(position)

def time_to_sector(folder, files_per_sector):
    
    blocks = folder.size / files_per_sector
    
    if blocks - np.floor(blocks) > 0.:
        raise Exception("-- Error: The files_per_sector argument was provided but " +
                 "the number of telecover files cannot be evenly divided by it! " +
                 "Please revise the files_per_sector value. If the number of " +
                 "files per sector was not constant during measurements then " +
                 "provide the telecover in individual folders per sector.")
    
    sec = [1, 2, 3, 4]
    
    sec_list = int(np.floor(blocks / 4.)) * sec
    sec_list.extend(sec[:((blocks - np.floor(blocks / 4.)) * 4.).astype(int)])
            
    sector = np.nan * np.zeros(folder.shape)
    
    for i in range(len(sec_list)):
        sector[i*files_per_sector:(i+1)*files_per_sector] = sec_list[i]
    
    return(sector)

def time_to_ring(folder, files_per_ring):
    
    blocks = folder.size / files_per_ring
    
    if blocks - np.floor(blocks) > 0.:
        raise Exception("-- Error: The files_per_ring argument was provided but " +
                 "the number of telecover files cannot be evenly divided by it! " +
                 "Please revise the files_per_ring value. If the number of " +
                 "files per ring was not constant during measurements then " +
                 "provide the telecover in individual folders per ring.")
    
    sec = [5, 6]
    
    sec_list = int(np.floor(blocks / 2.)) * sec
    sec_list.extend(sec[:((blocks - np.floor(blocks / 2.)) * 2.).astype(int)])
            
    ring = np.nan * np.zeros(folder.shape)
    
    for i in range(len(sec_list)):
        ring[i*files_per_ring:(i+1)*files_per_ring] = sec_list[i]
    
    return(ring)
