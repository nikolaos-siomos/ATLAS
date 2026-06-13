"""
@author: Nikos Siomos
"""
import re
import numpy as np
import xarray as xr
from datetime import datetime
import io, os, glob, contextlib, pickle, shutil

from typing import Callable
from typing import Any, Dict, List, Tuple

from utils.error_classes import FileReaderError
from utils.printouts import print_header, print_subsection, endpoint

from readers.read_scc import read_dataset as reader_scc
from readers.read_licel_lazy import read_dataset as reader_licel
from readers.read_tamarin_dimfix import read_dataset as reader_tamarin
from readers.read_polly_xt import read_dataset as reader_polly_xt
from readers.read_licel_matlab import read_dataset as reader_licel_matlab
from readers.read_polly_xt_first import read_dataset as reader_polly_xt_first

# Example registry of formats -> reading functions
READERS: dict[str, Callable[[str], object]] = {
    "scc": reader_scc,
    "licel": reader_licel,
    "polly_xt": reader_polly_xt,
}

reader_menu: dict[str, Callable[[str], object]] = {
    "scc": reader_scc,
    "licel": reader_licel,
    "tamarin": reader_tamarin,
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

def downcast_float_safe(da: xr.DataArray, tol: float = 1e-6) -> xr.DataArray:
    """
    Downcast float64 DataArray to float32 when possible.

    For NumPy-backed arrays, a safety check is performed before downcasting.
    For Dask-backed arrays, the downcast is returned lazily, without computing
    the full array. This preserves the RAM benefit of the lazy reader.
    """

    if not (np.issubdtype(da.dtype, np.floating) and da.dtype == np.float64):
        return da

    # Dask-backed arrays expose ``chunks``. Performing a full allclose check here
    # would force a full read/compute of the raw signal, defeating lazy loading.
    if da.chunks is not None:
        return da.astype(np.float32)

    arr = da.values
    arr32 = arr.astype(np.float32)

    if np.allclose(arr, arr32.astype(np.float64), rtol=tol, atol=tol, equal_nan=True):
        return da.astype(np.float32)

    return da


def _raw_cache_root(d: Dict[str, Any]) -> str:
    """Return the folder used for temporary raw-lidar Zarr caches.

    The cache is always placed directly under caller_info["output_folder"] / "cache".
    """

    output_folder = d.get("output_folder")

    if output_folder is None:
        raise FileReaderError(
            "caller_info['output_folder'] is required for raw-lidar caching"
        )

    return os.path.join(output_folder, "cache")

def _safe_cache_key(meas_key: str) -> str:
    """Convert a measurement key into a filesystem-safe cache key."""

    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(meas_key))


def _cache_entry_dir(cache_root: str, meas_key: str) -> str:
    return os.path.join(cache_root, _safe_cache_key(meas_key))


def _cache_signal_path(cache_root: str, meas_key: str) -> str:
    return os.path.join(_cache_entry_dir(cache_root, meas_key), "data.zarr")


def _cache_metadata_path(cache_root: str, meas_key: str) -> str:
    return os.path.join(_cache_entry_dir(cache_root, meas_key), "metadata.pkl")


def _cache_entry_exists(cache_root: str, meas_key: str) -> bool:
    """Check whether one measurement cache entry is complete enough to read."""

    return (
        os.path.isdir(_cache_signal_path(cache_root, meas_key))
        and os.path.isfile(_cache_metadata_path(cache_root, meas_key))
    )


def _cache_is_full(cache_root: str, physical_paths: List[Tuple[str, str]]) -> bool:
    """Return True only if every physical measurement has a cache entry."""

    if len(physical_paths) == 0:
        return False

    return all(_cache_entry_exists(cache_root, meas_key) for meas_key, _ in physical_paths)


def _clear_raw_cache(cache_root: str) -> None:
    """Remove the complete raw cache folder and recreate it empty."""

    if os.path.isdir(cache_root):
        shutil.rmtree(cache_root)

    os.makedirs(cache_root, exist_ok=True)


def _write_cache_entry(
    cache_root: str,
    meas_key: str,
    system_info: Any,
    channel_info: Any,
    time_info: Any,
    sig_raw: xr.DataArray,
    shots: xr.DataArray,
    meas_type: str,
) -> None:
    """
    Write one measurement to cache.

    The signal and shots are stored in Zarr. If sig_raw is Dask-backed, this
    computes chunk-by-chunk and does not assemble the full raw signal in RAM.
    Small metadata objects are stored in pickle next to the Zarr store.
    """

    entry_dir = _cache_entry_dir(cache_root, meas_key)

    if os.path.isdir(entry_dir):
        shutil.rmtree(entry_dir)

    os.makedirs(entry_dir, exist_ok=True)

    ds = xr.Dataset(
        data_vars={
            "signal": sig_raw,
            "shots": shots,
        }
    )

    ds.to_zarr(_cache_signal_path(cache_root, meas_key), mode="w")

    metadata = {
        "system_info": system_info,
        "channel_info": channel_info,
        "time_info": time_info,
        "mtype": meas_type,
    }

    with open(_cache_metadata_path(cache_root, meas_key), "wb") as f:
        pickle.dump(metadata, f, protocol=pickle.HIGHEST_PROTOCOL)


def _read_cache_entry(cache_root: str, meas_key: str) -> Tuple[Any, Any, Any, xr.DataArray, xr.DataArray, str]:
    """Read one cached measurement lazily from Zarr plus its pickled metadata."""

    if not _cache_entry_exists(cache_root, meas_key):
        raise FileReaderError(f"Raw lidar cache entry is missing or incomplete for {meas_key}")

    ds = xr.open_zarr(_cache_signal_path(cache_root, meas_key))

    with open(_cache_metadata_path(cache_root, meas_key), "rb") as f:
        metadata = pickle.load(f)

    return (
        metadata["system_info"],
        metadata["channel_info"],
        metadata["time_info"],
        ds["signal"],
        ds["shots"].compute(),
        metadata["mtype"],
    )

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

    if station_id == "tam":
        return "tamarin"

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

    Cache behavior is controlled by caller_info["quick_run"]:

    A1. quick_run is False:
        Read raw files normally and write/replace the raw-lidar cache.

    B1. quick_run is True and the cache is complete:
        Read profiles and metadata directly from the cache. Do not write again.

    B2. quick_run is True and the cache is missing/incomplete:
        Read raw files normally and write a fresh cache.

    The cache stores signals and shots as Zarr, so reopening from cache remains
    lazy and RAM-safe. Metadata are stored as small pickle files.
    """

    print_header("Reading lidar signals...")

    file_format = d["raw_file_format"]
    quick_run = bool(d.get("quick_run", False))
    cache_root = _raw_cache_root(d)

    profiles: Dict[str, Any] = {}

    metadata: Dict[str, Any] = {
        "system_info": {},
        "channel_info": {},
        "time_info": {},
        "shots": {},
        "mtype": {},
    }

    physical_paths = _iter_physical_paths(d)

    if len(physical_paths) == 0:
        endpoint(1)

    read_from_cache = quick_run and _cache_is_full(cache_root, physical_paths)

    if read_from_cache:
        print(f"-- quick_run = True and raw cache is complete. Reading from cache:\n   {cache_root}")

        for meas_key, _ in physical_paths:
            print_subsection(f"{meas_key} dataset")

            system_info, channel_info, time_info, sig_raw, shots, meas_type = _read_cache_entry(
                cache_root=cache_root,
                meas_key=meas_key,
            )

            profiles[meas_key] = sig_raw
            metadata["system_info"][meas_key] = system_info
            metadata["channel_info"][meas_key] = channel_info
            metadata["time_info"][meas_key] = time_info
            metadata["shots"][meas_key] = shots
            metadata["mtype"][meas_key] = meas_type

        if profiles == {}:
            endpoint(1)

        return profiles, metadata

    if quick_run:
        print(f"-- quick_run = True but raw cache is missing or incomplete. Reading raw files and creating cache:\n   {cache_root}")
    else:
        print(f"-- quick_run = False. Reading raw files and replacing raw cache:\n   {cache_root}")

    # For normal runs or incomplete quick-run cache, rebuild the cache so that
    # future quick runs are consistent and complete.
    _clear_raw_cache(cache_root)

    for meas_key, filepath in physical_paths:
        meas_type = measurement_type(meas_key)

        print_subsection(f"{meas_key} dataset")

        system_info, channel_info, time_info, sig_raw, shots = reader_menu[file_format](
            filepath,
            meas_type=meas_type,
        )

        if isinstance(sig_raw, list):
            continue

        if chunk and sig_raw.chunks is None:
            time_chunks = min(5, max(10, sig_raw.sizes["time"]))
            channel_chunks = 5
            
            chunks = {
                "time": time_chunks,
                "channel": channel_chunks,
                "bins": -1,
            }

            sig_raw = sig_raw.chunk(chunks)

        if downscale:
            sig_raw = downcast_float_safe(sig_raw)

        profiles[meas_key] = sig_raw

        metadata["system_info"][meas_key] = system_info
        metadata["channel_info"][meas_key] = channel_info
        metadata["time_info"][meas_key] = time_info
        metadata["shots"][meas_key] = shots
        metadata["mtype"][meas_key] = meas_type

        _write_cache_entry(
            cache_root=cache_root,
            meas_key=meas_key,
            system_info=system_info,
            channel_info=channel_info,
            time_info=time_info,
            sig_raw=sig_raw,
            shots=shots,
            meas_type=meas_type,
        )

        # Reopen from cache so downstream computations use the fast chunked
        # temporary store rather than the original raw-file delayed tasks.
        system_info, channel_info, time_info, sig_raw_cached, shots_cached, meas_type = _read_cache_entry(
            cache_root=cache_root,
            meas_key=meas_key,
        )

        profiles[meas_key] = sig_raw_cached
        metadata["system_info"][meas_key] = system_info
        metadata["channel_info"][meas_key] = channel_info
        metadata["time_info"][meas_key] = time_info
        metadata["shots"][meas_key] = shots_cached
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
