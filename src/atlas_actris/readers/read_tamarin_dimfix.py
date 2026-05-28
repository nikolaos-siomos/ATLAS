"""
@authors: Siomos and Paschou

================
Input:
    arg 1: full filename from current running directory, should be a string scalar
    arg 2: the code (first letters) of the measurements of the lidar
Returns:
    arg 1-2: the signals if xarray format with dimensions (time, channel, bins)
    arg 3-4: info from header per channel (e.g. pmt type, laser polarization, bins, resolution, shots, channel type(p,s, total), ADC_range, ADC_bit, wavelength)
    arg 5-6: the ground altitude and the measurement angle(off-zenith)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import glob
import datetime as dt
import xarray as xr
from readers.check_file_format import detect_netcdf
from utils.error_classes import FileReaderError
from utils.time_conversions import datetimes_to_iso
from utils.error_classes import CustomWarning


# Read measurement
def read_dataset(dir_meas: str, meas_type: str):

    """
    dir_meas: Measurement folder (can contain more than one netcdf files)
    meas_type: If set to drk then the dark signals will be extracted
    """

    # Setting sig, info, and time as empty lists in the beginning
    sig_raw = []
    shots = []
    time_info = []

    system_info = []
    channel_info = []

    list_sig = []
    list_time = []
    list_shots = []

    input_path = Path(dir_meas)

    if not input_path.exists():
        CustomWarning("The folder for reading signals does not exist! " +
              f"Check the input directory! \n Given folder: {dir_meas}")
        return(system_info, channel_info, time_info, sig_raw, shots)

    if input_path.is_file():
        mfiles = [input_path]
    else:
        mfiles = [p for p in input_path.glob("*.*") if p.is_file()]

    if len(mfiles) == 0:
        CustomWarning(f"No files to read in: {dir_meas}")
        return(system_info, channel_info, time_info, sig_raw, shots)

    if detect_netcdf(mfiles[0]) == None:
        raise FileReaderError(f"--QA test folder contains non netcdf files: {dir_meas}")

    print(f'-- Reading {len(mfiles)} file(s)!')

    for k in range(len(mfiles)):

        if detect_netcdf(mfiles[k]) == None:
            raise FileReaderError(f"--The following file is not in netcdf format: {mfiles[k]}")

        if is_tamarin_file(mfiles[k]):
            time_info_f = get_time_info_tamarin(mfiles[k],
                                                meas_type = meas_type,
                                                filename = mfiles[k].name)
            if time_info_f.empty:
                continue

            if len(list_sig) == 0:
                system_info = read_meas_tamarin(mfiles[k], meas_type = meas_type)
                channel_info = read_channels_tamarin(mfiles[k], meas_type = meas_type)

            sig_raw_f = read_signals_tamarin(mfiles[k],
                                             time = time_info_f.index,
                                             channels = channel_info.index,
                                             meas_type = meas_type)

            is_analog = xr.DataArray(
                channel_info["acquisition_mode"].values == "a",
                coords={"channel": sig_raw_f.channel.values},
                dims=["channel"],
                )
            
            sig_raw_f = sig_raw_f.where(~is_analog, 500. - sig_raw_f)

            shots_f = read_shots_tamarin(mfiles[k],
                                         time = time_info_f.index,
                                         channels = channel_info.index,
                                         meas_type = meas_type)

        else:
            raw_data = xr.open_dataset(mfiles[k])

            if "Measurement_ID" not in raw_data.attrs:
                raise FileReaderError("Measurement_ID parameter not found in the netcdf file. This is not a scc raw file")

            # Reading the scc file metadata
            time_info_f = get_time_info(raw_data,
                                        meas_type = meas_type,
                                        filename = mfiles[k].name)
            if time_info_f.empty:
                continue

            if len(list_sig) == 0:
                system_info = read_meas(raw_data = raw_data)
                channel_info = read_channels(raw_data = raw_data)

            # Reading the licel signals
            sig_raw_f = read_signals(raw_data,
                                     time = time_info_f.index,
                                     channels = channel_info.index,
                                     meas_type = meas_type)

            shots_f = read_shots(raw_data,
                                 time = time_info_f.index,
                                 channels = channel_info.index,
                                 meas_type = meas_type)

        channel_info["bins"] = channel_info.index.size * [int(sig_raw_f.shape[-1] + 1)]

        # Append the arrays to list in order to concatenate later
        list_sig.append(sig_raw_f)
        list_shots.append(shots_f)
        list_time.append(time_info_f)

    if len(list_sig) > 0:
        # Append in the time dimension all the time frames
        sig_raw = xr.concat(list_sig, dim='time')
        shots = xr.concat(list_shots, dim='time')
        time_info = pd.concat(list_time)

        # Transpose the dims in [time, channel, bins]
        sig_raw = sig_raw.transpose('time','channel','bins')
        shots = shots.transpose('time','channel')

        # Sort by time
        sig_raw = sig_raw.sortby('time').copy()
        shots = shots.sortby('time').copy()
        time_info = time_info.sort_index()

    return(system_info, channel_info, time_info, sig_raw, shots)


def read_meas(raw_data):

    system_info = pd.Series(dtype = object)

    return(system_info)


def read_channels(raw_data):

    ch_index = raw_data["channel_ID"].values.astype(str)

    channel_info = pd.DataFrame(index = ch_index)
    channel_info["data_acquisition_range"] = raw_data["DAQ_Range"].values

    return(channel_info)


def get_time_info(raw_data, meas_type, filename):

    time_info = pd.Series(dtype = object)

    if meas_type == 'drk':
        if "Raw_Bck_Start_Time" not in list(raw_data.variables):
            print("--Raw_Bck_Start_Time parameter not found. No dark profile embedded in SCC file -> skipping")
            return(time_info)
        if "Raw_Bck_Stop_Time" not in list(raw_data.variables):
            print("--Raw_Bck_Stop_Time parameter not found. No dark profile embedded in SCC file -> skipping")
            return(time_info)
        sdate = raw_data.RawBck_Start_Date
        stime = raw_data.RawBck_Start_Time_UT
        start_time_sec = raw_data.Raw_Bck_Start_Time[:,0].values.astype(float)
        stop_time_sec = raw_data.Raw_Bck_Stop_Time[:,0].values.astype(float)
        filenames = np.empty(raw_data.time_bck.size, dtype = object)

    else:
        if "Raw_Data_Start_Time" not in list(raw_data.variables):
            print("--Raw_Data_Start_Time parameter not found. No profile other than dark embedded in SCC file -> skipping")
            return(time_info)
        if "Raw_Data_Stop_Time" not in list(raw_data.variables):
            print("--Raw_Data_Stop_Time parameter not found. No profile other than dark embedded in SCC file -> skipping")
            return(time_info)
        sdate = raw_data.RawData_Start_Date
        stime = raw_data.RawData_Start_Time_UT
        start_time_sec = raw_data.Raw_Data_Start_Time[:,0].values.astype(float)
        stop_time_sec = raw_data.Raw_Data_Stop_Time[:,0].values.astype(float)
        filenames = np.empty(raw_data.time.size, dtype = object)

    # Convert and store start time
    sdt = dt.datetime.strptime(sdate + ' ' + stime, "%Y%m%d %H%M%S") # start meas

    start_time_arr = np.array([sdt + dt.timedelta(seconds = t) for t in start_time_sec])
    end_time_arr = np.array([sdt + dt.timedelta(seconds = t) for t in stop_time_sec])

    filenames[:] = filename

    tdata = np.array([filenames,
                      datetimes_to_iso(start_time_arr),
                      datetimes_to_iso(end_time_arr)],
                     dtype = object)

    time_info = pd.DataFrame(tdata.T,
                             index = start_time_arr,
                             columns = ['filename', 'start_time', 'end_time'])

    return(time_info)


def read_signals(raw_data, time, channels, meas_type):

    if meas_type == 'drk':
        if "Background_Profile" not in list(raw_data.variables):
            raise FileReaderError("--Background_Profile parameter not found. Is this really a dark measurement file?")

        sig_arr = raw_data["Background_Profile"].values
    else:
        if "Raw_Lidar_Data" not in list(raw_data.variables):
            raise FileReaderError("--Raw_Lidar_Data parameter not found. Is this really a non-dark measurement file?")

        sig_arr = raw_data["Raw_Lidar_Data"].values

    sig_arr[sig_arr >= 9.96e+36] = np.nan

    bins = 1. + np.arange(0, sig_arr.shape[-1])

    sig_raw = xr.DataArray(sig_arr,
                           coords=[time, channels, bins], #range_sig
                           dims=['time', 'channel', 'bins']) #'range'

    # Sort by time
    sig_raw = sig_raw.copy().sortby('time')

    return(sig_raw)


def read_shots(raw_data, time, channels, meas_type):

    if meas_type == "drk":
        shots = np.tile(np.median(raw_data["Laser_Shots"].values, axis = 0), (len(time), 1))
    else:
        shots = raw_data["Laser_Shots"].values

    shots = xr.DataArray(shots,
                         coords=[time, channels], #range_sig
                         dims=['time', 'channel']) #'range'

    # Sort by time
    shots = shots.copy().sortby('time')

    return(shots)


# -----------------------------------------------------------------------------
# Tamarin NetCDF group support through xarray only
# -----------------------------------------------------------------------------

def _decode_nc_value(value):
    """Return a Python scalar/string for NetCDF scalar attributes/values."""

    if isinstance(value, bytes):
        return value.decode().strip()

    if isinstance(value, np.bytes_):
        return bytes(value).decode().strip()

    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _decode_nc_value(value.item())
        if value.size == 1:
            return _decode_nc_value(value.reshape(-1)[0])

    return value


def _open_tamarin_dataset(nc_file, group = None):
    """Open a Tamarin NetCDF group with xarray.

    Note: hierarchical NetCDF4 groups require an xarray backend such as
    netcdf4 or h5netcdf. This avoids importing h5py directly.
    """

    kwargs = {}
    if group is not None:
        kwargs["group"] = group

    try:
        return xr.open_dataset(nc_file, **kwargs)
    except ValueError as exc:
        raise FileReaderError(
            "--Could not open this NetCDF4 group with xarray. "
            "Install an xarray NetCDF4-capable backend, e.g. netCDF4 or h5netcdf, "
            "or use the h5py-based reader."
        ) from exc


def _try_open_tamarin_group(nc_file, group_path):
    try:
        return _open_tamarin_dataset(nc_file, group = group_path)
    except Exception:
        return None


def is_tamarin_file(nc_file):
    """Detect the Tamarin hierarchical NetCDF flavor using xarray only."""

    try:
        root = _open_tamarin_dataset(nc_file)
    except Exception:
        return False

    system = str(_decode_nc_value(root.attrs.get("System", ""))).upper()
    if system == "TAMARIN":
        return True

    # Conservative fallback: a baseline calibration group without SCC attributes.
    if "Measurement_ID" not in root.attrs:
        baseline = _try_open_tamarin_group(nc_file, "calibrations/baseline")
        if baseline is not None:
            return True

    return False


def _group_has_signal_variables(ds):
    if ds is None:
        return False
    if "analog_mean" in ds.variables:
        return True
    return any(str(name).startswith("count_") for name in ds.variables)


def _tamarin_candidate_groups(meas_type):
    if meas_type == "drk":
        return ["calibrations/baseline"]

    # Include both singular/plural variants and their lidar_signal subgroup.
    return [
        "atmospheric_signals/lidar_signal",
        "atmospheric_signal/lidar_signal",
        "atmospheric_signals",
        "atmospheric_signal",
    ]


def _tamarin_group_path(nc_file, meas_type):
    """Return the Tamarin group containing the requested profiles."""

    for group_path in _tamarin_candidate_groups(meas_type):
        ds = _try_open_tamarin_group(nc_file, group_path)
        if _group_has_signal_variables(ds):
            return group_path

    return None


def _tamarin_count_variable(ds):
    """Find the photon-counting profile variable, if present."""

    preferred = ["count_vec", "count_mean", "count_sum", "count_max", "count_min", "count_std"]
    for name in preferred:
        if name in ds.variables:
            return name

    count_names = [name for name in ds.variables if str(name).startswith("count_")]
    if len(count_names) == 0:
        return None

    profile_count_names = []
    for name in count_names:
        if ds[name].ndim == 3:
            profile_count_names.append(name)

    if len(profile_count_names) > 0:
        return sorted(profile_count_names)[0]

    return sorted(count_names)[0]


def _read_string_array(values):
    arr = np.asarray(values)
    out = []
    for val in arr.reshape(-1):
        val = _decode_nc_value(val)
        out.append(str(val))
    return np.asarray(out, dtype = object)


def _tamarin_base_channel_names(root):
    if "channel_name" in root.variables:
        names = _read_string_array(root["channel_name"].values)
    elif "channels" in root.variables:
        names = _read_string_array(root["channels"].values)
    elif "channel_ID" in root.variables:
        names = _read_string_array(root["channel_ID"].values)
    elif "DAQ_Range" in root.variables:
        names = np.asarray([str(i + 1) for i in range(root["DAQ_Range"].shape[0])], dtype = object)
    elif "channel_voltage_for_detector_gain" in root.variables:
        names = np.asarray([str(i + 1) for i in range(root["channel_voltage_for_detector_gain"].shape[0])], dtype = object)
    else:
        raise FileReaderError("--Could not identify Tamarin channels")

    return names.astype(str)


def _tamarin_count_variable_name(nc_file, meas_type):
    """Return the Tamarin photon-counting profile variable name, if present."""

    group_path = _tamarin_group_path(nc_file, meas_type)
    if group_path is None:
        return None

    ds = _open_tamarin_dataset(nc_file, group = group_path)
    return _tamarin_count_variable(ds)


def _tamarin_signal_variables(nc_file, meas_type):
    """Return the fixed Tamarin output signal layout.

    The output channel layout must be identical for normal and dark profiles:
    first all analog channels, then the corresponding photon channels. If a
    Tamarin group does not contain photon profiles, read_signals_tamarin fills
    those photon profiles with zeros.
    """

    group_path = _tamarin_group_path(nc_file, meas_type)
    if group_path is None:
        return []

    ds = _open_tamarin_dataset(nc_file, group = group_path)
    count_name = _tamarin_count_variable(ds)

    variables = []
    variables.append(("ana", "analog_mean" if "analog_mean" in ds.variables else None))
    variables.append(("cpt", count_name))

    return variables


def _make_unique_channel_names(names):
    """Return unique string channel names while preserving order."""

    names = [str(name) for name in names]
    seen = {}
    unique = []

    for name in names:
        if name not in seen:
            seen[name] = 0
            unique.append(name)
        else:
            seen[name] += 1
            unique.append(f"{name}_{seen[name] + 1}")

    return np.asarray(unique, dtype = object)


def _tamarin_channel_names(nc_file, meas_type):
    """Return stable SCC-style ascending channel indices for Tamarin output.

    The exported layout is fixed for every measurement type:

        1..N       analog channels
        N+1..2N    photon channels

    This keeps normal/ray and dark measurements on the same channel coordinate,
    even when dark files contain only analog profiles.
    """

    root = _open_tamarin_dataset(nc_file)
    n_base_channels = _tamarin_base_channel_names(root).size
    n_output_channels = n_base_channels * _tamarin_output_multiplier(nc_file, meas_type)

    return np.asarray([str(i + 1) for i in range(n_output_channels)], dtype = object)


def _tamarin_output_multiplier(nc_file, meas_type):
    """Number of output blocks per physical channel for Tamarin.

    Tamarin output is always analog + photon so that normal and dark profiles
    expose the same channel set.
    """

    return 2


def _tamarin_daq_ranges(nc_file, meas_type):
    root = _open_tamarin_dataset(nc_file)

    if "DAQ_Range" in root.variables:
        daq = np.asarray(root["DAQ_Range"].values)
    elif "channel_acquisition_input_range" in root.variables:
        daq = np.asarray(root["channel_acquisition_input_range"].values)
    elif "channel_voltage_for_detector_gain" in root.variables:
        daq = np.asarray(root["channel_voltage_for_detector_gain"].values)
    else:
        daq = np.full(_tamarin_base_channel_names(root).size, np.nan)

    return np.tile(daq, _tamarin_output_multiplier(nc_file, meas_type))


def _repeat_tamarin_channel_values(values, nc_file, meas_type):
    """Repeat per-physical-channel metadata for analog and photon channels."""

    values = np.asarray(values)
    return np.tile(values, _tamarin_output_multiplier(nc_file, meas_type))


def _tamarin_group_scalar(nc_file, meas_type, name, default = np.nan):
    """Read a scalar value from the Tamarin signal/baseline group."""

    group_path = _tamarin_group_path(nc_file, meas_type)
    if group_path is None:
        return default

    group = _open_tamarin_dataset(nc_file, group = group_path)

    if name in group.variables:
        return _decode_nc_value(group[name].values)

    if name in group.attrs:
        return _decode_nc_value(group.attrs[name])

    return default


def _tamarin_root_scalar(root, name, default = np.nan):
    """Read a scalar value from the Tamarin root group variables or attributes."""

    if name in root.variables:
        return _decode_nc_value(root[name].values)

    if name in root.attrs:
        return _decode_nc_value(root.attrs[name])

    return default


def _tamarin_first_scalar(value, default = np.nan):
    """Reduce a scalar or array-like NetCDF value to one Python scalar.

    Some Tamarin metadata, such as Laser_Repetition_Rate, is stored per
    channel even though the SCC system_info field expects one scalar. In that
    case, use the first non-NaN value.
    """

    value = _decode_nc_value(value)

    if isinstance(value, np.ndarray):
        arr = value.reshape(-1)
    else:
        arr = np.asarray(value).reshape(-1)

    if arr.size == 0:
        return default

    for item in arr:
        item = _decode_nc_value(item)
        try:
            if pd.isna(item):
                continue
        except Exception:
            pass
        return item

    return default


def _tamarin_laser_repetition_rate(nc_file, meas_type):
    """Read Laser_Repetition_Rate and reduce it to one scalar value."""

    root = _open_tamarin_dataset(nc_file)

    if "Laser_Repetition_Rate" in root.variables:
        return _tamarin_first_scalar(root["Laser_Repetition_Rate"].values)

    if "Laser_Repetition_Rate" in root.attrs:
        return _tamarin_first_scalar(root.attrs["Laser_Repetition_Rate"])

    return _tamarin_first_scalar(
        _tamarin_group_scalar(nc_file, meas_type, "Laser_Repetition_Rate", default = np.nan)
    )


def _tamarin_acquisition_modes(nc_file, meas_type):
    """Return SCC-style channel type labels: 'a' for analog, 'p' for photon."""

    root = _open_tamarin_dataset(nc_file)
    n_channels = _tamarin_base_channel_names(root).size
    return np.asarray((n_channels * ["a"]) + (n_channels * ["p"]), dtype = object)


def read_meas_tamarin(nc_file, meas_type = None):
    root = _open_tamarin_dataset(nc_file)

    system_info = pd.Series(dtype = object)

    system_info.loc["laser_A_repetition_rate"] = _tamarin_laser_repetition_rate(nc_file, meas_type)
    system_info.loc["laser_B_repetition_rate"] = np.nan
    system_info.loc["laser_C_repetition_rate"] = np.nan

    system_info.loc["lidar_name"] = _tamarin_root_scalar(root, "System")
    system_info.loc["station_name"] = _tamarin_root_scalar(root, "Location")
    system_info.loc["station_altitude"] = _tamarin_root_scalar(root, "Altitude_meter_asl")
    system_info.loc["station_latitude"] = _tamarin_root_scalar(root, "Latitude_degrees_north")
    system_info.loc["station_longitude"] = _tamarin_root_scalar(root, "Longitude_degrees_east")
    system_info.loc["zenith_angle"] = _tamarin_root_scalar(root, "Laser_Pointing_Angle")

    return(system_info)


def read_channels_tamarin(nc_file, meas_type):
    root = _open_tamarin_dataset(nc_file)
    ch_index = _tamarin_channel_names(nc_file, meas_type).astype(str)
    n_output_channels = ch_index.size

    channel_info = pd.DataFrame(index=ch_index)

    # existing fields...
    channel_info["data_acquisition_range"] = _tamarin_daq_ranges(nc_file, meas_type)

    # add this
    laser_repetition_rate = _tamarin_laser_repetition_rate(nc_file, meas_type)
    channel_info.loc[:, "laser_repetition_rate"] = n_output_channels * [laser_repetition_rate]

    if "Detected_Wavelength" in root.variables:
        detected_wavelength = root["Detected_Wavelength"].values
        detected_wavelength = _repeat_tamarin_channel_values(detected_wavelength, nc_file, meas_type)
    elif "channel_wavelength" in root.variables:
        detected_wavelength = root["channel_wavelength"].values
        detected_wavelength = _repeat_tamarin_channel_values(detected_wavelength, nc_file, meas_type)
    else:
        detected_wavelength = np.full(n_output_channels, np.nan)

    if "Emitted_Wavelength" in root.variables:
        emitted_wavelength = root["Emitted_Wavelength"].values
        emitted_wavelength = _repeat_tamarin_channel_values(emitted_wavelength, nc_file, meas_type)
    else:
        emitted_wavelength = np.full(n_output_channels, np.nan)

    channel_info.loc[:, "detected_wavelength"] = detected_wavelength
    channel_info.loc[:, "emitted_wavelength"] = emitted_wavelength
    channel_info.loc[:, "acquisition_mode"] = _tamarin_acquisition_modes(nc_file, meas_type)
    channel_info.loc[:, "channel_type"] = _tamarin_acquisition_modes(nc_file, meas_type)

    range_resolution = _tamarin_group_scalar(nc_file, meas_type, "_distance_resolution", default = np.nan)
    channel_info.loc[:, "range_resolution"] = n_output_channels * [range_resolution]

    channel_info.loc[:, "laser"] = n_output_channels * ["1"]

    return(channel_info)


def _normalise_dim_name(dim):
    """Return a lower-case NetCDF dimension name for robust comparisons."""

    return str(dim).strip().lower()


def _find_dim(dims, accepted, role):
    """Find exactly one dimension by name, not by length."""

    matches = [dim for dim in dims if _normalise_dim_name(dim) in accepted]

    if len(matches) != 1:
        raise FileReaderError(
            f"--Could not identify exactly one {role} dimension from {dims}. "
            f"Expected one of {sorted(accepted)}"
        )

    return matches[0]


def _profile_to_time_channel_bins(profile_var, n_channels, n_times):
    """Convert a Tamarin profile DataArray to a NumPy array with axes
    (time, channel, bins), using NetCDF dimension names.

    Tamarin profile variables are expected to use dimensions equivalent to
    (time, channels, points). The output reader convention is
    (time, channel, bins), so the NetCDF dimension named "points" is mapped to
    output "bins".
    """

    if profile_var.ndim != 3:
        raise FileReaderError(
            f"--Expected a 3D Tamarin signal variable, got dims "
            f"{profile_var.dims} and shape {profile_var.shape}"
        )

    dims = profile_var.dims
    time_dim = _find_dim(dims, {"time", "time_bck"}, "time")
    channel_dim = _find_dim(dims, {"channels", "channel", "channel_id"}, "channel")
    bin_dim = _find_dim(dims, {"points", "point", "bins", "bin", "range"}, "bin")

    arr = profile_var.transpose(time_dim, channel_dim, bin_dim).values.astype(float)

    if arr.shape[0] != n_times:
        raise FileReaderError(
            f"--Tamarin signal time dimension has length {arr.shape[0]}, "
            f"expected {n_times}. Variable dims: {profile_var.dims}"
        )

    if arr.shape[1] != n_channels:
        raise FileReaderError(
            f"--Tamarin signal channel dimension has length {arr.shape[1]}, "
            f"expected {n_channels}. Variable dims: {profile_var.dims}"
        )

    return arr


def _as_time_channel(arr, n_channels, n_times):
    """Convert a per-profile/per-channel array to (time, channel).

    This helper is kept for generic arrays. For Tamarin Laser_Shots,
    read_shots_tamarin uses _laser_shots_to_time_channel so the values are
    read directly from the group variable, using the variable dimensions when
    available.
    """

    arr = np.asarray(arr)
    if arr.ndim == 1:
        if arr.size == n_channels:
            return np.tile(arr, (n_times, 1))
        if arr.size == n_times:
            return np.tile(arr.reshape(-1, 1), (1, n_channels))

    if arr.ndim != 2:
        raise FileReaderError(f"--Expected a 1D or 2D Tamarin shots array, got shape {arr.shape}")

    if arr.shape == (n_times, n_channels):
        return arr
    if arr.shape == (n_channels, n_times):
        return arr.T

    raise FileReaderError(f"--Could not convert Tamarin shots shape {arr.shape} to (time, channel)")


def _laser_shots_to_time_channel(shots_var, n_channels, n_times):
    """Read Tamarin Laser_Shots directly as a (time, channel) array.

    Tamarin stores Laser_Shots in the same group as the signal, typically with
    dimensions such as (time, points). This function does not take medians,
    sums, or totals. It only reorders axes when needed to match the SCC reader
    output format.
    """

    arr = np.asarray(shots_var.values).astype(float)
    dims = tuple(str(dim).lower() for dim in shots_var.dims)

    if arr.ndim != 2:
        raise FileReaderError(
            f"--Expected Tamarin Laser_Shots to be 2D with dimensions like "
            f"(time, points), got shape {arr.shape}"
        )

    time_axis = None
    channel_axis = None

    for ax, dim in enumerate(dims):
        if dim in ["time", "time_bck"] or "time" in dim:
            time_axis = ax
        if dim in ["channel", "channels", "channel_id"] or "channel" in dim:
            channel_axis = ax

    if time_axis is None:
        matches = [ax for ax, size in enumerate(arr.shape) if size == n_times]
        if len(matches) == 1:
            time_axis = matches[0]

    if channel_axis is None:
        matches = [ax for ax, size in enumerate(arr.shape) if size == n_channels and ax != time_axis]
        if len(matches) == 1:
            channel_axis = matches[0]

    if time_axis is None or channel_axis is None or time_axis == channel_axis:
        raise FileReaderError(
            f"--Could not identify time/channel axes for Tamarin Laser_Shots "
            f"with dims {shots_var.dims} and shape {arr.shape}"
        )

    arr = np.transpose(arr, (time_axis, channel_axis))

    if arr.shape != (n_times, n_channels):
        raise FileReaderError(
            f"--Tamarin Laser_Shots has shape {arr.shape} after transpose, "
            f"expected ({n_times}, {n_channels})"
        )

    return arr


def _get_root_time_attr(root, name):
    if name in root.attrs:
        return _decode_nc_value(root.attrs[name])
    if name in root.variables:
        return _decode_nc_value(root[name].values)
    return None


def get_time_info_tamarin(nc_file, meas_type, filename):
    time_info = pd.Series(dtype = object)
    group_path = _tamarin_group_path(nc_file, meas_type)

    if group_path is None:
        if meas_type == "drk":
            print("--calibrations/baseline group not found. No Tamarin dark profile -> skipping")
        else:
            print("--atmospheric_signals group not found. No Tamarin atmospheric profile -> skipping")
        return(time_info)

    group = _open_tamarin_dataset(nc_file, group = group_path)
    root = _open_tamarin_dataset(nc_file)

    if "Raw_Data_Start_Time" not in group.variables:
        print(f"--Raw_Data_Start_Time not found in Tamarin group {group_path} -> skipping")
        return(time_info)

    if "Raw_Data_Stop_Time" not in group.variables:
        print(f"--Raw_Data_Stop_Time not found in Tamarin group {group_path} -> skipping")
        return(time_info)

    sdate = _get_root_time_attr(root, "RawData_Start_Date")
    stime = _get_root_time_attr(root, "RawData_Start_Time_UT")

    if sdate is None or stime is None:
        raise FileReaderError("--RawData_Start_Date or RawData_Start_Time_UT not found in Tamarin file attributes")

    start_time_sec = np.asarray(group["Raw_Data_Start_Time"].values)
    stop_time_sec = np.asarray(group["Raw_Data_Stop_Time"].values)

    start_time_sec = start_time_sec.reshape(start_time_sec.shape[0], -1)[:,0].astype(float)
    stop_time_sec = stop_time_sec.reshape(stop_time_sec.shape[0], -1)[:,0].astype(float)

    sdt = dt.datetime.strptime(str(sdate) + ' ' + str(stime), "%Y%m%d %H%M%S")

    start_time_arr = np.array([sdt + dt.timedelta(seconds = t) for t in start_time_sec])
    end_time_arr = np.array([sdt + dt.timedelta(seconds = t) for t in stop_time_sec])

    filenames = np.empty(start_time_arr.size, dtype = object)
    filenames[:] = filename

    tdata = np.array([filenames,
                      datetimes_to_iso(start_time_arr),
                      datetimes_to_iso(end_time_arr)],
                     dtype = object)

    time_info = pd.DataFrame(tdata.T,
                             index = start_time_arr,
                             columns = ['filename', 'start_time', 'end_time'])

    return(time_info)


def read_signals_tamarin(nc_file, time, channels, meas_type):
    group_path = _tamarin_group_path(nc_file, meas_type)
    variables = _tamarin_signal_variables(nc_file, meas_type)

    if group_path is None or len(variables) == 0:
        raise FileReaderError("--No Tamarin profile variables found")

    group = _open_tamarin_dataset(nc_file, group = group_path)
    root = _open_tamarin_dataset(nc_file)
    n_channels = _tamarin_base_channel_names(root).size
    n_times = len(time)

    template = None
    sig_blocks = []

    for mode, var_name in variables:
        if var_name is not None and var_name in group.variables:
            sig_arr = _profile_to_time_channel_bins(group[var_name],
                                                    n_channels = n_channels,
                                                    n_times = n_times)
            template = sig_arr
        else:
            if template is None:
                # Find any existing 3D profile to infer the number of bins.
                for _, candidate in variables:
                    if candidate is not None and candidate in group.variables:
                        template = _profile_to_time_channel_bins(
                            group[candidate],
                            n_channels = n_channels,
                            n_times = n_times)
                        break

            if template is None:
                raise FileReaderError("--Could not infer Tamarin profile shape for synthetic channels")

            # Missing photon channels, for example in dark/baseline profiles,
            # are exported as synthetic zero profiles so all meas_types have the
            # same channel set.
            sig_arr = np.zeros_like(template, dtype = float)

        sig_blocks.append(sig_arr)

    sig_arr = np.concatenate(sig_blocks, axis = 1)
    sig_arr[sig_arr >= 9.96e+36] = np.nan

    bins = 1. + np.arange(0, sig_arr.shape[-1])

    sig_raw = xr.DataArray(sig_arr,
                           coords=[time, channels, bins],
                           dims=['time', 'channel', 'bins'])

    sig_raw = sig_raw.copy().sortby('time')

    return(sig_raw)


def read_shots_tamarin(nc_file, time, channels, meas_type):
    group_path = _tamarin_group_path(nc_file, meas_type)

    if group_path is None:
        raise FileReaderError("--No Tamarin group found for requested measurement type")

    group = _open_tamarin_dataset(nc_file, group = group_path)
    root = _open_tamarin_dataset(nc_file)
    n_channels = _tamarin_base_channel_names(root).size
    n_times = len(time)

    if "Laser_Shots" not in group.variables:
        shots_arr = np.full((n_times, n_channels), np.nan)
    else:
        shots_arr = _laser_shots_to_time_channel(group["Laser_Shots"],
                                                 n_channels = n_channels,
                                                 n_times = n_times)

    # Tamarin output is always analog + photon. Duplicate the directly-read
    # per-profile shots so the shots DataArray keeps the exact same
    # (time, channel) shape as sig_raw, including synthetic photon channels.
    shots_arr = np.tile(shots_arr, (1, _tamarin_output_multiplier(nc_file, meas_type)))

    shots = xr.DataArray(shots_arr,
                         coords=[time, channels],
                         dims=['time', 'channel'])

    shots = shots.copy().sortby('time')

    return(shots)
