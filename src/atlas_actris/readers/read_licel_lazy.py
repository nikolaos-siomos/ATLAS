import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime as dt
from datetime import timedelta
import xarray as xr
import dask.array as da
from dask import delayed
from readers.check_file_format import is_licel_header
from utils.error_classes import FileReaderError
from utils.time_conversions import datetimes_to_iso
from utils.error_classes import CustomWarning


HEADER_BODY_SEPARATOR = b"\r\n\r\n"
SIGNAL_BLOCK_SEPARATOR = b"\r\n"
SIGNAL_DTYPE = "<u4"  # legacy int.from_bytes(..., byteorder="little") was unsigned


def find_sep_bytes(raw, fname):
    """Return the byte offset where the Licel header/body separator starts."""

    sep = raw.find(HEADER_BODY_SEPARATOR)
    if sep < 0:
        raise FileReaderError(
            f"Could not find header/data separator. Is {fname} a Licel file?"
        )
    return sep


def find_sep_file(buffer, fname):
    """Identify the header/body separator field: CRLF CRLF."""

    old_pos = buffer.tell()
    buffer.seek(0)
    raw = buffer.read()
    buffer.seek(old_pos)
    return find_sep_bytes(raw, fname)


def read_body_from_bytes(raw, channel_info, output_channels, fname, dtype=SIGNAL_DTYPE):
    """Fast binary body reader using explicit header/body and channel separators.

    The Licel body is parsed for *all* channels listed in ``channel_info`` so
    that auxiliary channels are accounted for and separators can be validated.
    Only ``output_channels`` are copied to the returned array. This keeps the
    public output restricted to BT/BC channels while preserving correct byte
    alignment over the full raw-file body.
    """

    sep = find_sep_bytes(raw, fname)
    body = raw[sep + len(HEADER_BODY_SEPARATOR):]

    bins_all = channel_info.loc[:, "bins"].values.astype(int)
    channels_all = np.asarray(channel_info.index.values)
    output_channels = list(output_channels)

    missing = [ch for ch in output_channels if ch not in channel_info.index]
    if missing:
        raise FileReaderError(
            f"--Requested output channels are missing from {fname}: {missing}"
        )

    output_mask = np.isin(channels_all, output_channels)
    output_bins = channel_info.loc[output_channels, "bins"].values.astype(int)

    sig_arr = np.nan * np.zeros(
        (len(output_channels), int(np.max(output_bins))),
        dtype=float,
    )

    pos = 0
    out_i = 0

    for ch_i, n_bins in enumerate(bins_all):
        n_bytes = 4 * int(n_bins)
        end = pos + n_bytes

        if end > len(body):
            raise FileReaderError(
                f"--Could not read {n_bins} bins for channel {channels_all[ch_i]} "
                f"in {fname}. The binary body ended earlier than expected."
            )

        if output_mask[ch_i]:
            vals = np.frombuffer(
                body,
                dtype=dtype,
                count=int(n_bins),
                offset=pos,
            )
            sig_arr[out_i, :int(n_bins)] = vals.astype(float)
            out_i += 1

        pos = end

        if body[pos:pos + len(SIGNAL_BLOCK_SEPARATOR)] != SIGNAL_BLOCK_SEPARATOR:
            raise FileReaderError(
                f"--Unexpected signal separator after channel {channels_all[ch_i]} "
                f"in {fname} at body byte position {pos}."
            )

        pos += len(SIGNAL_BLOCK_SEPARATOR)

    if pos != len(body):
        trailing = body[pos:]
        if trailing.strip(b"\x00\r\n\t "):
            raise FileReaderError(
                f"--Unexpected body size in {fname}. Parsed {pos} bytes, "
                f"but body contains {len(body)} bytes."
            )

    return sig_arr


def read_body(buffer, channel_info, output_channels, fname, dtype=SIGNAL_DTYPE):
    """Read selected signal channels from an open binary file object."""

    buffer.seek(0)
    raw = buffer.read()
    return read_body_from_bytes(
        raw=raw,
        channel_info=channel_info,
        output_channels=output_channels,
        fname=fname,
        dtype=dtype,
    )


def read_header(buffer):
    """Retrieves system and timing information from the Licel header."""

    system_info = pd.Series()

    # Just skip the first line
    buffer.readline()

    # Get the information from the 2nd line
    secondline = str(buffer.readline(), encoding="utf-8").split()
    file_format = "new"

    # Check if the file has the old 2 rack format
    if len(secondline) == 1:
        buffer.readline()
        secondline = str(buffer.readline(), encoding="utf-8").split()
        file_format = "old"

    start_date = secondline[1]
    start_time = secondline[2]

    end_date = secondline[3]
    end_time = secondline[4]

    stime = dt.strptime(start_date + " " + start_time, "%d/%m/%Y %H:%M:%S")
    etime = dt.strptime(end_date + " " + end_time, "%d/%m/%Y %H:%M:%S")

    system_info["station_altitude"] = float(secondline[5])
    system_info["station_latitude"] = np.round(float(secondline[6]), 4)
    system_info["station_longitude"] = np.round(float(secondline[7]), 4)

    if len(secondline) > 8:
        system_info["zenith_angle"] = float(secondline[8])

    if len(secondline) > 9:
        system_info["azimuth_angle"] = float(secondline[9])

    # Get the information from the 3rd line
    thirdline = str(buffer.readline(), encoding="utf-8").split()

    num_channels = int(thirdline[4])

    system_info["laser_A_repetition_rate"] = float(thirdline[1])

    if len(thirdline) > 2:
        system_info["laser_B_repetition_rate"] = float(thirdline[3])
    else:
        system_info["laser_B_repetition_rate"] = np.nan

    if len(thirdline) > 5:
        system_info["laser_C_repetition_rate"] = float(thirdline[6])
    else:
        system_info["laser_C_repetition_rate"] = np.nan

    return system_info, stime, etime, num_channels, file_format


def read_channels(buffer, num_channels, file_format):
    """Collect channel-specific information from the Licel header.

    All channel metadata rows are parsed here. Filtering to BT/BC happens only
    after the full binary body layout is known, so auxiliary channels do not
    disturb byte alignment. Auxiliary channels may differ between files.
    """

    if file_format == "new":
        cols = [
            "active",
            "acquisition_mode",
            "laser_id",
            "bins",
            "laser_polarization",
            "pmt_high_voltage",
            "range_resolution",
            "wave_pol",
            "unk1",
            "unk2",
            "unk3",
            "unk4",
            "analog_to_digital_resolution",
            "shots",
            "data_acquisition_range",
            "recorder_channel_id",
        ]
    else:
        cols = [
            "active",
            "acquisition_mode",
            "laser_id",
            "bins",
            "laser_polarization",
            "pmt_high_voltage",
            "range_resolution",
            "wave_pol",
            "analog_to_digital_resolution",
            "shots",
            "data_acquisition_range",
            "recorder_channel_id",
            "unk1",
            "unk2",
            "unk3",
        ]

    temp_info = []

    for _ in range(num_channels):
        linevars = str(buffer.readline(), encoding="utf-8").split()
        if len(linevars) < len(cols):
            raise FileReaderError(
                "--Malformed Licel channel header line. "
                f"Expected at least {len(cols)} columns but got {len(linevars)}."
            )
        temp_info.append(linevars[:len(cols)])

    temp_info = pd.DataFrame(
        temp_info,
        index=np.arange(len(temp_info)),
        columns=cols,
    )

    header_channel_id = temp_info.loc[:, "recorder_channel_id"].values
    laser_id = temp_info.loc[:, "laser_id"].values

    # Combine from the recorder channel ID and laser number only if needed.
    if len(header_channel_id) == len(set(header_channel_id)):
        recorder_channel_id = header_channel_id
    else:
        recorder_channel_id = [
            f"{ch_id}_L{lr_id}" for ch_id, lr_id in zip(header_channel_id, laser_id)
        ]

    if len(recorder_channel_id) != len(set(recorder_channel_id)):
        raise FileReaderError(
            "-- Error: At least two of the Licel channels have both the same id "
            "and laser number. Please correct this in the recorder settings"
        )

    channel_info = pd.DataFrame(index=recorder_channel_id)

    info_columns = [
        "acquisition_mode",
        "laser_id",
        "bins",
        "range_resolution",
        "shots",
        "data_acquisition_range",
        "analog_to_digital_resolution",
        "recorder_channel_id",
    ]

    channel_info.loc[:, info_columns] = temp_info.loc[:, info_columns].copy().values.astype(object)

    mask_an = channel_info.loc[:, "acquisition_mode"].values == "0"

    channel_info.loc[mask_an, "data_acquisition_range"] = (
        1000.0 * channel_info.loc[mask_an, "data_acquisition_range"].astype(float)
    )
    channel_info.loc[~mask_an, "data_acquisition_range"] = None

    wave = np.array(
        list(np.char.split(temp_info.wave_pol.values.astype("str"), sep="."))
    )[:, 0].astype(float)

    channel_info.loc[:, "detected_wavelength"] = wave
    channel_info.loc[:, "dead_time_correction_type"] = 0.0

    channel_info["acquisition_mode"] = channel_info["acquisition_mode"].astype("object")
    channel_info.loc[mask_an, "acquisition_mode"] = "a"
    channel_info.loc[~mask_an, "acquisition_mode"] = "p"

    return channel_info


def output_channel_names(channel_info):
    """Return BT/BC channels to be exposed by the reader."""

    channels = np.asarray(channel_info.index.values)
    keep = np.array([str(ch).startswith(("BT", "BC")) for ch in channels])

    if not keep.any():
        raise FileReaderError("--No BT/BC signal channels found in the Licel header")

    return list(channels[keep])


def read_buffer(fname):
    """Reads the binary file as a single byte sequence."""

    with open(fname, "rb") as f:
        buffer = f.read()

    return buffer


def unit_conv_bits_to_mV(channel_info, signal, shots):
    """Convert analog signals from bits to mV without eager signal loading."""

    if len(signal) == 0:
        return signal

    factors = xr.ones_like(shots, dtype=float)

    mask_an = channel_info.acquisition_mode.values == "a"
    channel_id_an = channel_info.index.values[mask_an]

    data_acquisition_range = channel_info.data_acquisition_range.astype(float)
    analog_to_digital_resolution = channel_info.analog_to_digital_resolution.astype(float)

    for ch in channel_id_an:
        ch_d = dict(channel=ch)
        factors.loc[ch_d] = (
            data_acquisition_range.loc[ch]
            / (
                shots.loc[ch_d]
                * (np.power(2, analog_to_digital_resolution.loc[ch]) - 1.0)
            )
        )

    return signal * factors


def read_single_file_signal(
    fname,
    expected_output_channels,
    expected_output_bins,
    expected_file_format,
):
    """Read one Licel file signal lazily.

    All channel blocks present in this file are parsed and separator-checked.
    Only the expected BT/BC channels are returned in the array. Auxiliary
    channels are allowed to differ between files.
    """

    with open(fname, "rb") as buffer:
        _, _, _, num_channels_i, file_format_i = read_header(buffer)

        if file_format_i != expected_file_format:
            raise FileReaderError(f"--Not all files have the same Licel format: {fname}")

        channel_info_i = read_channels(
            buffer=buffer,
            num_channels=num_channels_i,
            file_format=file_format_i,
        )

        output_channels_i = output_channel_names(channel_info_i)

        if list(output_channels_i) != list(expected_output_channels):
            raise FileReaderError(
                f"--Not all files have the same BT/BC output channels: {fname}"
            )

        bins_i = channel_info_i.loc[expected_output_channels, "bins"].values.astype(int)
        if not np.array_equal(bins_i, expected_output_bins):
            raise FileReaderError(
                f"--Not all files include BT/BC channels with the same number of bins: {fname}"
            )

        buffer.seek(0)
        raw = buffer.read()

    return read_body_from_bytes(
        raw=raw,
        channel_info=channel_info_i,
        output_channels=expected_output_channels,
        fname=fname,
    )


# Read measurement
def read_dataset(dir_meas, meas_type=None, lazy=True, chunks=None):
    """Read raw Licel files.

    Returned ``channel_info``, ``shots``, and ``sig_raw`` include only BT/BC
    channels. The binary body is still parsed over all raw-file channels to keep
    the byte alignment correct and to validate the per-channel CRLF separators.

    Parameters
    ----------
    dir_meas : str
        Folder containing raw Licel files.
    meas_type : optional
        Kept for compatibility with existing calls.
    lazy : bool, optional
        If True, ``sig_raw`` is Dask-backed and each raw file is read only when
        needed. If False, signal bodies are loaded eagerly.
    chunks : dict or None, optional
        Reader-level xarray chunks. By default: time chunks between 10 and 50,
        all channels together, and 4096 bins.
    """

    sig_raw = []
    shots = []
    start_time_arr = []
    end_time_arr = []
    filename = []

    system_info = []
    channel_info = []
    time_info = []

    if not os.path.exists(dir_meas):
        print(
            "---- Warning : The folder for reading signals does not exist! "
            f"Check the input directory! \n Given folder: {dir_meas}"
        )
        return system_info, channel_info, time_info, sig_raw, shots

    mfiles = glob.glob(os.path.join(dir_meas, "*.*"))
    mfiles = [file for file in mfiles if os.path.basename(file) != "temp.dat"]

    if len(mfiles) == 0:
        CustomWarning(f"No files to read in: {dir_meas}")
        return system_info, channel_info, time_info, sig_raw, shots

    print(f"-- Reading {len(mfiles)} file(s)!")

    if not is_licel_header(mfiles[0]):
        raise FileReaderError(f"--QA test folder contains non Licel files: {dir_meas}")

    with open(mfiles[0], "rb") as buffer:
        system_info, stime, etime, num_channels, file_format = read_header(buffer)
        all_channel_info = read_channels(
            buffer=buffer,
            num_channels=num_channels,
            file_format=file_format,
        )

    all_channels = all_channel_info.index
    channels = output_channel_names(all_channel_info)
    channel_info = all_channel_info.loc[channels, :].copy()

    bins = channel_info.loc[:, "bins"].values.astype(int)

    # Add repetition rate to the returned BT/BC channel_info only.
    for ch in channels:
        if channel_info.loc[ch, "laser_id"] == "1" and system_info["laser_A_repetition_rate"] is not None:
            channel_info.loc[ch, "laser_repetition_rate"] = system_info["laser_A_repetition_rate"]
        if channel_info.loc[ch, "laser_id"] == "2" and system_info["laser_B_repetition_rate"] is not None:
            channel_info.loc[ch, "laser_repetition_rate"] = system_info["laser_B_repetition_rate"]
        if channel_info.loc[ch, "laser_id"] == "3" and system_info["laser_C_repetition_rate"] is not None:
            channel_info.loc[ch, "laser_repetition_rate"] = system_info["laser_C_repetition_rate"]

    bins_arr = np.arange(1, max(bins) + 1, 1)

    start_time_arr = np.nan * np.zeros(len(mfiles), dtype=object)
    end_time_arr = np.nan * np.zeros(len(mfiles), dtype=object)
    shots_arr = np.nan * np.zeros((len(mfiles), len(channels)), dtype=object)
    filename = np.empty(len(mfiles), dtype=object)

    lazy_signal_blocks = []
    eager_signal_blocks = []

    for k, fname in enumerate(mfiles):
        filename[k] = os.path.basename(fname)

        if not is_licel_header(fname):
            raise FileReaderError(f"--The following file is not in raw Licel format: {fname}")

        with open(fname, "rb") as buffer:
            _, stime_i, etime_i, num_channels_i, file_format_i = read_header(buffer)

            if file_format_i != file_format:
                raise FileReaderError(
                    f"--Not all files have the same Licel format: {fname}\nCompare with: {mfiles[0]}"
                )

            channel_info_i = read_channels(
                buffer=buffer,
                num_channels=num_channels_i,
                file_format=file_format_i,
            )

            output_channels_i = output_channel_names(channel_info_i)
            if list(output_channels_i) != list(channels):
                raise FileReaderError(
                    f"--Not all files have the same BT/BC output channels: {fname}\nCompare with: {mfiles[0]}"
                )

            if not channel_info_i.loc[channels, "bins"].equals(channel_info.loc[:, "bins"]):
                raise FileReaderError(
                    f"--Not all files include BT/BC channels with the same number of bins: {fname}\nCompare with: {mfiles[0]}"
                )

            shots_arr[k, :] = channel_info_i.loc[channels, "shots"].values

            if not lazy:
                sig = read_body(
                    buffer=buffer,
                    channel_info=channel_info_i,
                    output_channels=channels,
                    fname=fname,
                )
                eager_signal_blocks.append(sig)

        if lazy:
            sig_delayed = delayed(read_single_file_signal)(
                fname=fname,
                expected_output_channels=channels,
                expected_output_bins=bins,
                expected_file_format=file_format,
            )

            sig_block = da.from_delayed(
                sig_delayed,
                shape=(len(channels), len(bins_arr)),
                dtype=np.float64,
            )

            lazy_signal_blocks.append(sig_block)

        start_time_arr[k] = stime_i

        if stime_i >= etime_i:
            end_time_arr[k] = etime_i + timedelta(milliseconds=500)
        else:
            end_time_arr[k] = etime_i

    if lazy:
        sig_arr = da.stack(lazy_signal_blocks, axis=0)
    else:
        sig_arr = np.stack(eager_signal_blocks, axis=0)

    sig_raw = xr.DataArray(
        sig_arr,
        coords=[start_time_arr, channels, bins_arr],
        dims=["time", "channel", "bins"],
    ).astype(float)

    shots = xr.DataArray(
        shots_arr,
        coords=[start_time_arr, channels],
        dims=["time", "channel"],
    ).astype(float)

    tdata = np.array(
        [
            filename,
            datetimes_to_iso(start_time_arr),
            datetimes_to_iso(end_time_arr),
        ],
        dtype=object,
    )

    properties = ["filename", "start_time", "end_time"]

    time_info = pd.DataFrame(
        tdata.T,
        index=start_time_arr,
        columns=properties,
        dtype=object,
    )

    # Sort by time. Avoid copy() so lazy signal remains lazy.
    sig_raw = sig_raw.sortby("time")
    shots = shots.sortby("time")
    time_info = time_info.sort_index()

    # Convert bits to mV relying solely on the returned BT/BC channel_info.
    # This remains lazy when sig_raw is Dask-backed.
    sig_raw = unit_conv_bits_to_mV(
        signal=sig_raw,
        shots=shots,
        channel_info=channel_info,
    )

    if lazy:
        if chunks is None:
            time_chunks = min(5, max(10, sig_raw.sizes["time"]))
            channel_chunks = 5
            
            chunks = {
                "time": time_chunks,
                "channel": channel_chunks,
                "bins": -1,
            }

        sig_raw = sig_raw.chunk(chunks)

    return system_info, channel_info, time_info, sig_raw, shots
