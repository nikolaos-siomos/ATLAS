"""
@authors: N. Siomos, P. Paschou, Ioannis Binietoglou 
based on SULA project (https://react-gitlab.space.noa.gr/ReACT/eve/data-processing)
and also on https://gitlab.com/ioannis_binietoglou/lidar-processing/

Processing routines for signals 

=================================
Signal in 3D xarray dataset with dimensions [time, channel, bin/range]

Fucntions
 -- average_by_time: Average signals across the timeframes
 -- background_calculation: Calculates the solar background per timeframe and channel 
 -- background_correction: Performs the background correction on signals
 -- dark_correction: Removes the dark signals from the normal ananlog signals
 -- dead time correction: Performs the dead time correction onphoton channels
 -- detect_saturation: Identifies regions where signals are saturated
 -- height_calculation: Calculates the height above the lidar values per bin and channel
 -- range calculation: Calculates the range above the lidar values per bin and channel
 -- range_correction: Performs the range correction on signals
 -- smoothing: Smooths the signals (sliding average)
 -- trigger_correction: Perform the trigger correction per channel
 -- trim_vertically: Trim channels up to a maximum altitude
 -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels

"""

import io
import copy
import contextlib
import numpy as np
import xarray as xr

from arc_actris import arc

from typing import Any, Dict
from utils.dataarray_utils import shallow_copy
from helper_functions.printouts import print_entry

def compute_molecular_calculations(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)

    height_asl = output_data["height_asl"]
    meteo = output_data["meteo"]
    channel_info = output_data["channel_info"]

    temperature_scale = np.arange(180.0, 330.0, 10.0)

    parameters = ["c_ext", "c_bsc", "ext", "bsc", "OD", "atten_bsc"]
    bin_chunk_size = 4096

    molecular = {}

    if not meteo:
        print_entry("No meteorological profiles found. Molecular calculations not performed")
        return output_data

    for key in meteo.keys():

        channels = channel_info[key].channel.values

        meteo_on_bins = (
            meteo[key]
            .interp(height_asl=height_asl[key])
            .transpose("channel", "parameters", "bins")
        )

        emitted_wavelength = channel_info[key].loc["emitted_wavelength"]
        detected_wavelength = channel_info[key].loc["detected_wavelength"]
        channel_bandwidth = channel_info[key].loc["channel_bandwidth"]

        molecular_channels = []

        for ch in channels:

            c_ext, c_bsc = get_optical_parameters(
                ch=ch,
                temperature_scale=temperature_scale,
                emitted_wavelength=emitted_wavelength,
                detected_wavelength=detected_wavelength,
                channel_bandwidth=channel_bandwidth,
            )

            T_on_bins = (
                meteo_on_bins
                .sel(channel=ch, parameters="T")
                .reset_coords(drop=True)
            )

            N_on_bins = (
                meteo_on_bins
                .sel(channel=ch, parameters="N")
                .reset_coords(drop=True)
            )

            z_on_bins = (
                height_asl[key]
                .sel(channel=ch)
                .reset_coords(drop=True)
            )

            c_ext_on_bins = (
                c_ext
                .interp(T=T_on_bins)
                .reset_coords(drop=True)
            )

            c_bsc_on_bins = (
                c_bsc
                .interp(T=T_on_bins)
                .reset_coords(drop=True)
            )

            ext_on_bins = (
                N_on_bins * c_ext_on_bins
            ).reset_coords(drop=True)

            bsc_on_bins = (
                N_on_bins * c_bsc_on_bins
            ).reset_coords(drop=True)

            # Cumulative optical depth:
            # OD(z_i) = integral from first bin to z_i of ext(z) dz
            dz = z_on_bins.diff("bins", label="upper")

            ext_upper = ext_on_bins.isel(bins=slice(1, None))

            ext_lower = (
                ext_on_bins
                .isel(bins=slice(None, -1))
                .assign_coords(bins=ext_upper.bins)
            )

            dz = dz.assign_coords(bins=ext_upper.bins)

            ext_mid = 0.5 * (ext_upper + ext_lower)

            od_intervals = (ext_mid * dz).cumsum("bins")

            zero = xr.zeros_like(ext_on_bins.isel(bins=0))

            od_on_bins = xr.concat(
                [zero, od_intervals],
                dim="bins",
            ).assign_coords(
                bins=ext_on_bins.bins
            ).reset_coords(drop=True)

            atten_bsc_on_bins = (
                bsc_on_bins * np.exp(-2.0 * od_on_bins)
            ).reset_coords(drop=True)

            molecular_ch = xr.concat(
                [
                    c_ext_on_bins,
                    c_bsc_on_bins,
                    ext_on_bins,
                    bsc_on_bins,
                    od_on_bins,
                    atten_bsc_on_bins,
                ],
                dim=xr.DataArray(
                    parameters,
                    dims="parameters",
                    name="parameters",
                ),
            )

            molecular_ch = molecular_ch.expand_dims(channel=[ch])

            molecular_channels.append(molecular_ch)

        molecular[key] = (
            xr.concat(molecular_channels, dim="channel")
            .transpose("channel", "parameters", "bins")
            .chunk({
                "channel": 1,
                "parameters": -1,
                "bins": bin_chunk_size,
            })
        )

        print_entry(f"Molecular calculations for {key} QA test complete")

    output_data["molecular"] = molecular

    return output_data

def get_optical_parameters(ch, temperature_scale, emitted_wavelength, 
                           detected_wavelength, channel_bandwidth):
    
    c_ext = np.nan * np.zeros_like(temperature_scale)
    c_bsc = np.nan * np.zeros_like(temperature_scale)
        
    laser_wavelength = float(emitted_wavelength.loc[ch].values)
    
    filter_parameters = {
        'central_wavelength':float(detected_wavelength.loc[ch].values),
        'bandwidth':float(channel_bandwidth.loc[ch].values)
        }
     
    if ch[5] == 'v' and ch[7] == 'n':
        mode = 'vibrational_raman_N2'
        filter_parameters['transmission_shape'] = 'Gaussian'
        normalize = False
        
        rrv = arc(
            laser_wavelength, 
            max_J = 40, 
            backscattering = True,
            mode = mode,
            )
        
        incident_wavelength = rrv.lamda_pol['N2']
        
    elif ch[5] == 'v' and ch[7] == 'o':
        mode = 'vibrational_raman_O2'
        filter_parameters['transmission_shape'] = 'Gaussian'
        normalize = False
        
        with contextlib.redirect_stdout(io.StringIO()):

            rrv = arc(
                laser_wavelength, 
                max_J = 40, 
                backscattering = True,
                mode = mode,
                )
            
            incident_wavelength = rrv.lamda_pol['O2']

        
    elif ch[5] in ['p','c','t','r']:
        mode = 'rotational_raman'
        incident_wavelength = laser_wavelength
        
        if ch[5] == 'r':
            filter_parameters['transmission_shape'] = 'Tophat'
            normalize = False
            
        else:
            filter_parameters['transmission_shape'] = 'Gaussian'
            normalize = True
            
    else: 
        mode = 'not_applicable'
        
    if mode != 'not_applicable':
        for i, T in enumerate(temperature_scale):
            with contextlib.redirect_stdout(io.StringIO()):

                rrb = arc(
                    incident_wavelength, 
                    temperature = T,
                    max_J = 40, 
                    backscattering = True,
                    mode = mode,
                    filter_parameters = filter_parameters
                    )
                
                rre = arc(
                    incident_wavelength, 
                    temperature = T,
                    max_J = 40, 
                    backscattering = False,
                    mode = "rotational_raman",
                    )
            
            with contextlib.redirect_stdout(io.StringIO()):
                c_ext[i] = rre.cross_section(cross_section_type = 'full')
                c_bsc[i] = rrb.cross_section(cross_section_type = 'full', normalize = normalize)
            
            if ch[5] == 'p':
                with contextlib.redirect_stdout(io.StringIO()):
                    mldr = rrb.mldr(mldr_type = 'full')
                c_bsc[i] = 1. / (1. + mldr) * c_bsc[i]
            elif ch[5] == 'c':
                with contextlib.redirect_stdout(io.StringIO()):
                    mldr = rrb.mldr(mldr_type = 'full')
                c_bsc[i] = mldr / (1. + mldr) * c_bsc[i]
    
    c_ext = xr.DataArray(c_ext, dims = ['T'], coords = [temperature_scale])
    c_bsc = xr.DataArray(c_bsc, dims = ['T'], coords = [temperature_scale])
    
    return(c_ext, c_bsc)