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
import contextlib
import numpy as np
import xarray as xr

from arc_actris import arc

from typing import Any, Dict
from utils.dataarray_utils import shallow_copy
from utils.printouts import print_entry

def compute_molecular_calculations(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)

    height_asl = output_data["height_asl"]
    meteo = output_data["meteo"]
    channel_info = output_data["channel_info"]

    temperature_scale = np.arange(180.0, 330.0, 10.0)

    opto_parameters = [
        "c_ext_f", 
        "c_ext_b", 
        "c_bsc", 
        "ext_f", 
        "ext_b", 
        "bsc", 
        "OD_f", 
        "OD_b", 
        "atten_bsc"
        ]
    
    channel_chunks = 5

    molecular = {}

    if not meteo:
        print_entry("No meteorological profiles found. Molecular calculations not performed")
        return output_data

    for key in meteo.keys():

        channels = channel_info[key].channel.values

        meteo_on_bins = (
            meteo[key]
            .interp(height_asl=height_asl[key])
            .transpose("channel", "atmo_parameters", "bins")
        )

        emitted_wavelength = channel_info[key].loc["emitted_wavelength"]
        detected_wavelength = channel_info[key].loc["detected_wavelength"]
        channel_bandwidth = channel_info[key].loc["channel_bandwidth"]

        molecular_channels = []

        for ch in channels:

            c_ext_f, c_ext_b, c_bsc = get_optical_parameters(
                ch=ch,
                temperature_scale=temperature_scale,
                emitted_wavelength=emitted_wavelength,
                detected_wavelength=detected_wavelength,
                channel_bandwidth=channel_bandwidth,
            )

            T_on_bins = (
                meteo_on_bins
                .sel(channel=ch, atmo_parameters="T")
                .reset_coords(drop=True)
            )

            N_on_bins = (
                meteo_on_bins
                .sel(channel=ch, atmo_parameters="N")
                .reset_coords(drop=True)
            )

            z_on_bins = (
                height_asl[key]
                .sel(channel=ch)
                .reset_coords(drop=True)
            )

            c_ext_f_on_bins = (
                c_ext_f
                .interp(T=T_on_bins)
                .reset_coords(drop=True)
            )
            
            c_ext_b_on_bins = (
                c_ext_b
                .interp(T=T_on_bins)
                .reset_coords(drop=True)
            )

            c_bsc_on_bins = (
                c_bsc
                .interp(T=T_on_bins)
                .reset_coords(drop=True)
            )

            ext_f_on_bins = (
                N_on_bins * c_ext_f_on_bins
            ).reset_coords(drop=True)

            ext_b_on_bins = (
                N_on_bins * c_ext_b_on_bins
            ).reset_coords(drop=True)

            bsc_on_bins = (
                N_on_bins * c_bsc_on_bins
            ).reset_coords(drop=True)

            # Cumulative optical depth:
            # OD(z_i) = integral from first bin to z_i of ext(z) dz
            dz = z_on_bins.diff("bins", label="upper")

            ext_upper_f = ext_f_on_bins.isel(bins=slice(1, None))
            ext_upper_b = ext_b_on_bins.isel(bins=slice(1, None))

            ext_lower_f = (
                ext_f_on_bins
                .isel(bins=slice(None, -1))
                .assign_coords(bins=ext_upper_f.bins)
            )
            
            ext_lower_b = (
                ext_b_on_bins
                .isel(bins=slice(None, -1))
                .assign_coords(bins=ext_upper_b.bins)
            )

            dz_f = dz.assign_coords(bins=ext_upper_f.bins)
            dz_b = dz.assign_coords(bins=ext_upper_b.bins)

            ext_mid_f = 0.5 * (ext_upper_f + ext_lower_f)
            ext_mid_b = 0.5 * (ext_upper_b + ext_lower_b)

            od_intervals_f = (ext_mid_f * dz_f).cumsum("bins")
            od_intervals_b = (ext_mid_b * dz_b).cumsum("bins")

            zero_f = xr.zeros_like(ext_f_on_bins.isel(bins=0))
            zero_b = xr.zeros_like(ext_b_on_bins.isel(bins=0))

            od_f_on_bins = xr.concat(
                [zero_f, od_intervals_f],
                dim="bins",
            ).assign_coords(
                bins=ext_f_on_bins.bins
            ).reset_coords(drop=True)

            od_b_on_bins = xr.concat(
                [zero_b, od_intervals_b],
                dim="bins",
            ).assign_coords(
                bins=ext_b_on_bins.bins
            ).reset_coords(drop=True)

            atten_bsc_on_bins = (
                bsc_on_bins * np.exp(-(od_f_on_bins + od_b_on_bins))
            ).reset_coords(drop=True)

            molecular_ch = xr.concat(
                [
                    c_ext_f_on_bins,
                    c_ext_b_on_bins,
                    c_bsc_on_bins,
                    ext_f_on_bins,
                    ext_b_on_bins,
                    bsc_on_bins,
                    od_f_on_bins,
                    od_b_on_bins,
                    atten_bsc_on_bins,
                ],
                dim=xr.DataArray(
                    opto_parameters,
                    dims="opto_parameters",
                    name="atmo_parameters",
                ),
            )

            molecular_ch = molecular_ch.expand_dims(channel=[ch])

            molecular_channels.append(molecular_ch)

        molecular[key] = (
            xr.concat(molecular_channels, dim="channel")
            .transpose("channel", "opto_parameters", "bins")
            .chunk({
                "channel": channel_chunks,
                "opto_parameters": -1,
                "bins": -1,
            })
        )

        print_entry(f"Molecular calculations for {key} QA test complete")

    output_data["molecular"] = molecular

    return output_data

def get_optical_parameters(ch, temperature_scale, emitted_wavelength, 
                           detected_wavelength, channel_bandwidth):
    
    c_ext_f = np.nan * np.zeros_like(temperature_scale)
    c_ext_b = np.nan * np.zeros_like(temperature_scale)
    c_bsc = np.nan * np.zeros_like(temperature_scale)
        
    forward_wavelength = float(emitted_wavelength.loc[ch].values)
    backward_wavelength = float(detected_wavelength.loc[ch].values)
    
    filter_parameters = {
        'central_wavelength':float(detected_wavelength.loc[ch].values),
        'bandwidth':float(channel_bandwidth.loc[ch].values)
        }
     
    if ch[5] == 'v' and ch[7] == 'n':
        mode = 'vibrational_raman_N2'
        filter_parameters['transmission_shape'] = 'Gaussian'
        normalize = False
        
    elif ch[5] == 'v' and ch[7] == 'o':
        mode = 'vibrational_raman_O2'
        filter_parameters['transmission_shape'] = 'Gaussian'
        normalize = False
        
    elif ch[5] in ['p','c','t','r']:
        mode = 'rotational_raman'
        
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
                
                rre_f = arc(
                        incident_wavelength = forward_wavelength, 
                        temperature = T,
                        max_J = 40, 
                        backscattering = False,
                        mode = "rotational_raman",
                        )
                
                rre_b = arc(
                        incident_wavelength = backward_wavelength, 
                        temperature = T,
                        max_J = 40, 
                        backscattering = False,
                        mode = "rotational_raman",
                        )
                
                if ch[5] in ['p','c','t','r']:
                    rrb = arc(
                        incident_wavelength = forward_wavelength, 
                        temperature = T,
                        max_J = 40, 
                        backscattering = True,
                        mode = mode,
                        filter_parameters = filter_parameters,
                        )
                    
                    c_bsc[i] = rrb.cross_section(cross_section_type = 'full', normalize = normalize)

                    mldr = rrb.mldr(mldr_type = 'full')

                    if ch[5] == 'p':
                        c_bsc[i] = 1. / (1. + mldr) * c_bsc[i]
                    elif ch[5] == 'c':
                        c_bsc[i] = mldr / (1. + mldr) * c_bsc[i]

                elif ch[5] in ['v']:
                    rrb = arc(
                        incident_wavelength = forward_wavelength, 
                        temperature = T,
                        max_J = 40, 
                        backscattering = True,
                        mode = mode,
                        filter_parameters = filter_parameters,
                        )
                    
                    c_bsc[i] = rrb.cross_section(cross_section_type = 'full', normalize = normalize)
                    
                else:
                    c_bsc[i] = np.nan
                
                c_ext_f[i] = rre_f.cross_section(cross_section_type = 'full')
                c_ext_b[i] = rre_b.cross_section(cross_section_type = 'full')
    
    c_ext_f = xr.DataArray(c_ext_f, dims = ['T'], coords = [temperature_scale])
    c_ext_b = xr.DataArray(c_ext_b, dims = ['T'], coords = [temperature_scale])
    c_bsc = xr.DataArray(c_bsc, dims = ['T'], coords = [temperature_scale])
    
    return(c_ext_f, c_ext_b, c_bsc)