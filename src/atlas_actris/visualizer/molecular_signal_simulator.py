#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simulate raw and range-corrected analog lidar molecular signals."""

import numpy as np

from utils import us_std
from arc_actris import arc
from typing import Optional
from scipy.integrate import cumulative_trapezoid

def get_molecular_profile(
    wavelength: float,
    sig_nr: float = 40.0,
    sig_d: float = 6.,
    alt_nr: float = 300.0,
    bsc_p_nr: float = 20.0,
    angstrom: float = 1.0,
    pldr: float = 0.2,
    alt_fr: float = 40000.0,
    aod: float = 0.7,
    alt_aerosol_end: float = 3000.,
    zenith_angle: float = 5.,
    n_altitude_points: int = 3000,
    channel_type: str = 't',
    channel_subtype: str = 'n',
    altitude_agl: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate baseline-corrected far-range analog signals.

    Parameters
    ----------
    sig_d
        Dark-signal baseline in mV.
    sig_nr
        Signal at the near-range peak in mV.
    wavelength
        Channel wavelength in nm.
    alt_nr
        Near-range peak height in m.
    bsc_p_nr
        Aerosol backscatter coefficient at 355 nm and the near-range peak,
        in Mm^-1 sr^-1.
    angstrom
        Aerosol Angstrom exponent.
    lidar_ratio
        Aerosol lidar ratio in sr. Retained for API compatibility; it is not
        used when ``aod`` is supplied directly between ``alt_nr`` and
        ``alt_aerosol_end``.
    pldr
        Particle linear depolarization ratio.
    alt_fr
        Far-range target molecular height in m.
    aod
        Vertical aerosol optical depth at 355 nm between ``alt_nr`` and
        ``alt_aerosol_end``. A value of zero represents an aerosol-free path
        between the two reference ranges.
    alt_aerosol_end
        Height above which there are no aerosols
    zenith_angle
        Zenith angle of the system
    n_altitude_points
        Number of points used for the internal atmospheric integration.
    altitude_agl
        Optional one-dimensional altitude-above-ground array in m on which
        the output profiles are requested. The calculations are performed on
        the internal grid and interpolated to these coordinates. Values below
        ``alt_nr`` or above ``alt_fr`` are returned as NaN.

    Returns
    -------
    z_range : np.ndarray
        Slant range in m on the requested output grid. When ``altitude_agl``
        is supplied, it has the same shape and bin ordering as that array;
        bins outside ``alt_nr`` to ``alt_fr`` are NaN. Empty for an
        unsupported channel.
    sig_an : np.ndarray
        Simulated baseline-corrected raw signal in same unit as sig_nr. When
        ``altitude_agl`` is supplied, it has the same shape and bin ordering
        as that array; bins outside ``alt_nr`` to ``alt_fr`` are NaN. Empty
        for an unsupported channel.
    sig_rc_an : np.ndarray
        Simulated baseline-corrected range-corrected signal.
        When ``altitude_agl`` is supplied, it has the same shape and bin
        ordering as that array; bins outside ``alt_nr`` to ``alt_fr`` are
        NaN. Empty for an unsupported channel.
    """


    vr_wavelength = None
    vr_mode = None
    
    if alt_fr <= alt_nr:
        raise ValueError("alt_fr must be greater than alt_nr.")
    if not alt_nr < alt_aerosol_end < alt_fr:
        raise ValueError(
            "alt_aerosol_end must lie between alt_nr and alt_fr."
        )
    if n_altitude_points < 2:
        raise ValueError("n_altitude_points must be at least 2.")
    if wavelength <= 0.0:
        raise ValueError("wavelength must be positive.")
    if pldr < 0.0:
        raise ValueError("pldr must be non-negative.")
    if not 0.0 <= zenith_angle < 90.0:
        raise ValueError(
            "zenith_angle must satisfy 0 <= zenith_angle < 90 degrees."
        )
    if aod < 0.0:
        raise ValueError(
            "The provided aod must be non-negative."
        )

    if altitude_agl is not None:
        altitude_agl = np.asarray(altitude_agl, dtype=float)

        if altitude_agl.ndim != 1:
            raise ValueError("altitude_agl must be a one-dimensional array.")

        if np.any(~np.isfinite(altitude_agl)):
            raise ValueError("altitude_agl must contain only finite values.")
        
    allowed_channel_types = {"t", "p", "c", "r", "v"}
    allowed_channel_subtypes = {"n", "o"}

    if channel_type not in allowed_channel_types:
        empty = np.array([], dtype=float)
        return empty, empty.copy(), empty.copy()
    
    if channel_type == "v" and channel_subtype not in allowed_channel_subtypes:
        empty = np.array([], dtype=float)
        return empty, empty.copy(), empty.copy()
        
    atmosphere = us_std.Atmosphere()
    boltzmann_constant = 1.380649e-23

    altitude = np.linspace(alt_nr, alt_fr, n_altitude_points)
    pressure = np.array([1e2 * atmosphere.pressure(value) for value in altitude])
    temperature = np.array([atmosphere.temperature(value) for value in altitude])
    number_density = pressure / (boltzmann_constant * temperature)

    cos_zenith = np.cos(np.deg2rad(zenith_angle))
    airmass_factor = 1.0 / cos_zenith
    
    z_range = altitude * airmass_factor
    
    end_index = int(
        np.argmin(np.abs(altitude - alt_aerosol_end))
    )
        
    number_density_ratio = number_density[end_index] / number_density[0]
    
    range_attenuation_ratio = (z_range[0] / z_range[end_index]) ** 2
    
    if channel_type == 'v':
        if channel_subtype == 'n':
            vr_mode = 'vibrational_raman_N2'
            vr_species = 'N2'
            species_fraction = 0.78084

        elif channel_subtype == 'o':
            vr_mode = 'vibrational_raman_O2'
            vr_species = 'O2'
            species_fraction = 0.20946

        raman_info = arc(
            incident_wavelength=wavelength,
            max_J=40,
            backscattering=True,
            mode=vr_mode,
        )
        vr_wavelength = raman_info.lamda_pol[vr_species]
        molecular_backscatter_cross_section_v = (
            raman_info.cross_section("main_line")
        )

    particle_backscatter_nr = bsc_p_nr * (wavelength / 355.0) ** (-angstrom)

    # ``aod`` is the vertical aerosol optical depth between ``alt_nr`` and
    # ``alt_aerosol_end`` at 355 nm. Convert it to the slant optical depth
    # and then scale it to the emitted and Raman-shifted wavelengths.
    delta_aod_up = (
        aod
        * airmass_factor
        * (wavelength / 355.0) ** (-angstrom)
    )

    if channel_type == 'v':
        delta_aod_dn = (
            aod
            * airmass_factor
            * (vr_wavelength / 355.0) ** (-angstrom)
        )
    else:
        delta_aod_dn = delta_aod_up

    rayleigh_bsc = arc(
        incident_wavelength=wavelength,
        max_J=40,
        backscattering=True,
        mode='rotational_raman',
    )

    rayleigh_sca = arc(
        incident_wavelength=wavelength,
        max_J=40,
        backscattering=False,
        mode='rotational_raman',
    )

    if channel_type == 'v':
        rayleigh_sca_vr = arc(
            incident_wavelength=vr_wavelength,
            max_J=40,
            backscattering=False,
            mode=vr_mode,
        )

    molecular_extinction_cross_section = rayleigh_sca.cross_section("full")
    
    molecular_backscatter_cross_section_t = rayleigh_bsc.cross_section("main_line")
    molecular_backscatter_cross_section_r = rayleigh_bsc.cross_section("O")
    molecular_depolarization_ratio = rayleigh_bsc.mldr("main_line")

    molecular_backscatter_cross_section_p = (
        molecular_backscatter_cross_section_t 
        / (1.0 + molecular_depolarization_ratio)
    )
    molecular_backscatter_cross_section_s = (
        molecular_backscatter_cross_section_t
        * molecular_depolarization_ratio
        / (1.0 + molecular_depolarization_ratio)
    )
    molecular_extinction = molecular_extinction_cross_section * number_density

    mod_cumulative_up = cumulative_trapezoid(
        molecular_extinction, 
        x=z_range, 
        initial=0.0
        )
    
    if channel_type == 'v':
        molecular_extinction_cross_section_vr = (
            rayleigh_sca_vr.cross_section("full")
        )
        molecular_extinction_v = (
            molecular_extinction_cross_section_vr
            * number_density
        )

        mod_cumulative_dn = cumulative_trapezoid(
            molecular_extinction_v,
            x=z_range,
            initial=0.0
            )
    else:
        mod_cumulative_dn = mod_cumulative_up
    
   
    
    if channel_type == 'p':
        particle_backscatter_p = particle_backscatter_nr / (1.0 + pldr)
        molecular_backscatter_p = molecular_backscatter_cross_section_p * number_density       
        scattering_ratio_nr = 1.0 + (
            1e-6 * particle_backscatter_p / molecular_backscatter_p[0]
        )
        molecular_backscatter = molecular_backscatter_cross_section_p * number_density

    elif channel_type == 'c':
        molecular_backscatter_s = molecular_backscatter_cross_section_s * number_density
        particle_backscatter_s = pldr * particle_backscatter_nr / (1.0 + pldr)
        scattering_ratio_nr = 1.0 + (
            1e-6 * particle_backscatter_s / molecular_backscatter_s[0]
        )
        molecular_backscatter = molecular_backscatter_cross_section_s * number_density

    elif channel_type == 't':
        molecular_backscatter_t = molecular_backscatter_cross_section_t * number_density
        scattering_ratio_nr = 1.0 + (
            1e-6 * particle_backscatter_nr / molecular_backscatter_t[0]
            )
        molecular_backscatter = molecular_backscatter_cross_section_t * number_density

    elif channel_type == 'v':
        scattering_ratio_nr = 1.0
        species_number_density = species_fraction * number_density
        molecular_backscatter = (
            molecular_backscatter_cross_section_v
            * species_number_density
        )

    elif channel_type == 'r':
        scattering_ratio_nr = 1.0
        molecular_backscatter = molecular_backscatter_cross_section_r * number_density

    delta_mod_up = mod_cumulative_up[end_index]

    if channel_type == 'v':
        delta_mod_dn = mod_cumulative_dn[end_index]
    else:
        delta_mod_dn = delta_mod_up

    od_up_term = np.exp(-delta_aod_up - delta_mod_up)
    
    if channel_type == 'v':
        od_dn_term = np.exp(-delta_aod_dn - delta_mod_dn)
    else:
        od_dn_term = od_up_term
        
    sig_an_end = (
        (sig_nr - sig_d)
        * range_attenuation_ratio
        / scattering_ratio_nr
        * number_density_ratio
        * od_up_term
        * od_dn_term
    )
    
    sig_mol = (
        molecular_backscatter
        * np.exp(
            -mod_cumulative_up
            -mod_cumulative_dn
        )
        / z_range**2
    )

    sig_an = sig_an_end * sig_mol / sig_mol[end_index]
    
    sig_rc_an = sig_an * z_range**2

    if altitude_agl is not None:
        valid_output = (
            (altitude_agl >= alt_nr)
            & (altitude_agl <= alt_fr)
        )

        output_z_range = np.full(altitude_agl.shape, np.nan, dtype=float)
        output_sig_an = np.full(altitude_agl.shape, np.nan, dtype=float)
        output_sig_rc_an = np.full(altitude_agl.shape, np.nan, dtype=float)

        output_z_range[valid_output] = (
            altitude_agl[valid_output] * airmass_factor
        )
        output_sig_an[valid_output] = np.interp(
            altitude_agl[valid_output],
            altitude,
            sig_an,
        )
        output_sig_rc_an[valid_output] = np.interp(
            altitude_agl[valid_output],
            altitude,
            sig_rc_an,
        )

        return output_z_range, output_sig_an, output_sig_rc_an

    return z_range, sig_an, sig_rc_an
