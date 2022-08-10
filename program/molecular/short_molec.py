"""
@author: N. Siomos & P. Paschou
"""

import os
from molecular import height_scales
from molecular.rayleigh_atm import elastic, depol, raman, \
    atmospheric_profiles
import xarray as xr
import pandas as pd
import numpy as np

def short_molec(signal, cfg, molec_map):
    
    print('-----------------------------------------')
    print('Calculating the molecular profiles...')
    print('-----------------------------------------')
    
    dirs = cfg.directories
    imol = cfg.imol
    gnd_alt = float(cfg.lidar.altitude)
    gnd_us = cfg.us
    SZA = float(cfg.angles.zenith)
            
    sig_range = signal.range.values
    dir_inp = os.path.join(dirs.input, 'atmosphere')
    
    alt_lr = height_scales.range_to_alt(sig_range, gnd_alt, SZA)
    
    alt_rs, P, T, RH, N, N_H2O = \
        atmospheric_profiles(dir_inp, imol, alt_lr, gnd_us, gnd_alt, SZA)
    
    m_range = height_scales.alt_to_range(alt_rs, gnd_alt, SZA)

    c_range, signal = height_scales.trim_range(signal, m_range)
    
    a_mol_el_e, b_mol_el, attn_el, map_el = \
        elastic(c_range, P, T, molec_map.elastic)
    a_mol_el_r = a_mol_el_e.copy()
    
    a_mol_dp_e, b_mol_dp, attn_dp, map_dp = depol(c_range, P, T, molec_map.depol)
    a_mol_dp_r = a_mol_dp_e.copy()

    a_mol_rm_e, a_mol_rm_r, b_mol_rm, attn_rm, map_rm = \
        raman(c_range, P, T, molec_map.raman)

    # Merge all xarrays
    a_mol_e = xr.concat([a_mol_el_e, a_mol_dp_e, a_mol_rm_e], dim = 'id')
    a_mol_r = xr.concat([a_mol_el_r, a_mol_dp_r, a_mol_rm_r], dim = 'id')
    b_mol = xr.concat([b_mol_el, b_mol_dp, b_mol_rm], dim = 'id')
    attn = xr.concat([attn_el, attn_dp, attn_rm], dim = 'id')
    
    optics = xr.concat([a_mol_e, a_mol_r, b_mol, attn], dim = 'types')
    optics = \
        optics.assign_coords({"types": ['ext_forth', 'ext_back', 'bsc', 'atten_bsc']})
        
    
    atm_ids = ['P', 'T', 'RH', 'N', 'N_H2O']
    atm = xr.DataArray(np.vstack((P,T,RH,N,N_H2O)), 
                       coords=[atm_ids, c_range], dims=['types', 'range'])
    
    molec_pack = dict()
    
    molec_pack['optics'] = optics
    molec_pack['atm'] = atm
    
    molec_map = pd.concat([map_el, map_dp, map_rm]).astype(float)
    
    print('-- Calculation of molecular properties complete!')
    print('-----------------------------------------')
    print('')

    return(signal, molec_pack, molec_map)