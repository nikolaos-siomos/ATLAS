"""
@author: N. Siomos & P. Paschou
"""

import numpy as np

def range_to_alt(x_range, ground, SZA):
    alt = ground + x_range*np.cos(SZA*np.pi/180.)
    return(alt)

def alt_to_range(altitude, ground, SZA):    
    x_range = (altitude - ground)/np.cos(SZA*np.pi/180.)
    x_range = np.round(x_range, decimals = 2)
    return(x_range)
    
def profiles_interp(alt, prs, tem, hum, scale):
    
    ulim = np.min([scale[-1], alt[-1]])

    scale = scale[scale <= ulim]

    prs = np.interp(scale, alt , prs)
    tem = np.interp(scale, alt, tem)
    hum = np.interp(scale, alt, hum)
    
    return(scale, prs, tem, hum)

def trim_range(signal, m_range):
        
    c_range = np.intersect1d(signal.range.values, m_range, 
                             assume_unique = True)
    
    rng_cursor = dict(range = c_range)
    signal = signal.loc[rng_cursor]
    
    return(c_range, signal)

def bin_conversions(info, sig, cfg):
    
    gnd_alt = float(cfg.lidar.altitude)
    SZA     = float(cfg.angles.zenith)

    x_range = sig.range.values
    
    alt = range_to_alt(x_range, gnd_alt, SZA)
    alt_resol = range_to_alt(info.resol,0.,SZA)
    
    for ch in cfg.ch_map:
        ref_idx = np.abs(cfg.ref_h.loc[ch] - alt)
        full_ovl_idx = np.abs(cfg.full_ovl.loc[ch] - x_range)
        info.loc[ch,'ref_idx'] = np.argmin(ref_idx)
        info.loc[ch,'full_ovl_idx'] = np.argmin(full_ovl_idx)
        info.loc[ch,'ref_hwin_rng'] = cfg.ref_hwin.loc[ch]/info.resol.loc[ch]
        info.loc[ch,'ref_hwin'] = cfg.ref_hwin.loc[ch]/alt_resol.loc[ch]
    
    # convert the half window bins dtype from float to int
    info.ref_hwin = info.ref_hwin.astype(int)
    info.ref_hwin_rng = info.ref_hwin_rng.astype(int)
    
    return(info)
