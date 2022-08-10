"""
@author: N. Siomos & P. Paschou

Determines the reference height and makes rayleigh fit plots
"""
import os, re
import pandas as pd
import numpy as np
from molecular import height_scales
from .rayleigh_fit import calc_window, update_cfg_ref 
from helper_functions.helper_functions import find_array_index
from readers.read_maps import plt_options
from readers.read_config import comma_split
from plotting.plot_rfit import plot_rfit_one, plot_rfit_two


def rayleigh_fit(signal, info, attn, cfg, maps_dir):
    
    '''
    Parameters
    ----------
    signal : 3D xarray [time, channel, range]
        The signals of the channels.
    info : Dataframe with mixed object/float/integer type values
        Contains info about the channels in the signal xarray.
    attn : 2D xarray [id, range]
        The attenuated backscatter coefficient for each channel (wavelength oriented).
    cfg : class with atrributes of multiple dtypes
        Contains all the atributes given in the config_file.

    '''
    
    print('-----------------------------------------')
    print('Proceeding to the Rayleigh fit test...')
    print('-----------------------------------------')
    
    # Extract configuration parameters      
    dirs      = cfg.directories
    vis_map   = cfg.vis_map # info.index.values
    gnd_alt   = float(cfg.lidar.altitude)
    SZA       = float(cfg.angles.zenith)
    molec_map = cfg.molec_map.values
    inter     = cfg.inter
    ref_idx   = info.ref_idx.loc[vis_map]
    ref_hwin  = info.ref_hwin.loc[vis_map]
    iray_fit  = cfg.iray_fit
    stime_clb = cfg.start_t_clb
    etime_clb = cfg.end_t_clb
    
    # Extract plotting options from plot_options.ini
    options = plt_options(maps_dir, 'plot_options.ini','Rayleigh_Fit').rfit
    plt_dif = int(options.plt_dif)

    # Create range array from range array
    alt = height_scales.range_to_alt(signal.range.values, gnd_alt, SZA)
    alt_dict = {'altitude':alt}    

    ch_dict =  {'channel':vis_map} 
    
    # Replace ids with signal channels
    temp = attn.sel(id = molec_map)
    sig_mol = temp.swap_dims({'range':'altitude'}).assign_coords(alt_dict).swap_dims({'id':'channel'}).assign_coords(ch_dict)   

    # Get only signals that will be normalized
    temp = signal.copy().loc[:,vis_map,:]
    sig_fit = temp.swap_dims({'range':'altitude'}).assign_coords(alt_dict)

    # Create the plotting directory if not exists    
    dir_plt = os.path.join(dirs.plots, 'pre_processing')
    os.makedirs(dir_plt, exist_ok=True)
    
    print('-- Starting Rayleigh fit Test')
    
    find_ref = 1
    while find_ref != 0: #loop to run for every time I set new ref_h
        
        ref_sig_mol, ref_sig_fit, alt_lims = \
            calc_window(vis_map, ref_idx, ref_hwin, sig_fit, sig_mol)
        
        # Calculate the normalization/calibration factor 
        norm_factor = ref_sig_mol/ref_sig_fit # dims = [channel, time]
                
        # adjust the norm_factor if specific timeframe is given
        if len(stime_clb)>0 and len(etime_clb)>0:
            t1 = np.datetime64(stime_clb)
            t2 = np.datetime64(etime_clb)                    
            norm_factor = norm_factor.loc[:,t1:t2].mean(dim='time')

        # Normalize the signal
        sig_nrm = sig_fit*norm_factor
        print('-- Normalization of rc signals(see vis_map) in the reference height complete!')


        if iray_fit == 1:
            print('-- Generating Rayleigh fit plots')
            
            #Create the plots
            if plt_dif == 0:
                # one block of plots in figure. Only the rayleigh fit
                plot_rfit_one(dir_plt, sig_nrm, sig_mol, alt_lims, cfg, options)            
            else:
                # Two block of plots in figure. The rayleigh fit and the relative difference
                plot_rfit_two(dir_plt, sig_nrm, sig_mol, alt_lims, cfg, options)            
            
            if inter == 1:
            # Interactive option on
            #After the last rayleigh fit plot, ask the user to decide to set new height or not
                find_ref, info, cfg = ask_user(info, sig_nrm, cfg)
                # re-assign the vars
                ref_idx   = info.ref_idx
                
            else:
            # Interactive option off
            # the reference height is the input value from the user in the config file
                find_ref = 0
        else:
            find_ref = 0
                
    # append new channels from vis_map in the reh_h and ref_hwin
    cfg = update_cfg_ref(cfg, ref_idx, alt)
        
    print('-- Rayleigh fit test complete!') 
    print('-----------------------------------------')
    print('')
    
    return(sig_nrm, sig_mol, info)

def ask_user(info, sig, cfg):
    '''
    Asks user if he want to set new reference height per wavelength of channels
    Returns the find_ref indicator, the info dataframe, and the cfg parameters 
    '''
    
    vis_map = cfg.vis_map
    
    print('')
    print('Set new reference height? (1: yes / 0: no)')

    find_ref = input()
    
    # check for non number input from user
    if find_ref.strip().isalpha():
        raise TypeError('Warning! User\'s given answer is not 0 or 1.\n Programm stopped!')

    else:    
        find_ref = int(find_ref)
    
    # Ask user for new ref_h
    if find_ref == 1:
        
        sep = ', '
        
        print('Set reference height (in meters) for each channel with comma separation:')
        print(f'channels: {sep.join(vis_map)}')
        temp_ref_h = input('height: ')
        
        # Split the user input val and create pandas Series with the ref_h per wl
        #temp_ref_h = comma_split(temp_ref_h, float)
        temp_ref_h = pd.Series(comma_split(temp_ref_h, float),
                               index = vis_map)

        # assign the new ref_height in the right channels per wl          
        for ch in vis_map:
            info.ref_idx.loc[ch] = find_array_index(sig.altitude.values, 
                                                    temp_ref_h.loc[ch], 
                                                    float)

    print('')        

    return(find_ref, info, cfg)