"""
@author: N.Siomos & P. Paschou

Functions for the retieval of the products:
- bsc coefficient (Klett)
"""
import numpy as np
from tools import prod_xarray
from elastic_retrievals import klett_backscatter_aerosol


old = np.seterr(divide='ignore',invalid='ignore')

def klett_bsc(prod, height_levels, info_prod, sig_ds, lr_const, sr_ref, ref_height, ref_hwindow):
    '''
    Retrieves the backscatter coefficient using the Klett inversion method

    Returns
    -------
    prod : x-array with 2 dimensions [(TIME), product_ids, bins)
        The prod x-array with the new retrieved products.
   
    height_levels : x-array with 2 dimensions [(TIME), product_ids, bins)
        The x-array with the corresponding height levels for each product.
    
    info_prod : DATAFRAME WITH OBJECT_TYPE VALUES
        The info about the new retrieved products.
    '''
   

    sig_wave = sig_ds.Detected_Wavelength
    
    # mask to fetch only the elastic channels
    mask_elastic = ((sig_wave >= 354.4) & (sig_wave <= 355.4)) \
            | ((sig_wave >= 531.4) & (sig_wave <= 533.4))\
            | ((sig_wave >= 1063.4) & (sig_wave <= 1065.4))
    
    
    # -- Retrieval --
    # Signal
    signal = sig_ds.Range_Corrected_Signals # DataArray 2-D (channel, bins)
        
    for j in signal.channel[mask_elastic.values].values:#range(signal.channel.values[mask_elastic.values].size):    
        
        
        # elastic signal id for the retrieval
        sig_id = j #signal.channel.values[j]
        
        ind_ch = dict(channel = sig_id)
        
        # wavelength
        wave_val = sig_wave.loc[sig_id].values #str(int(info.wave.loc[sig_id]))
        
        # Height and vertical resolution
        height = sig_ds.Height_levels.loc[ind_ch]
        step = sig_ds.Raw_Data_Range_Resolution.loc[sig_id].values
        
        pd_id = f'parBsc_klett_{sig_id}'
        # mol_id = sig_ds.channel.values[j] # molecular ids are same as channels
        
        # Create/append prod (DataArray 2-D; dims are channel, bins) and info_prod (Dataframe)
        prod, height_levels, info_prod = prod_xarray(prod, height_levels, info_prod, [pd_id], signal)
        
        # index of reference height and half window bins number -> ref_idx - ref_hwin, ref_idx + ref_hwin
        ref_idx = np.argmin(np.abs(ref_height - height.values))
        ref_hwin = int(ref_hwindow / step)
        
        
        # Molecular profiles (ext, bsc, lr)
        a_mol = sig_ds.Extinction_Coefficient_Forward.loc[ind_ch].values
        b_mol = sig_ds.Backscatter_Coefficient.loc[ind_ch].values
        lr_mol = a_mol/b_mol
                
        # for t in range(prod.time.size):
            # ind_t = dict(time = t)
        sig_tmp = signal.loc[ind_ch].values
        if not all(np.isnan(sig_tmp)):
            #Retrieves the product only if the signal in the time frame exists
            bsc = klett_backscatter_aerosol(sig_tmp, lr_const, b_mol, ref_idx, 
                                            ref_hwin, sr_ref, step, lr_mol)  
            
            prod.loc[dict(product = pd_id)] = bsc
            height_levels.loc[dict(product = pd_id)] = height.values
            
            info_prod.loc[pd_id, 'wave'] = wave_val
            info_prod.loc[pd_id, 'polarization'] = 'linear'
            info_prod.loc[pd_id, 'type'] = 'aerosol_backscatter_coefficient'
            info_prod.loc[pd_id, 'type_short'] = 'bsc'
            info_prod.loc[pd_id, 'method'] = 'klett'
            info_prod.loc[pd_id, 'units'] = r'$m^{-1}sr^{-1}$'
            info_prod.loc[pd_id, 'process'] = 'retrieved'
            info_prod.loc[pd_id, 'retrieving_info'] = f'Fixed lidar ratio: {lr_const}sr; Reference region: {ref_height-ref_hwindow} - {ref_height+ref_hwindow} m' 
    
    return(prod, height_levels, info_prod)

