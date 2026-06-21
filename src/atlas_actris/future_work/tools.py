# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 12:45:55 2025

@author: Peristera
"""
import numpy as np
import xarray as xr
import pandas as pd

def prod_xarray(prod, height_levels, info_prod, prod_id, signal=[]):
# Create an xarray or append new indexes
    
    if len(prod) > 0: # Append the new indexes of products
        add_idx = list(prod_id)#prod_map.prod_id.values)
        existing_idx = prod.product.values
        if not any(np.isin(existing_idx, add_idx)):
        #if the products does not already exist, then append
            prod = prod.reindex(product = list(existing_idx) + add_idx)    
            height_levels = height_levels.reindex(product = list(existing_idx) + add_idx)    
            
    if len(prod) == 0 and len(signal) > 0: 
        prod = \
            xr.DataArray(coords=[#signal.time.values, 
                                 prod_id,#prod_map.prod_id.values, 
                                 signal.bins.values],
                         dims=['product', 'bins'],  #'time',
                         name='products')
                                         #channel   #altitude
        height_levels = \
            xr.DataArray(coords=[#signal.time.values, 
                                 prod_id,#prod_map.prod_id.values, 
                                 signal.bins.values],
                         dims=['product', 'bins'],  #'time',
                         name='products')
                                         
    if len(info_prod) == 0:
        info_prod = pd.DataFrame(index = prod.product.values, 
                                 columns=['wave', 'type', 'type_short',
                                          'units', 'method', 'polarization','process'])
    
    return(prod, height_levels, info_prod)
