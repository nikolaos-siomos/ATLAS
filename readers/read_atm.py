"""
@author: N.Siomos and P. Paschou

Functions for reading the temperature and pressure profiles from
- radiosonde
- NWP model (WRF)
"""
import numpy as np
import glob
import os
import datetime as dt
from helper_functions.helper_functions import is_var_number

def rsonde(dir_path):  
    
    alt=[]
    prs=[]
    tem=[]
    hum=[]
    imol = 3

    try:
		#fname = glob.glob(os.path.join(dir_path, 'SO*'))
        fname = glob.glob(os.path.join(dir_path, '*.txt'))
        
        if len(fname) > 1:
            print('-- Warning! Too many radiosonde files found in folder --> \n'
				  'Using U.S. standard atmosphere instead..')
            imol=1 # change the imol indicator for US standard
            
        if len(fname) == 0:
            print('---- Warning! No radiosonde file found in folder --> \n'
				  'Using U.S. standard atmosphere instead..')

            imol=1 # change the imol indicator for US standard

        if len(fname) == 1:
            data = np.loadtxt(fname[0],skiprows = 7, delimiter='~', dtype=str)			

            prs = np.nan*np.zeros(len(data) + 1)
            alt = np.nan*np.zeros(len(data) + 1)
            tem = np.nan*np.zeros(len(data) + 1)
            hum = np.nan*np.zeros(len(data) + 1)
			
            for j in range(1,len(data) + 1):
                jj = j-1
                prs[j] = float(data[jj][0:8])
                alt[j] = float(data[jj][8:15])
                tem[j] = float(data[jj][15:20]) + 273.15
                hum[j] = float(data[jj][28:35])
			
        # assign params also for zero alt        
        prs[0] = prs[1]
        alt[0] = 0.
        tem[0] = tem[1]    
        hum[0] = hum[1]    
        
        mask_empty = (prs == prs) & (alt == alt) & (tem == tem) & (hum == hum)
        
        prs = prs[mask_empty]
        alt = alt[mask_empty]
        tem = tem[mask_empty]
        hum = hum[mask_empty]
    except:
        raise OSError('Please check the /Atmosphere/rsonde/ directory! Something is missing or doesn\'t exist...')

    return(alt,prs,tem,hum,imol)

def model(dir_path):
    
    alt=[]
    prs=[]
    tem=[]
    hum=[]
    imol = 2
    
    try:
        # fname = glob.glob(os.path.join(dir_path, '*.txt'))
        fname = os.listdir(dir_path)
        
        if len(fname) > 1:
            print('---- Warning! Too many model files found in folder --> \n'
                  'Using U.S. standard atmosphere instead..')
            imol=1 # change the imol indicator for US standard
        
        if len(fname) == 0:
            print('---- Warning! No model file found in folder --> \n'
                  'Using U.S. standard atmosphere instead..')
            imol=1 # change the imol indicator for US standard
            
        if len(fname) == 1:
            data = np.loadtxt(os.path.join(dir_path,fname[0]), skiprows = 0, 
                              delimiter='~', dtype=str)
		
            # checks if the first row is a header. If it is, then removes first row
            data = check_header(data)
            
            # create arrays for P, T, RH, alt 
            prs = np.nan*np.zeros(len(data) + 1)
            alt = np.nan*np.zeros(len(data) + 1)
            tem = np.nan*np.zeros(len(data) + 1)
            hum = np.nan*np.zeros(len(data) + 1)
            
            col = len(data[0].split(',')) # columns size
            
            for j in range(1,len(data)+1):
                split_data = data[j-1].split(',')
                
                # The columns are: [datetime; optional], P, T, Dew Point, RH, alt
                prs[j] = float(split_data[col-5].strip()) # pick 1st or 2nd element for P
                alt[j] = float(split_data[col-1].strip()) # pick last element for alt
                tem[j] = float(split_data[col-4].strip()) + 273.15 # pick 1st or 2nd col for T
                hum[j] = float(split_data[col-2].strip()) # pick one before last element for RH
            
            # assign params also for zero alt
            prs[0] = prs[1]
            alt[0] = 0.
            tem[0] = tem[1]    
            hum[0] = hum[1]    
            
            # mask to remove empty bins
            mask_empty = (prs == prs) & (alt == alt) & (tem == tem) & (hum == hum)
    
            prs = prs[mask_empty]
            alt = alt[mask_empty]
            tem = tem[mask_empty]
            hum = hum[mask_empty]
    except:
        raise OSError('Please check the /Atmosphere/model/ directory! Something is missing or doesn\'t exist...')
        
    return(alt,prs,tem,hum,imol)

def extract_rsonde_info(dir_path):
    fname = glob.glob(os.path.join(dir_path, '*.txt'))

    if len(fname) == 1:
        data = np.loadtxt(fname[0], max_rows=1, dtype=object)
        st_id = data[0]
        loc = data[2]
        t = f'{data[-4][:2]}UTC'
        m = str(dt.datetime.strptime(data[-2],'%b').month)
        date = f'{data[-3]}.{m}.{data[-1][-2:]}'
        info = ', '.join([loc,st_id,date,t])
    else:
        info = 'location, station_ID, dd.mm.yy, hhUTC'
        print('---- Warning! No or Too many radiosonde file(s) found in folder!\n'
              f'Returns no info about {info}')
        
    return(info)

def check_header(arr):
    # check if the first row is a header. If it is then remove first row
    # splits the 1st row in its elements
    row1 = arr[0].split(',')
    
    isnum = is_var_number(row1[-1]) # check the last element
    
    if not isnum:
        
        arr = np.delete(arr,0)

    return(arr)