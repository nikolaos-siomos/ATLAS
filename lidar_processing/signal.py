"""
@authors: N. Siomos & P. Paschou

Processing routines for signals

=================================
Signal in 3D xarray dataset with dimensions [time, channel, bin/range]

fucntions
 --f1:  unit conversion
 --f2:  average by time
 --f3:  resolution change
 --f4:  signal dimensions conversion
 --f5:  dead time correction
 --f6:  trigger delay correction
 --f7:  background calculation 
 --f8:  background subtraction
 --f9:  range correction
 --f10: dark signal structure correction
 --f11: signal upper trimming
 --f12: construction of total signal
 --f13: total signal depolarization correction 
"""
import numpy as np
from molecular import height_scales
from lidar_processing.construct import add_param_err
from helper_functions.helper_functions import datetime_unit, check_timedelta
import xarray as xr

#--f1 a)
def unit_conv_bits_to_mV(signal, info):
    if len(signal)>0:
        for j in info.index:
            ch = dict(channel = j)
            if info.ch_mode[j] == 0.0: 
                # analog conversion (to mV)
                signal.loc[ch] = signal.loc[ch].values*info.ADC_range[j]/(info.shots[j]*(np.power(2,info.ADC_bit[j]) -1))
        
    return(signal) 
#--f1 b)
def unit_conv_mV_to_summed_mV(signal, info):

    if len(signal)>0:
        for ind in info.index:
                
            shots = info.shots.loc[ind]
            
            ch = dict(channel = ind)
                
            signal.loc[ch] = signal.copy().loc[ch] * shots
    
    return(signal) 
#--f1 c)
def unit_conv(signal, info):

    if len(signal)>0:
        for j in info.index:
            ch = dict(channel = j)
            if info.ch_mode[j] == 0.0: 
                # analog conversion (to mV)
                signal.loc[ch] = signal.loc[ch].values*info.ADC_range[j]/(info.shots[j]*(np.power(2,info.ADC_bit[j]) -1))
            else:
                #photon counting conversion (to MHz)
                signal.loc[ch] = signal.loc[ch].values * info.sampl_rate[j]/info.shots[j]
    return(signal) 
#--f1 d)
def sqd_to_std(signal_sqd, info):
    
    signal_std = np.nan*signal_sqd.copy()
   
    for j in info.index:
        ch = dict(channel = j)
        sqr = signal_sqd.loc[ch].values
        signal_std.loc[ch] = info.shots[j]*sqr/(np.sqrt(info.shots[j]*(info.shots[j]-1)))
        
    return(signal_std)
#--f1 e)
def unit_conv_counts_to_MHz(signal, info):
    if len(signal)>0:    
        for j in info.index:
            if info.ch_mode[j] == 1.0: 
                ch = dict(channel = j)
                #photon counting conversion (Counts to MHz)
                signal.loc[ch] = signal.loc[ch].values * \
                    info.sampl_rate[j]/info.shots[j]
    return(signal) 
#--f1 f)
def unit_conv_MHz_to_counts(signal, info):
    if len(signal)>0:
        for j in info.index:
            if info.ch_mode[j] == 1.0:
                ch = dict(channel = j)
                #photon counting conversion (MHz to counts)
                signal.loc[ch] = signal.loc[ch].values * \
                    info.shots[j]/info.sampl_rate[j]
    return(signal)       

#--f2.1
def average_by_time(signal, timescale):

    adj_timescale = timescale

    time_arr = signal.time.values
    # average the signal only if the timescale is not empty and the signal time dim is >1
    if timescale != '' and len(time_arr) > 1:
                
        # calculate the time duration (min) between first and last measurement
        dt = time_arr[-1] - time_arr[0]
        meas_dur = np.ceil(dt / np.timedelta64(1, 'm'))
            
        # convert timescale in minutes
        t_avg = int(''.join(filter(str.isdigit, timescale)))
        freq = (datetime_unit(timescale)).lower()
        
        if freq == 'h':
            t_avg = t_avg * 60.
        if freq == 'd':
            t_avg = t_avg * 24. * 60.
        if freq == 's':
            t_avg = t_avg / 60.
        
        freq = 'm'
        
        #if the measurement time < given timescale, adjust the returned timescale 
        adj_timescale = check_timedelta(signal.time.values, timescale, freq)        
        
        # define how many averaged timeframes will be created
        avg_num = int(np.ceil(meas_dur/t_avg))
        
        # Calculate the new time index based on the timescale of averaging.
        # -> It is also the starting time in the time slice for averaging
        t_index = [(time_arr[0]+ii*np.timedelta64(t_avg,freq)) for ii in range(avg_num)]
        # print(f't_index: {t_index}')
        
        # Calculate the ending time in the time slice for averaging
        t_slice = [(time_arr[0]+(ii+1)*np.timedelta64(t_avg,freq)) for ii in range(avg_num)]
        
        # Create a new xarray with reduced size of time dimension based on number of averaged timeframes
        sig_avg = signal.reindex(time=t_index, fill_value=np.nan)
        
        count_occur = []
        
        for i in range(avg_num):
            # select a slice of timeframes to be averaged in one timeframe
            sel_tslice = dict(time=slice(t_index[i],t_slice[i]))
            t_idx = dict(time = i)
            
            # counting the occurences in each slice of timeframe
            count_occur.append(signal.loc[sel_tslice].time.size)
            
            # Insert in the new xarray the averaged signal from the sliced timeframes
            sig_avg[t_idx] = signal.loc[sel_tslice].mean(dim='time').values #3-d array [iters, channel, bins]
            
        frames = xr.DataArray(data = count_occur, 
                              dims = ['time'], 
                              coords = [sig_avg.time.values])
        
        # Mask the averaged timeframes for time occureneces below 95% of the max count occurence
        mask_time = (count_occur < np.max(count_occur)*0.95)

        #Drop the timeframes with occurences below the 95% of the max count occurence
        signal = sig_avg.drop_sel(time = sig_avg.time.values[mask_time])
        frames = frames.drop_sel(time = frames.time.values[mask_time])

    else:
    # If no timescale ('') is given then returns the raw temporal resolution
       
        # calc the timedelta between timeframes to find raw temporal resolution   
        deltat = [(time_arr[t+1]-time_arr[t])/np.timedelta64(1,'s') for t in range(time_arr.size-1)]

        if len(deltat) > 0:
            # minimum time delta -> raw temporal resolution (+- 1 sec)
            t_resol = int(min(deltat))
            adj_timescale = f'{t_resol}s'     
         
        frames = xr.DataArray(data = np.ones(time_arr.size), 
                  dims = ['time'], 
                  coords = [time_arr])
    # elif timescale == '' and len(time_arr) > 1:
    # # If no timescale ('') is given then returns the raw temporal resolution
       
    #     # calc the timedelta between timeframes to find raw temporal resolution   
    #     deltat = [(time_arr[t+1]-time_arr[t])/np.timedelta64(1,'s') for t in range(time_arr.size-1)]
    #     if len(deltat) > 0:
    #         # minimum time delta -> raw temporal resolution (+- 1 sec)
    #         t_resol = int(min(deltat))
    #         adj_timescale = f'{t_resol}s'     
         
    #         frames = xr.DataArray(data = np.ones(time_arr.size), 
    #                   dims = ['time'], coords = [time_arr])
    # else: 
    #     frames = xr.DataArray(data = np.ones(time_arr.size), 
    #                           dims = ['time'], coords = [time_arr])      

    return(signal, frames, adj_timescale)

#--f2.2 (implementing time resampling)    
def average_by_time_resampling(signal, timescale):

    adj_timescale = timescale

    # average the signal if the timescale is not empty
    if timescale != '':
        # Sampling the measurements considering the time of the first measurement 
        indicator = datetime_unit(timescale)
    
        # #Set the base_val to start sampling with respect to the first measurement            
        # base_val = int(signal.time[0].dt.strftime(f'%{indicator}').values)

        count_occur = signal.time.resample(time=timescale, skipna=True, #base=base_val,
                                            closed='left', label='left').count(dim='time').values

        #if the measurement time is less than the given timescale, adjust the timescale
        adj_timescale = check_timedelta(signal.time.values, timescale, indicator)        

        # Mask the resampled timeframes for time occureneces below 75% of the max count occurence
        mask_time = (count_occur < np.max(count_occur)*0.75)
        
        # Time resampling and averaging of the signals inside each timeframe
        signal = signal.resample(time=timescale, skipna=True, #base=base_val,
                                 closed='left', label='left').mean(dim='time', skipna=True)
        
        #Drop the timeframes with occurences below the 75% of the max count occurence
        signal = signal.drop_sel(time = signal.time.values[mask_time])                                                                   
        
    return(signal, adj_timescale)

#--f3
def change_resol(signal, info, resol):  

    bins = signal.bins.values
    sig = signal.values
    step = info.resol.values
    
    iters = signal.iters.values
    time = signal.time.values
    channel = signal.channel.values
        
    if resol == '':
        resol = np.max(info.resol.values)  
    
    alt_l = bins[-1]*np.max(step)

    # x_com = np.arange(0., alt_l + resol, resol)
    x_com = np.arange(resol, alt_l + resol, resol)
    y_com = np.nan*np.zeros([sig.shape[0], sig.shape[1], 
                             sig.shape[2], x_com.size])
    
    # interpolate the signal for common spatial resolution in all channels.
    # for common raw resolution in all channels the signal remains intacted
    for j in range(sig.shape[2]):
        x = bins*step[j]
        for k in range(sig.shape[1]):
            for n in range(sig.shape[0]):
                y = signal[n, k, j, :]
                y_com[n, k, j, :] = np.interp(x_com, x, y)
    
    
    #convert the signal xarray dim 'bins' to 'range' => dims [iters, time, channel, range]
    signal = xr.DataArray(y_com, 
                          coords=[iters, time, channel, x_com],
                          dims=['iters', 'time', 'channel', 'range'])

    info.resol[:] = resol

    return(signal, info, resol)    

#--f4
def change_dim(signal, dim_name, new_dim_name, new_array):
    if len(signal)>0:
        signal[dim_name] = new_array
        rename_dim = {dim_name:new_dim_name}
        signal = signal.rename(rename_dim)
        #signal[new_dim_name] = signal[new_dim_name].reset_index(new_dim_name)
    return(signal)

#--f5
def dead_time_correction_MHz(signal, info, dead_time, ipmt):
    
    for j in info.index:
        if info.ch_mode.loc[j] == 1.: # only for pc channels

            if ipmt == 1: # Paralyzable
                print('Warning! Dead time correction for Paralyzable system not implemented yet!')
                print('--> Applyinf Dead time correction for non- Paralyzable system intead...')
                ipmt = 0
                #temp_sig = pre_pro.correct_dead_time_paralyzable(temp_sig, time_interval[j], dead_time)

            #signal.values[signal.values > 1000./dead_time] = 0.
            #TODO: run a test to check for saturated signals (like in Giannis function for deadtime)    

            if ipmt == 0: # Non-paralyzable
                ch = dict(channel = j)
                signal.loc[ch] = signal.loc[ch] / (1. - signal.loc[ch] * dead_time.loc[j] * 1e-3)

            
    return(signal)

#--f5
def dead_time_correction_counts(signal, info, dead_time, ipmt):
    
    for j in info.index:
        if info.ch_mode.loc[j] == 1.: # only for pc channels
            if ipmt == 0: # Non-paralyzable
                ch = dict(channel = j)
                T = 1E3 / info.sampl_rate[j]  # in nanoseconds (like the dead time)
                shots = info.shots[j]
                dt = dead_time.loc[j]
                # The counts in the denominator must be per shot
                # The counts in the nomitor are the summed ones (no need to devide with the shots)
                signal.loc[ch] = signal.loc[ch] / \
                    (1. - signal.loc[ch] * dt / (T * shots))
            if ipmt == 1: # Paralyzable
                print('Warning! Dead time correction for Paralyzable system not implemented yet!')
                #temp_sig = pre_pro.correct_dead_time_paralyzable(temp_sig, time_interval[j], dead_time)
            
    return(signal)

#--f6
def trigger_delay(signal, info, trigger_delay_bins):
    '''
    Corrects the signal in case of pre-triggering or trigger delays
    
    signal.coords = [time, channel, range]
    trigger_dealy = array of the trigger delay bins or sec (can defined in itr) per channel
    '''

    
    for t in range (signal.time.size):
        ind_t = dict(time = t)
        
        for j in signal.channel.values:
           ind_ch = dict(channel = j)
           signal[ind_t].loc[ind_ch] = correct_trigger_delay_bins(signal[ind_t].loc[ind_ch].values, 
                                                             trigger_delay_bins.loc[j])
    return(signal)

#--f7
def background(signal, llim_arr, ulim_arr):
    '''
    Calculation of the background value in the pre-trigger region or in the end of signal region
    signal = 3D xarray [time, channel, bin/range]
    llim = lower limit (bin) for the background calculation
    ulim = upper limit (bin) for the background calculation
    '''
    bc = np.nan*signal.mean(dim = 'bins')
    
    for ch in signal.channel.values:
        ch_ind = dict(channel = ch)
        llim = llim_arr.loc[ch]
        ulim = ulim_arr.loc[ch]
        bc.loc[ch_ind] = signal.loc[ch_ind].isel(bins = slice(llim,ulim)).mean(dim = 'bins', skipna = True)
        
    return(bc)

#--f8        
def background_correction(signal, bc):
    signal = signal - bc
    return(signal)

#--f9
def range_correction(signal):
    signal = signal*np.power(signal.range.values, 2)
    return(signal)

#--f10
def dark_correction(signal, info, dark): 
    
    an_ch = dict(channel = info.index[info.ch_mode == 0].values) #analog channel indexes
    
    drk_mean = dark.mean(dim='time', skipna = True).loc[an_ch]
    
    for t in range (signal.time.size):
        ind_t = dict(time = t)      
        signal[ind_t].loc[an_ch] = signal[ind_t].loc[an_ch] - drk_mean

    return(signal)
  
#--f11
def trim(signal, ground, SZA, alt_lim):
        
    alt = height_scales.range_to_alt(signal.range.values, ground = ground, SZA = SZA)
    
    if alt_lim > alt[-1]:
        alt_lim = alt[-1]
        print(('---- Warning! Given trim altitude too high --> \r\n '+\
              'Selecting the end of signal as the uppermost limit'))
    
    alt = alt[alt <= alt_lim]
    x_range = height_scales.alt_to_range(alt, ground = ground, SZA = SZA)
    
    signal = signal.sel(range=x_range)
  
    return(signal, alt_lim)

#--f12
def construct_total(signal, info, tot):
    '''

    Parameters
    ----------
    signal : 3-D xarray dataset
        Contains the signals with dimensions [time, channel, range].
    info : object type dataframe
        Contains info about each channel.
    itot : integer
        Indicates whether to construct the total signal.
    tot : class with dataframes as attributes
        Contains the info (map and parameters) from the total_map.ini file.


    Returns
    -------
    None.

    * Itot ~ cal_f*H_r*I_t - H_t*I_r
    
    The formula is based on eq. 65 from Freudenthaler et al. 2016
    We assume the caf_f is for the ratio R/T regardless which is the parallel/cross channel
    '''
    
    # extract the info from the total_map.ini file    
    id_tot = tot.map.id_new
    id_r   = tot.map.id_r
    id_t   = tot.map.id_t

    cal_f  = tot.params.cal_f
    cal_f_err = tot.params.cal_f_err
    H_r    = tot.params.h_r
    H_r_err    = tot.params.h_r_err
    H_t    = tot.params.h_t
    H_t_err    = tot.params.h_t_err
        
    for j in id_tot.index:
        #Add new index in the dimension of channel for the total signals
        signal = signal.reindex(channel = list(signal.channel.values)+[id_tot[j]])    
        
        #Add the info about the total signal in the info_array
        info.loc[id_tot[j]] = info.loc[id_r[j]].copy()
        info.loc[id_tot[j], 'full_ovl_idx'] = \
            max([info.loc[id_t[j], 'full_ovl_idx'],info.loc[id_r[j], 'full_ovl_idx']])
        info.loc[id_tot[j], 'ch_pol'] = 'u' # for calibrated sum instead of 'o'

        print(f'-> Calibrated sum: Adding channel {id_tot[j]}')


        # Add statistical error to params
        c_f = add_param_err(cal_f[j], cal_f_err[j], signal.iters.size) # 1D array with size = iters
        H_R = add_param_err(H_r[j], H_r_err[j], signal.iters.size)
        H_T = add_param_err(H_t[j], H_t_err[j], signal.iters.size)
        
        #Construct the total signal for all time dims        
        for n in range(signal.iters.size):
            ind = dict(iters = n)
            
            for t in range(signal.time.size):                
                ind['time'] = t
                
                temp_sig = c_f[n] * H_R[n] * signal[ind].loc[dict(channel = id_t[j])].values -\
                    (H_T[n]*signal[ind].loc[dict(channel = id_r[j])].values)
                
                signal[ind].loc[dict(channel = id_tot[j])] = temp_sig


    return(signal, info)

#--f13
def correct_total_sig(signal, info, prod, info_prod, vdr_mp):
    '''
    Correction of the detected total signal using the retrieved volume depol ratio profile.
    
    The total signal is detected in the Reflected path
    
    References
    -----    
    Engelmann et al., 2016, Atmos. Meas. Tech., https://doi.org/10.5194/amt-9-1767-2016
    Mattis et al., 2009, Appl. Opt, https://doi.org/10.1364/AO.48.002742
    
    '''
    icor=False
    if len(prod)>0 and len(vdr_mp)>0 and len(signal)>0:

        n_size = signal.iters.size
    
        for j in vdr_mp.index:

            id_vldr = vdr_mp.prod_id.values[j]
            id_tot = vdr_mp.ir_id.values[j]
            
            # if vldr exists, have the same wl with signal amd signal is total then procced to the correction
            if (id_vldr in prod.product.values) and (info.ch_pol.loc[id_tot] =='o') \
                and (float(info.wave.loc[id_tot])==float(info_prod.wave.loc[id_vldr])):
                
                icor=True
                
                tr_tot = float(vdr_mp.tr_r.values[j])
                tr_tot_err = float(vdr_mp.tr_r_err.values[j])
                TR_tot = add_param_err(tr_tot, tr_tot_err, n_size)
        
                v_ind = dict(product=id_vldr)
                s_ind = dict(channel=id_tot)
                
                for n in range(signal.iters.size):
                    v_ind['iters' ]= n
                    s_ind['iters' ]= n
                    vdr_temp = prod.loc[v_ind].values
                    sig_temp = signal.loc[s_ind].copy().values
                
                    cor_factor = (TR_tot[n] * vdr_temp + 1) / (vdr_temp + 1) 
                    signal.loc[s_ind] = sig_temp * cor_factor
            
    return(signal, icor)

#-- additional functions

def correct_trigger_delay_bins(signal, trigger_delay_bins):
    """
    Shifts the signal for one channel by a specified number of bins.

    Parameters
    ----------
    signal: float array
        The lidar signal profile. Can be either 1D or 2D with dimensions (time, range).
    trigger_delay_bins: int
        Number of bins to shift the signal

    Returns
    -------
    signal: float array
        Corrected lidar signal profile (shifted)
    
    Notes
    -----
    The number of bins to be shifted should be positive. For negative values
    the shift will be done in reverse.

    The function is obtained from:
    https://gitlab.com/ioannis_binietoglou/lidar-processing/
    
    #################################################################################
    # The MIT License (MIT)
    # 
    # Copyright (c) 2015, Ioannis Binietoglou
    # 
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    # 
    # The above copyright notice and this permission notice shall be included in all
    # copies or substantial portions of the Software.
    # 
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    # SOFTWARE.
    #################################################################################
    
    """
    # Force 2D data
    if signal.ndim == 1:
        signal = signal[np.newaxis, :]
        dimension_added = True
    else:
        dimension_added = False

    signal = np.roll(signal, trigger_delay_bins, axis=1)

    # Remove rolled bins
    if trigger_delay_bins >= 0:
        signal[:, :trigger_delay_bins] = np.nan
    else:
        signal[:, trigger_delay_bins:] = np.nan

    #signal = np.ma.masked_invalid(signal)

    # Reduce back to 1D if required
    if dimension_added:
        signal = signal[0, :]

    return (signal)