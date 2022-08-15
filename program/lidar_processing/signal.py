"""
@authors: N. Siomos, P. Paschou, Ioannis Binietoglou 
based on SULA project (https://react-gitlab.space.noa.gr/ReACT/eve/data-processing)
and also on https://gitlab.com/ioannis_binietoglou/lidar-processing/

Processing routines for signals 

=================================
Signal in 3D xarray dataset with dimensions [time, channel, bin/range]

fucntions
 -- background calculation: Calculated the solar background per timeframe and channel 
 -- dead time correction: Applied to the photon channels to account for the dead time effect
 -- trigger_delay: Perform the trigger correction per channel
 -- trim_vertically: Trim channels up to a maximum altitude
 -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels
 --f2:  average by time
 --f3:  resolution change
 --f4:  signal dimensions conversion
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
# from helper_functions.helper_functions import datetime_unit, check_timedelta
import xarray as xr
import sys

def average_by_time(sig, timescale):

    """
    General:
        Averages the lidar signals across the time deminsion 
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the following
            dimensions: (time, channel, bins). The photon signals must not be
            dead time corrected yet
            
        timescale: 
            The temporal avraging step, in hours. E.g. if set to 1 then 1h
            averages will be generated from the original signals

              
    Returns:
        
        sig_out: 
            An xarray in the same dimensions with sig but with the time dim
            reduced by averaging, according to the provided timescale
        
        frames:
            An xarray with a single dimension ('time') which coincides with the
            sig_out time dimension. It provied the number of averaged profiles
            per temporal timeframe of sig_out
            
    """
    
    time_arr = sig.time.values

    # Calculate the time duration (hours) between first and last measurement
    dt = time_arr[-1] - time_arr[0]
    
    meas_dur = np.ceil(dt / np.timedelta64(60,'m'))
    
    # average the sig only if the timescale is not empty and the sig time dim is >1
    if np.isnan(timescale) or \
        (not np.isnan(timescale) and timescale >= meas_dur):
                
        m_time = time_arr[0] + dt / 2.
        
        sig_avg = sig.mean(dim = 'time').values
        
        sig_out = sig.reindex(time = [m_time], fill_value = np.nan)

        sig_out.loc[dict(time = m_time)] = sig_avg
        
        # Store the time frame information per slice to use it for the error simulation
        frames = xr.DataArray(data = time_arr.size, 
                              dims = ['time'], 
                              coords = [sig_out.time.values])
        

    else:
        halfscale_m = int(timescale * 60/ 2)

        timescale_m = int(timescale * 60)

        # Define how many averaged timeframes will be created
        avg_num = int(np.ceil(meas_dur/timescale))
        
        # Calculate the new time index based on the timescale of averaging.
        # -> It is also the middle time in the time slice for averaging
        m_time = [(time_arr[0] + ii * np.timedelta64(halfscale_m,'m')) 
                  for ii in range(avg_num)]

        # Calculate the starting time in the time slice for averaging
        s_time = [(time_arr[0] + ii * np.timedelta64(timescale_m,'m')) 
                  for ii in range(avg_num)]

        # Calculate the ending time in the time slice for averaging
        e_time = [(time_arr[0] + (ii + 1) * np.timedelta64(timescale_m,'m')) 
                  for ii in range(avg_num)]        

        # Create a new xarray with reduced size of time dimension based on number of averaged timeframes
        sig_out = sig.reindex(time = m_time, fill_value = np.nan)
        
        count_occur = np.zeros(avg_num)
        
        for i in range(avg_num):
            # Select a slice of timeframes to be averaged in one timeframe
            tslice = dict(time = slice(s_time[i], e_time[i]))
            
            t_idx = dict(time = i)
            
            # Counting the occurences in each slice of timeframe
            count_occur[i] = sig.loc[tslice].time.size
            
            # Insert in the new xarray the averaged sig from the sliced timeframes
            sig_out[t_idx] = sig.loc[tslice].mean(dim='time').values

        # Store the time frame information per slice to use it for the error simulation
        frames = xr.DataArray(data = count_occur, 
                              dims = ['time'], 
                              coords = [sig_out.time.values])
        
        # Mask the averaged timeframes for time occureneces below 75% of the max count occurence
        mask_time = (count_occur < np.max(count_occur) * 0.75)

        #Drop the timeframes with occurences below the 75% of the max count occurence
        sig_out = sig_out.drop_sel(time = sig_out.time.values[mask_time])

        #Drop the timeframes with occurences below the 75% of the max count occurence        
        frames = frames.drop_sel(time = frames.time.values[mask_time])

    return(sig_out, frames)

def background(sig, lower_bin, upper_bin):

    """
    General:
        Calculates the solar background based on the provided metadata
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the following
            dimensions: (time, channel, bins). The photon signals must not be
            dead time corrected yet
            
        lower_bin: 
            A pandas series with the lower bin applied for the background
            correction. It can be in the pretrigger or at the end of the
            signal profile. Bin numbering starts at 1            

        upper_bin: 
            A pandas series with the upper bin applied for the background
            correction. It can be in the pretrigger or at the end of the
            signal profile. Bin numbering starts at 1              
    Returns:
        
        bgr: 
            An xarray with the following dimensions (time, channel) similar to
            sig. It contains the solar background signal
            
    """
    channels = sig.channel.values

    bgr = np.nan*sig.mean(dim = 'bins').copy()
    
    for ch in channels:
        
        ch_d = dict(channel = ch)
        
        llim = lower_bin.loc[ch]
        
        ulim = upper_bin.loc[ch]
        
        bin_d = dict(bins = slice(llim,ulim))
        
        bgr.loc[ch_d] = sig.loc[ch_d].loc[bin_d].mean(dim = 'bins', skipna = True)
        
    return(bgr)

def dead_time_correction(sig, dead_time, dead_time_cor_type):
   
    """
    General:
        Applies the dead time correction for all photon channels
        
    Input:
        sig: 
            A 2D or 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, ...). 
            
        dead_time: 
            A pandas series with the dead time per channel in nanoseconds. 
            The index should correspond to the channel dimension of sig            

        dead_time_cor_type: 
            A pandas series with the dead time correction type per channel 
            (0 for non paralyzable or 1 for non paralyzable)
            The index should correspond to the channel dimension of sig            
            
    Returns:
        
        sig_out: 
            An xarray in the same shape as sig with the dead time corrected
            signals
            
    """
    
    # Check if the Paralyzable Deadtime Correction Type is 1 (Paralyzable)--> if not provide warning and switch to Non Paralyzable (0)
    if (dead_time_cor_type.values == 1).any():
        print('-- Warning: A Dead time correction for Paralyzable system not implemented yet!')
        print('--> Applying Dead time correction for non- Paralyzable system intead...')
               
    channels = sig.channel.values
    
    sig_out = np.nan * sig.copy()
    
    for ch in channels:

        ch_d = dict(channel = ch)

        if ch[2] == 'p': #3rd digit of channel name is the acquisition mode (a or p)

            if (sig.loc[ch_d].values > 1000./dead_time[ch]).any():
                print(f"-- Warning: Channel {ch} - A countrate value above the maximum allowed value was detected! Consider revising the input file ")

            if dead_time_cor_type[ch] == 0: # Non-paralyzable
                sig_out.loc[ch_d] = sig.loc[ch_d] / \
                    (1. - sig.loc[ch_d] * dead_time.loc[ch] * 1e-3)

        if ch[2] == 'a': #3rd digit of channel name is the acquisition mode (a or p)
            
            sig_out.loc[ch_d] = sig.loc[ch_d]
            
    return(sig_out)

def trigger_delay(sig, trigger_delay_bins):

    """
    General:
        Shift the bin axis to the left (pre-trigger) or to the 
        right (trigger delay) in order to bring all channels to the same 
        ground bin. Empty bins are filled with nans
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, bins). 
            
        trigger_delay_bins: 
            A pandas series with the trigger delay bins per channel. 
            Negative values correspond to pretriggering --> shift to the left
            The index should correspond to the channel dimension of sig            

    Returns:
        
        sig_out: 
            An xarray in the same shape as sig with the trigger corrected
            signals
            
    """

    channels = sig.channel.values

    sig_out = np.nan * sig.copy()
        
    for ch in channels:
        
        ch_d = dict(channel = ch)
        
        shift_bins = trigger_delay_bins.loc[ch]
        
        sig_out.loc[ch_d] = sig.loc[ch_d].shift(bins = shift_bins).values
        
    return(sig_out)

def trim_vertically(sig, ground_alt, zenith_angle, alt_lim, resol):
 
    """
    General:
        Removes all the bins above the maximum allowed altitude. If the signals
        have different vertical resolution the highest resolution is used that 
        correpsonds to the maximum allowed bins. In order for this correction
        to make sense the signals must be corrected for the trigger delay
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, bins). 
            
        ground_alt: 
            A float with the ground altitude of the lidar in m
            
        zenith_angle:
            A float with the zenith angle where the lidar points to
            
        alt_lim:
            A float with the maximum altitude above which no computations will
            be performed
            
        resol:
            A pandas series with the range resolution in m per channel. 
            The index should correspond to the channel dimension of sig
            
    Returns:
        
        sig_out: 
            An xarray in the same shape as sig with the signals trimmed above
            the maximum allowed altitude
            
    """
       
    sig_out = np.nan * sig.copy()

    rng_lim = (alt_lim - ground_alt)/np.cos(zenith_angle)
    
    bin_lim = np.int(rng_lim / np.min(resol))
    
    bins = sig.bins.values
    
    bin_d = dict(bins = slice(bins[0], bin_lim))
    
    if bin_lim < np.max(bins):
        
        sig_out = sig.loc[bin_d]
  
    return(sig_out)


def unit_conv_counts_to_MHz(sig, shots, resol):
    
    """
    General:
        Converts photon units from summed counts (licel default) to MHz
        
    Input:
        sig: 
            A 2D or 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, ...). The units of the photon 
            channels must be raw counts
            
        shots: 
            A 2D xarray with the laser shots per channel and timeframe.
            It should include the following dimensions: (time, channel, ...) 
            The index should correspond to the channel dimension of sig            

        dead_time_cor_type: 
            A pandas series with the dead time correction type per channel 
            (0 for non paralyzable or 1 for non paralyzable)
            The index should correspond to the channel dimension of sig            
            
    Returns:
        
        sig_out: 
            An xarray in the same shape as sig. The units of the photon 
            channels are converted from raw counts to MHz
            
    """
    
    channels = sig.channel.values

    sig_out = np.nan * sig.copy()
    
    for ch in channels:
        
        ch_d = dict(channel = ch)
    
        if ch[2] == 'p': #3rd digit of channel name is the acquisition mode (a or p)
            
            sampl_rate = 150. / resol.loc[ch]
            
            sig_out.loc[ch_d] = sig.loc[ch_d] * sampl_rate / shots.loc[:,ch]
        
        if ch[2] == 'a': #3rd digit of channel name is the acquisition mode (a or p)
            
            sig_out.loc[ch_d] = sig.loc[ch_d]
    
    return(sig_out) 



# def average_by_time(signal, timescale):

#     adj_timescale = timescale

#     time_arr = signal.time.values
#     # average the signal only if the timescale is not empty and the signal time dim is >1
#     if timescale != '' and len(time_arr) > 1:
                
#         # calculate the time duration (min) between first and last measurement
#         dt = time_arr[-1] - time_arr[0]
#         meas_dur = np.ceil(dt / np.timedelta64(1, 'm'))
            
#         # convert timescale in minutes
#         t_avg = int(''.join(filter(str.isdigit, timescale)))
#         freq = (datetime_unit(timescale)).lower()
        
#         if freq == 'h':
#             t_avg = t_avg * 60.
#         if freq == 'd':
#             t_avg = t_avg * 24. * 60.
#         if freq == 's':
#             t_avg = t_avg / 60.
        
#         freq = 'm'
        
#         #if the measurement time < given timescale, adjust the returned timescale 
#         adj_timescale = check_timedelta(signal.time.values, timescale, freq)        
        
#         # define how many averaged timeframes will be created
#         avg_num = int(np.ceil(meas_dur/t_avg))
        
#         # Calculate the new time index based on the timescale of averaging.
#         # -> It is also the starting time in the time slice for averaging
#         t_index = [(time_arr[0]+ii*np.timedelta64(t_avg,freq)) for ii in range(avg_num)]
#         # print(f't_index: {t_index}')
        
#         # Calculate the ending time in the time slice for averaging
#         t_slice = [(time_arr[0]+(ii+1)*np.timedelta64(t_avg,freq)) for ii in range(avg_num)]
        
#         # Create a new xarray with reduced size of time dimension based on number of averaged timeframes
#         sig_avg = signal.reindex(time=t_index, fill_value=np.nan)
        
#         count_occur = []
        
#         for i in range(avg_num):
#             # select a slice of timeframes to be averaged in one timeframe
#             sel_tslice = dict(time=slice(t_index[i],t_slice[i]))
#             t_idx = dict(time = i)
            
#             # counting the occurences in each slice of timeframe
#             count_occur.append(signal.loc[sel_tslice].time.size)
            
#             # Insert in the new xarray the averaged signal from the sliced timeframes
#             sig_avg[t_idx] = signal.loc[sel_tslice].mean(dim='time').values #3-d array [iters, channel, bins]
            
#         frames = xr.DataArray(data = count_occur, 
#                               dims = ['time'], 
#                               coords = [sig_avg.time.values])
        
#         # Mask the averaged timeframes for time occureneces below 95% of the max count occurence
#         mask_time = (count_occur < np.max(count_occur)*0.95)

#         #Drop the timeframes with occurences below the 95% of the max count occurence
#         signal = sig_avg.drop_sel(time = sig_avg.time.values[mask_time])
#         frames = frames.drop_sel(time = frames.time.values[mask_time])

#     else:
#     # If no timescale ('') is given then returns the raw temporal resolution
       
#         # calc the timedelta between timeframes to find raw temporal resolution   
#         deltat = [(time_arr[t+1]-time_arr[t])/np.timedelta64(1,'s') for t in range(time_arr.size-1)]

#         if len(deltat) > 0:
#             # minimum time delta -> raw temporal resolution (+- 1 sec)
#             t_resol = int(min(deltat))
#             adj_timescale = f'{t_resol}s'     
         
#         frames = xr.DataArray(data = np.ones(time_arr.size), 
#                   dims = ['time'], 
#                   coords = [time_arr])
#     # elif timescale == '' and len(time_arr) > 1:
#     # # If no timescale ('') is given then returns the raw temporal resolution
       
#     #     # calc the timedelta between timeframes to find raw temporal resolution   
#     #     deltat = [(time_arr[t+1]-time_arr[t])/np.timedelta64(1,'s') for t in range(time_arr.size-1)]
#     #     if len(deltat) > 0:
#     #         # minimum time delta -> raw temporal resolution (+- 1 sec)
#     #         t_resol = int(min(deltat))
#     #         adj_timescale = f'{t_resol}s'     
         
#     #         frames = xr.DataArray(data = np.ones(time_arr.size), 
#     #                   dims = ['time'], coords = [time_arr])
#     # else: 
#     #     frames = xr.DataArray(data = np.ones(time_arr.size), 
#     #                           dims = ['time'], coords = [time_arr])      

#     return(signal, frames, adj_timescale)

# #--f2.2 (implementing time resampling)    
# def average_by_time_resampling(signal, timescale):

#     adj_timescale = timescale

#     # average the signal if the timescale is not empty
#     if timescale != '':
#         # Sampling the measurements considering the time of the first measurement 
#         indicator = datetime_unit(timescale)
    
#         # #Set the base_val to start sampling with respect to the first measurement            
#         # base_val = int(signal.time[0].dt.strftime(f'%{indicator}').values)

#         count_occur = signal.time.resample(time=timescale, skipna=True, #base=base_val,
#                                             closed='left', label='left').count(dim='time').values

#         #if the measurement time is less than the given timescale, adjust the timescale
#         adj_timescale = check_timedelta(signal.time.values, timescale, indicator)        

#         # Mask the resampled timeframes for time occureneces below 75% of the max count occurence
#         mask_time = (count_occur < np.max(count_occur)*0.75)
        
#         # Time resampling and averaging of the signals inside each timeframe
#         signal = signal.resample(time=timescale, skipna=True, #base=base_val,
#                                  closed='left', label='left').mean(dim='time', skipna=True)
        
#         #Drop the timeframes with occurences below the 75% of the max count occurence
#         signal = signal.drop_sel(time = signal.time.values[mask_time])                                                                   
        
#     return(signal, adj_timescale)

# #--f3
# def change_resol(signal, info, resol):  

#     bins = signal.bins.values
#     sig = signal.values
#     step = info.resol.values
    
#     iters = signal.iters.values
#     time = signal.time.values
#     channel = signal.channel.values
        
#     if resol == '':
#         resol = np.max(info.resol.values)  
    
#     alt_l = bins[-1]*np.max(step)

#     # x_com = np.arange(0., alt_l + resol, resol)
#     x_com = np.arange(resol, alt_l + resol, resol)
#     y_com = np.nan*np.zeros([sig.shape[0], sig.shape[1], 
#                              sig.shape[2], x_com.size])
    
#     # interpolate the signal for common spatial resolution in all channels.
#     # for common raw resolution in all channels the signal remains intacted
#     for j in range(sig.shape[2]):
#         x = bins*step[j]
#         for k in range(sig.shape[1]):
#             for n in range(sig.shape[0]):
#                 y = signal[n, k, j, :]
#                 y_com[n, k, j, :] = np.interp(x_com, x, y)
    
    
#     #convert the signal xarray dim 'bins' to 'range' => dims [iters, time, channel, range]
#     signal = xr.DataArray(y_com, 
#                           coords=[iters, time, channel, x_com],
#                           dims=['iters', 'time', 'channel', 'range'])

#     info.resol[:] = resol

#     return(signal, info, resol)    

# #--f4
# def change_dim(signal, dim_name, new_dim_name, new_array):
#     if len(signal)>0:
#         signal[dim_name] = new_array
#         rename_dim = {dim_name:new_dim_name}
#         signal = signal.rename(rename_dim)
#         #signal[new_dim_name] = signal[new_dim_name].reset_index(new_dim_name)
#     return(signal)






# #--f8        
# def background_correction(signal, bc):
#     signal = signal - bc
#     return(signal)

# #--f9
# def range_correction(signal):
#     signal = signal*np.power(signal.range.values, 2)
#     return(signal)

# #--f10
# def dark_correction(signal, info, dark): 
    
#     an_ch = dict(channel = info.index[info.ch_mode == 0].values) #analog channel indexes
    
#     drk_mean = dark.mean(dim='time', skipna = True).loc[an_ch]
    
#     for t in range (signal.time.size):
#         ind_t = dict(time = t)      
#         signal[ind_t].loc[an_ch] = signal[ind_t].loc[an_ch] - drk_mean

#     return(signal)
  
# #--f12
# def construct_total(signal, info, tot):
#     '''

#     Parameters
#     ----------
#     signal : 3-D xarray dataset
#         Contains the signals with dimensions [time, channel, range].
#     info : object type dataframe
#         Contains info about each channel.
#     itot : integer
#         Indicates whether to construct the total signal.
#     tot : class with dataframes as attributes
#         Contains the info (map and parameters) from the total_map.ini file.


#     Returns
#     -------
#     None.

#     * Itot ~ cal_f*H_r*I_t - H_t*I_r
    
#     The formula is based on eq. 65 from Freudenthaler et al. 2016
#     We assume the caf_f is for the ratio R/T regardless which is the parallel/cross channel
#     '''
    
#     # extract the info from the total_map.ini file    
#     id_tot = tot.map.id_new
#     id_r   = tot.map.id_r
#     id_t   = tot.map.id_t

#     cal_f  = tot.params.cal_f
#     cal_f_err = tot.params.cal_f_err
#     H_r    = tot.params.h_r
#     H_r_err    = tot.params.h_r_err
#     H_t    = tot.params.h_t
#     H_t_err    = tot.params.h_t_err
        
#     for j in id_tot.index:
#         #Add new index in the dimension of channel for the total signals
#         signal = signal.reindex(channel = list(signal.channel.values)+[id_tot[j]])    
        
#         #Add the info about the total signal in the info_array
#         info.loc[id_tot[j]] = info.loc[id_r[j]].copy()
#         info.loc[id_tot[j], 'full_ovl_idx'] = \
#             max([info.loc[id_t[j], 'full_ovl_idx'],info.loc[id_r[j], 'full_ovl_idx']])
#         info.loc[id_tot[j], 'ch_pol'] = 'u' # for calibrated sum instead of 'o'

#         print(f'-> Calibrated sum: Adding channel {id_tot[j]}')


#         # Add statistical error to params
#         c_f = add_param_err(cal_f[j], cal_f_err[j], signal.iters.size) # 1D array with size = iters
#         H_R = add_param_err(H_r[j], H_r_err[j], signal.iters.size)
#         H_T = add_param_err(H_t[j], H_t_err[j], signal.iters.size)
        
#         #Construct the total signal for all time dims        
#         for n in range(signal.iters.size):
#             ind = dict(iters = n)
            
#             for t in range(signal.time.size):                
#                 ind['time'] = t
                
#                 temp_sig = c_f[n] * H_R[n] * signal[ind].loc[dict(channel = id_t[j])].values -\
#                     (H_T[n]*signal[ind].loc[dict(channel = id_r[j])].values)
                
#                 signal[ind].loc[dict(channel = id_tot[j])] = temp_sig


#     return(signal, info)

# #--f13
# def correct_total_sig(signal, info, prod, info_prod, vdr_mp):
#     '''
#     Correction of the detected total signal using the retrieved volume depol ratio profile.
    
#     The total signal is detected in the Reflected path
    
#     References
#     -----    
#     Engelmann et al., 2016, Atmos. Meas. Tech., https://doi.org/10.5194/amt-9-1767-2016
#     Mattis et al., 2009, Appl. Opt, https://doi.org/10.1364/AO.48.002742
    
#     '''
#     icor=False
#     if len(prod)>0 and len(vdr_mp)>0 and len(signal)>0:

#         n_size = signal.iters.size
    
#         for j in vdr_mp.index:

#             id_vldr = vdr_mp.prod_id.values[j]
#             id_tot = vdr_mp.ir_id.values[j]
            
#             # if vldr exists, have the same wl with signal amd signal is total then procced to the correction
#             if (id_vldr in prod.product.values) and (info.ch_pol.loc[id_tot] =='o') \
#                 and (float(info.wave.loc[id_tot])==float(info_prod.wave.loc[id_vldr])):
                
#                 icor=True
                
#                 tr_tot = float(vdr_mp.tr_r.values[j])
#                 tr_tot_err = float(vdr_mp.tr_r_err.values[j])
#                 TR_tot = add_param_err(tr_tot, tr_tot_err, n_size)
        
#                 v_ind = dict(product=id_vldr)
#                 s_ind = dict(channel=id_tot)
                
#                 for n in range(signal.iters.size):
#                     v_ind['iters' ]= n
#                     s_ind['iters' ]= n
#                     vdr_temp = prod.loc[v_ind].values
#                     sig_temp = signal.loc[s_ind].copy().values
                
#                     cor_factor = (TR_tot[n] * vdr_temp + 1) / (vdr_temp + 1) 
#                     signal.loc[s_ind] = sig_temp * cor_factor
            
#     return(signal, icor)

# #-- additional functions
