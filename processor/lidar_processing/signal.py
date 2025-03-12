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
 -- height_calculation: Calculates the height above the lidar values per bin and channel
 -- range calculation: Calculates the range above the lidar values per bin and channel
 -- range_correction: Performs the range correction on signals
 -- smoothing: Smooths the signals (sliding average)
 -- trigger_correction: Perform the trigger correction per channel
 -- trim_vertically: Trim channels up to a maximum altitude
 -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels

"""

import numpy as np

import xarray as xr

import pandas as pd

def average_by_group(sig, time_info, grouper, start_time, stop_time):

    """
    General:
        Averages the lidar signals across the time deminsion 
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the following
            dimensions: (time, channel, bins). The photon signals must not be
            dead time corrected yet
        
        time_info : 
                A pandas Dataframe with all the necessary metadata for each 
                temporal frame fetched from the QA file
                
        grouper: 
            A string with the column name of time_info that contains the 
            flag used for grouping the data (e.g. telecover sector)

        start_time: 
            A string with the column name of time_info that contains the 
            starting time in seconds since the beginning of the measurement 
            per timeframe
        
        stop_time: 
            A string with the column name of time_info that contains the 
            ending time in seconds since the beginning of the measurement 
            per timeframe
              
    Returns:
        sig_out: 
            An xarray in the same dimensions with sig but with the time dim
            reduced by averaging. An average is generated each time the group
            changes
            
        time_info_out : 
                A pandas Dataframe with all the necessary metadata for each 
                temporal frame reformed according to the new temporal scale
                
            
    """
    
    grp = time_info.loc[:,grouper].values
    stime = time_info.loc[:,start_time].values
    etime = time_info.loc[:,stop_time].values

    start_dt = np.datetime64(sig.time.values[0], 's')

    sep = np.hstack([0, np.where(grp[:-1] != grp[1:])[0] + 1, grp.size])

    group_stime = np.array([stime[sep[i]] for i in range(0,sep.size-1)])
    group_etime = np.array([etime[sep[i]-1] for i in range(1,sep.size)])
    group_dtime = np.floor(group_stime + \
                           (group_etime - group_stime)/2.).astype(int)
    group_mtime = np.array([start_dt + np.timedelta64(gr_t, 's')
                            for gr_t in group_dtime])
    frames = np.array([sep[i+1] - sep[i] for i in range(0,sep.size-1)])
    
    group_grouper = np.array([grp[sep[i]] for i in range(0,sep.size-1)])
    
    sig_out = np.nan * sig[dict(time = range(len(group_mtime)))].copy()\
        .drop("time").assign_coords({'time' : group_mtime})
    
    for i in range(sep.size-1):
        t_d = dict(time = slice(sep[i], sep[i+1]))
        sig_out[dict(time = i)] =  sig[t_d].copy().mean(dim = 'time').values
    
    cols = [start_time, stop_time, grouper, 'Timeframes']

    data = np.array([group_stime, group_etime, group_grouper, frames],
                    dtype = object).T

    time_info_out = pd.DataFrame(data, index = group_mtime, columns = cols,
                                 dtype = object)
    
    return(sig_out, time_info_out)

def average_by_time(sig, time_info, timescale, start_time, stop_time):

    """
    General:
        Averages the lidar signals across the time deminsion 
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the following
            dimensions: (time, channel, bins). The photon signals must not be
            dead time corrected yet

        time_info : 
                A pandas Dataframe with all the necessary metadata for each 
                temporal frame fetched from the QA file
                
        timescale: 
            The temporal avraging step, in seconds. E.g. if set to 3600 then 1h
            averages will be generated from the original signals. If set to -1
            (default) then all timeframes will be averaged
            
        start_time: 
            A string with the column name of time_info that contains the 
            starting time in seconds since the beginning of the measurement 
            per timeframe
        
        stop_time: 
            A string with the column name of time_info that contains the 
            ending time in seconds since the beginning of the measurement 
            per timeframe

              
    Returns:
        
        sig_out: 
            An xarray in the same dimensions with sig but with the time dim
            reduced by averaging, according to the provided timescale
        
        time_info_out : 
                A pandas Dataframe with all the necessary metadata for each 
                temporal frame reformed according to the new temporal scale
                
            
    """
    
    time_arr = sig.time.values
    stime = time_info.loc[:,start_time].values
    etime = time_info.loc[:,stop_time].values
    
    # Calculate the time duration (hours) between first and last measurement
    dt = time_arr[-1] - time_arr[0]
    
    meas_dur = np.ceil(dt / np.timedelta64(1,'s'))
    
    if timescale == -1:
        timescale = meas_dur
    
    # average the sig only if the timescale is not empty and the sig time dim is >1
    if timescale >= meas_dur:
                
        m_time = np.datetime64(time_arr[0] + dt / 2.,'us').item()
        
        sig_avg = sig.mean(dim = 'time').values
        
        sig_out = sig.reindex(time = [m_time], fill_value = np.nan)

        sig_out.loc[dict(time = m_time)] = sig_avg
        
        # Creat the time_info object according to the new timescale
        cols = [start_time, stop_time, 'Timeframes']

        group_stime = np.array([stime[0]])
        group_etime = np.array([etime[-1]])
        group_count = np.array([time_arr.size])
        
        data = np.array([group_stime, group_etime, group_count],
                        dtype = object).T
    
        time_info_out = pd.DataFrame(data, 
                                     index = [m_time],
                                     columns = cols,
                                     dtype = object)
        

    else:
        halfscale_m = int(timescale / 2)

        timescale_m = int(timescale)

        # Define how many averaged timeframes will be created
        avg_num = int(np.ceil(meas_dur/timescale))
        
        # Calculate the new time index based on the timescale of averaging.
        # -> It is also the middle time in the time slice for averaging
        m_time = [(time_arr[0] + ii * np.timedelta64(halfscale_m,'s')) 
                  for ii in range(avg_num)]

        # Calculate the starting time in the time slice for averaging
        s_time = [(time_arr[0] + ii * np.timedelta64(timescale_m,'s')) 
                  for ii in range(avg_num)]

        # Calculate the ending time in the time slice for averaging
        e_time = [(time_arr[0] + (ii + 1) * np.timedelta64(timescale_m,'s')) 
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
            
        # Creat the time_info object according to the new timescale
        cols = [start_time, stop_time, 'Timeframes']

        group_stime = np.array([stime[np.sum(count_occur[:i])] for i in range(0,count_occur.size-1)])
        group_etime = np.array([etime[np.sum(count_occur[:i])-1] for i in range(1,count_occur.size)])
        
        data = np.array([group_stime, group_etime, count_occur],
                        dtype = object).T
    
        time_info_out = pd.DataFrame(data, index = m_time, columns = cols,
                                     dtype = object)
        
        # Mask the averaged timeframes for time occureneces below 75% of the max count occurence
        mask_time = (count_occur < np.max(count_occur) * 0.75)

        #Drop the timeframes with occurences below the 75% of the max count occurence
        sig_out = sig_out.drop_sel(time = sig_out.time.values[mask_time])
        
        #Drop the timeframes with occurences below the 75% of the max count occurence
        time_info_out = time_info_out[mask_time,:]

    return(sig_out, time_info_out)

def background_calculation(sig, lower_bin, upper_bin):

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

def background_correction(sig, bgr):
    
    """
    General:
        Performs the solar background subtraction 
        
    Input:
        sig: 
            An xarray with the lidar signals. It must include a time demension
            
        bgr: 
            An xarray with the background signals. It must be the same 
            dimensions as sig, excluding only time
              
    Returns:
        
        sig_out: 
            An xarray in the same dimensions with sig with the background
            corrected signals
            
    """
    
    sig_out = (sig.copy() - bgr.copy()).transpose('time','channel','bins')
    
    return(sig_out)

def dark_correction(sig, drk): 
   
    """
    General:
        Removes the dark signal from the normal signal for all analog channels
        
    Input:
        sig: 
            An xarray with the lidar signals, it should include the 
            following dimensions: (..., time, channel, ...). 
            
        drk: 
            An xarray with the dark signals, it should include the 
            following dimensions: (..., time, channel, ...). It should have 
            the same coordinates as sig with the exception of the time
            dimension             

    Returns:
        
        sig_out: 
            An xarray in the same coordinates as sig with the dark corrected
            signals
            
    """
            
    drk_mean = drk.mean(dim='time', skipna = True).copy()

    channels = sig.channel.values
    
    sig_out = sig.copy()
    
    for ch in channels:

        ch_d = dict(channel = ch)

        if ch[6] == 'a': #7th digit of channel name is the acquisition mode (a or p)
     
            sig_out.loc[ch_d] = sig_out.loc[ch_d] - drk_mean.loc[ch_d]

            
    return(sig_out)

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

        if ch[6] == 'p': #7th digit of channel name is the acquisition mode (a or p)

            if (sig.loc[ch_d].values > 1000./dead_time[ch]).any():
                print(f"-- Warning: Channel {ch} - A countrate value above the maximum allowed value was detected! Consider revising the input file ")

            if dead_time_cor_type[ch] in [0,1]: # Change after including non-paralyzable correction
                sig_out.loc[ch_d] = sig.loc[ch_d] / \
                    (1. - sig.loc[ch_d] * dead_time.loc[ch] * 1e-3)

        if ch[6] == 'a': #7th digit of channel name is the acquisition mode (a or p)
            
            sig_out.loc[ch_d] = sig.loc[ch_d]
            
    return(sig_out)

def height_calculation(bins, resol, zenith_angle):
 
    """
    General:
        Creates an xarray with the height (altitude) values per bin and channel
        
    Input:
        bins: 
            An 1D numpy array with the signal bins
            
        zenith_angle:
            A float with the zenith angle where the lidar points to
            
        resol:
            A pandas series with the range resolution in m per channel. 
            The index should correspond to the channel dimension of sig
            
    Returns:
        
        heights: 
            An xarray with the height values per bin and channel in meters
            
    """
    
    bins = xr.DataArray(bins, 
                        dims = ['bins'], 
                        coords = [bins])
    
    resolution = xr.DataArray(resol.values,
                              dims = ['channel'],
                              coords = [resol.index.values]).astype(float)
    
    zenith_angle_rad = np.deg2rad(zenith_angle)
    
    heights = resolution * (bins + 0.5) * np.cos(zenith_angle_rad)
    
    return(heights)

def range_calculation(bins, resol):
 
    """
    General:
        Creates an xarray with the altitude values per bin and channel
        
    Input:
        bins: 
            An 1D numpy array with the signal bins
                
        resol:
            A pandas series with the range resolution in m per channel. 
            The index should correspond to the channel dimension of sig
            
    Returns:
        
        ranges: 
            An xarray with the range values per bin and channel in meters
            
    """
    
    bins = xr.DataArray(bins, 
                        dims = ['bins'], 
                        coords = [bins])
    
    resolution = xr.DataArray(resol.values,
                              dims = ['channel'],
                              coords = [resol.index.values]).astype(float)
    ranges = resolution * (bins + 0.5)

    return(ranges)

def range_correction(sig, ranges):

    """
    General:
        Performs the range correction 
        
    Input:
        sig: 
            An xarray with the lidar signals. It must include at least the 
            following dimensions (channel, bins)
            
        ranges: 
            An xarray with the ranges per channel and bin. It must include at 
            least the following dimensions (channel, bins)
              
    Returns:
        
        sig_out: 
            An xarray in the same dimensions with sig with the range
            corrected signals
            
    """
    
    sig_out = (sig.copy() * np.power(ranges.copy(), 2))\
        .transpose('time','channel','bins')
    
    return(sig_out)

def smoothing(sig, smoothing_window, smoothing_sbin, smoothing_ebin):

    """
    General:
        A simple sliding average smoothing routine
        
    Input:
        sig: 
            An xarray with the lidar signals. It must include at least the 
            following dimensions (bins)
            
        window: 
            The number of bins that will be used for the averaging
        
        smoothing_sbin:
            An integer with the starting bin for smoothing. No smoothing will 
            be applied before it. If set to 0 the first bin will be used
            
        smoothing_ebin:
            An integer with the ending bin for smoothing. No smoothing will 
            be applied after it. If set to -1 the last bin will be used
              
    Returns:
        
        sig_out: 
            An xarray in the same dimensions with sig with the smoothed signals
            
    """   
    
    sig_out = sig.copy()
    std_out = sig.copy()
        
    if smoothing_ebin == -1:
        smoothing_ebin = sig.bins.size

    bins_d = dict(bins = slice(smoothing_sbin, smoothing_ebin))

    sig_out.loc[bins_d] = sig.copy()\
        .rolling(bins = smoothing_window, center = True).mean().loc[bins_d]

    std_out.loc[bins_d] = sig.copy()\
        .rolling(bins = smoothing_window, center = True).std().loc[bins_d]
                
    return(sig_out, std_out)

# def trim_clouds(sig, daq_trigger_offset):

#     time = sig.time.values
#     bins = sig.bins.values    
#     channels = sig.channel.values

#     mask_clouds = np.zeros()
    
#     skip_bins = 100
    
#     diff_bins = 25
    
#     smooth_bins = 10
    
#     sig_sm = 
    
#     for ch in channels:
        
#         ch_d = dict(channel = ch)
        
#         zero_bins = daq_trigger_offset.loc[ch]
        
#         start = dict(bins = slice(zero_bins + skip_bin ,bins.size - ))
        
#         delta_sig = sig.loc[].copy() - \
#             sig.loc[dict(bins = slice(100,10000))].copy()
        
#     delta_sig = np.abs(delta_sig.max(dim = 'bins')) / np.abs(delta_sig.min(dim = 'bins'))
    
#     mask_clouds = (delta_sig > 100.).any(dim = 'channel').values
    
#     sig = sig.copy().loc[dict(time = time[~mask_clouds])]
                         
#     return(sig)

def trigger_correction(sig, daq_trigger_offset):

    """
    General:
        Shift the bin axis to the left (pre-trigger) or to the 
        right (trigger delay) in order to bring all channels to the same 
        ground bin. Empty bins are filled with nans
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, bins). 
            
        daq_trigger_offset: 
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
        
        shift_bins = daq_trigger_offset.loc[ch]
        
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
            A float with the maximum altitude in m above which no computations 
            will be performed
            
        resol:
            A pandas series with the range resolution in m per channel. 
            The index should correspond to the channel dimension of sig
            
    Returns:
        
        sig_out: 
            An xarray in the same shape as sig with the signals trimmed above
            the maximum allowed altitude
            
    """
       
    sig_out = np.nan * sig.copy()
    
    zenith_angle_rad = np.deg2rad(zenith_angle)

    rng_lim = (float(alt_lim) - float(ground_alt))/np.cos(zenith_angle_rad)
    
    bin_lim = (rng_lim / np.min(resol)).astype(int)
    
    bins = sig.bins.values
    
    bin_d = dict(bins = slice(bins[0], bin_lim))

    if bin_lim < np.max(bins):
        
        sig_out = sig.loc[bin_d]
    
    else:
        sig_out = sig.copy()
  
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
    
        if ch[6] == 'p': #7th digit of channel name is the acquisition mode (a or p)
            
            sampl_rate = 150. / resol.loc[ch]

            sig_out.loc[ch_d] = sig.loc[ch_d] * sampl_rate / shots.loc[:,ch]
        
        if ch[6] == 'a': #7th digit of channel name is the acquisition mode (a or p)
            
            sig_out.loc[ch_d] = sig.loc[ch_d]
    
    return(sig_out) 
