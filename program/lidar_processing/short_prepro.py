"""
@authors: N. Siomos & P. Paschou
"""

from lidar_processing import signal
 
def rayleigh(sig_raw, shots, meas_info, channel_info, time_info, 
             external_info, sig_drk = []):
    '''
    Perform the signal preprocessing for the rayleigh fit in the following order
    
     -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels
    
     -- dead time correction: Performs the dead time correction onphoton channels
     
     -- average_by_time: Average signals across the timeframes
     
     -- background_calculation: Calculates the solar background per timeframe and channel 
     
     -- trigger_delay: Perform the trigger correction per channel
     
     -- trim_vertically: Trim channels up to a maximum altitude
     
     -- height_calculation: Calculates the height values per bin and channel
     
     -- range calculation: Calculates the range values per bin and channel
     
     -- background_correction: Performs the background correction on signals

     -- dark_correction: Removes the dark signal structure from analog channels
     
     -- range_correction: Performs the range correction on signals
     
     -- smoothing: Smooths the signals (sliding average)

    Returns:
        
    - sig:
        A 4D xarray with dimensions (iters, time, channel, bins) containing the
        signal profiles after preprocessing
    
    - pack_out:
        A dictionary that contains at least the following:
            
         -- ranges (range per bin and channel)
            
         -- heights(range per bin and channel)
         
         -- bgr(background per averaged timeframe and channel)
        
        and also the signal in all preprocessing stages if the debug argument 
        is set to True
            
    '''   
        
    print('-----------------------------------------')
    print('Start rayleigh signal preprocessing...')
    print('-----------------------------------------')
    
    idrk = external_info['dark_subtraction']
    idtc = external_info['dead_time_correction']
    iflt = external_info['signal_smoothing']
    
    alt_lim = external_info['trim_vertically']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
    smoothing_window = external_info['smoothing_window']
    smoothing_sbin = external_info['smoothing_sbin']
    smoothing_ebin = external_info['smoothing_ebin']
    trd_bins = channel_info.Trigger_Delay_Bins
    zenith_angle = meas_info.Laser_Pointing_Angle
    
    sig = sig_raw.copy()
    
    pack_out = dict()

    # --------------------------------------------------
    # Unit conversion - raw counts to MHz for the photon channels
    # --------------------------------------------------
    sig = signal.unit_conv_counts_to_MHz(sig = sig.copy(), 
                                         shots = shots.copy(), 
                                         resol = resol)
    
    if external_info['debug']: pack_out['sig_puc'] = sig.copy()
    
    
    print('-- Unit conversion from raw counts to MHz for the pc channels complete!')


    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if idtc:
        sig = signal.dead_time_correction(sig = sig.copy(), 
                                          dead_time = dead_time, 
                                          dead_time_cor_type = dead_time_cor_type)

        if external_info['debug']: pack_out['sig_dtc'] = sig.copy()

        
        print('-- Dead time correction for pc channels complete!')
        
    else:
        print('-- Warning: Dead time correction for pc channels deactivated!')


    # --------------------------------------------------
    # Temporal averaging 
    # --------------------------------------------------   
    sig, frames = signal.average_by_time(sig = sig.copy(), 
                                         timescale = None)
    
    if external_info['debug']: pack_out['sig_avg'] = sig.copy()
    
    print(f'-- Temporal averaging complete! [..{sig.time.size} averaged timeframe(s)]')            
   

    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    bgr = signal.background_calculation(sig = sig.copy(), 
                                        lower_bin = bg_low,
                                        upper_bin = bg_high)

    pack_out['bgr'] = bgr.copy()

    print('-- Solar background succesfully calculated!')
    
    
    # --------------------------------------------------
    # Remove the pre-triggering region (or correct triger delays)
    # --------------------------------------------------
    sig = signal.trigger_delay(sig = sig.copy(), 
                               trigger_delay_bins = trd_bins)
    
    if external_info['debug']: pack_out['sig_trc'] = sig.copy()

    print('-- Triggering correction complete!')

    # --------------------------------------------------
    # Trim signal up to an upper altitude limit
    # --------------------------------------------------
    if alt_lim:
        sig = signal.trim_vertically(sig = sig.copy(), 
                                     ground_alt = ground_alt,
                                     zenith_angle = zenith_angle, 
                                     alt_lim = alt_lim,
                                     resol = resol)
    
        if external_info['debug']: pack_out['sig_trm'] = sig.copy()
        
        print(f'-- Range bins sucessfully trimmed above {alt_lim}m distance!')

    # --------------------------------------------------
    # Calculation of the signals ranges
    # --------------------------------------------------
    ranges = signal.range_calculation(bins = sig.copy().bins.values, 
                                      resol = resol)
    
    pack_out['ranges'] = ranges.copy()

        
    print('-- Ranges calculated per signal bin and channel!')


    # --------------------------------------------------
    # Calculation of the signals heights
    # --------------------------------------------------
    heights = signal.height_calculation(bins = sig.copy().bins.values, 
                                        ground_alt = ground_alt,
                                        zenith_angle = zenith_angle,                                       
                                        resol = resol)
    
    pack_out['heights'] = heights.copy()
        
    print('-- Heights calculated per signal bin and channel!')


    # # --------------------------------------------------
    # # Error simulation and xarray expansion to 4D
    # # --------------------------------------------------
    # if ierr:     
    #     if not isinstance(dark_pack,list): # normal signal with dark signal available
    #         sig, info_sig = error(sig.copy(), cfg, info_sig, frames, 
    #                               isdark, drk_inp = dark_avg.copy()) 

    #     if isinstance(dark_pack,list): # signal without dark signal availability
    #         sig, info_sig = error(sig.copy(), cfg, info_sig, frames, isdark) 

    #     sig_mc = sig.copy()
 
    #     print(f'-- Error simulations [{mc_sim}] and 4D array expansion complete!')


    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')


    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    if external_info['debug']: pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')


    # --------------------------------------------------
    # Smoothing
    # --------------------------------------------------
    if iflt == 1:
        sig = signal.smoothing(sig = sig, 
                               smoothing_window = smoothing_window,
                               smoothing_sbin = smoothing_sbin,
                               smoothing_ebin = smoothing_ebin)
            
        if external_info['debug']: pack_out['sig_flt'] = sig.copy()
        
        print('-- Signal smoothing complete!')

    # --------------------------------------------------
    # Dark correction
    # --------------------------------------------------
    if idrk and not isinstance(sig_drk,list):
        
        sig = signal.dark_correction(sig = sig.copy(), 
                                     drk = sig_drk.copy())
            
        if external_info['debug']: pack_out['sig_drc'] = sig.copy()
        
        print('-- Dark signal structure succesfully removed!')
        
    
    return(sig, pack_out)

def telecover(sig_raw, shots, meas_info, channel_info, time_info, 
              external_info, sig_drk = []):
    '''
    Perform the signal preprocessing for the rayleigh fit in the following order
    
     -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels
    
     -- dead time correction: Performs the dead time correction onphoton channels
          
     -- background_calculation: Calculates the solar background per timeframe and channel 
     
     -- trigger_delay: Perform the trigger correction per channel
     
     -- trim_vertically: Trim channels up to a maximum altitude
     
     -- height_calculation: Calculates the height values per bin and channel
     
     -- range calculation: Calculates the range values per bin and channel
     
     -- background_correction: Performs the background correction on signals

     -- dark_correction: Removes the dark signal structure from analog channels
     
     -- range_correction: Performs the range correction on signals
     
     -- smoothing: Smooths the signals (sliding average)

    Returns:
        
    - sig:
        A 4D xarray with dimensions (iters, time, channel, bins) containing the
        signal profiles after preprocessing
    
    - pack_out:
        A dictionary that contains at least the following:
            
         -- ranges (range per bin and channel)
            
         -- heights(range per bin and channel)
         
         -- bgr(background per averaged timeframe and channel)
        
        and also the signal in all preprocessing stages if the debug argument 
        is set to True
            
    '''   
        
    print('-----------------------------------------')
    print('Start rayleigh signal preprocessing...')
    print('-----------------------------------------')
    
    idrk = external_info['dark_subtraction']
    idtc = external_info['dead_time_correction']
    iflt = external_info['signal_smoothing']

    
    alt_lim = external_info['trim_vertically']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
    smoothing_window = external_info['smoothing_window']
    smoothing_sbin = external_info['smoothing_sbin']
    smoothing_ebin = external_info['smoothing_ebin']
    trd_bins = channel_info.Trigger_Delay_Bins
    zenith_angle = meas_info.Laser_Pointing_Angle
    
    sig = sig_raw.copy()
    
    pack_out = dict()

    # --------------------------------------------------
    # Unit conversion - raw counts to MHz for the photon channels
    # --------------------------------------------------
    sig = signal.unit_conv_counts_to_MHz(sig = sig.copy(), 
                                         shots = shots.copy(), 
                                         resol = resol)
    
    if external_info['debug']: pack_out['sig_puc'] = sig.copy()
    
    
    print('-- Unit conversion from raw counts to MHz for the pc channels complete!')


    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if idtc:
        sig = signal.dead_time_correction(sig = sig.copy(), 
                                          dead_time = dead_time, 
                                          dead_time_cor_type = dead_time_cor_type)

        if external_info['debug']: pack_out['sig_dtc'] = sig.copy()

        
        print('-- Dead time correction for pc channels complete!')
        
    else:
        print('-- Warning: Dead time correction for pc channels deactivated!')


    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    bgr = signal.background_calculation(sig = sig.copy(), 
                                        lower_bin = bg_low,
                                        upper_bin = bg_high)

    pack_out['bgr'] = bgr.copy()

    print('-- Solar background succesfully calculated!')
    
    
    # --------------------------------------------------
    # Remove the pre-triggering region (or correct triger delays)
    # --------------------------------------------------
    sig = signal.trigger_delay(sig = sig.copy(), 
                               trigger_delay_bins = trd_bins)
    
    if external_info['debug']: pack_out['sig_trc'] = sig.copy()

    print('-- Triggering correction complete!')

    # --------------------------------------------------
    # Trim signal up to an upper altitude limit
    # --------------------------------------------------
    if alt_lim:
        sig = signal.trim_vertically(sig = sig.copy(), 
                                     ground_alt = ground_alt,
                                     zenith_angle = zenith_angle, 
                                     alt_lim = alt_lim,
                                     resol = resol)
    
        if external_info['debug']: pack_out['sig_trm'] = sig.copy()
        
        print(f'-- Range bins sucessfully trimmed above {alt_lim}m distance!')

    # --------------------------------------------------
    # Calculation of the signals ranges
    # --------------------------------------------------
    ranges = signal.range_calculation(bins = sig.copy().bins.values, 
                                      resol = resol)
    
    pack_out['ranges'] = ranges.copy()

        
    print('-- Ranges calculated per signal bin and channel!')


    # --------------------------------------------------
    # Calculation of the signals heights
    # --------------------------------------------------
    heights = signal.height_calculation(bins = sig.copy().bins.values, 
                                        ground_alt = ground_alt,
                                        zenith_angle = zenith_angle,                                       
                                        resol = resol)
    
    pack_out['heights'] = heights.copy()
        
    print('-- Heights calculated per signal bin and channel!')


    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')


    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    if external_info['debug']: pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')


    # --------------------------------------------------
    # Smoothing
    # --------------------------------------------------
    if iflt == 1:
        sig = signal.smoothing(sig = sig, 
                               smoothing_window = smoothing_window,
                               smoothing_sbin = smoothing_sbin,
                               smoothing_ebin = smoothing_ebin)
            
        if external_info['debug']: pack_out['sig_flt'] = sig.copy()
        
        print('-- Signal smoothing complete!')

    # --------------------------------------------------
    # Dark correction
    # --------------------------------------------------
    if idrk and not isinstance(sig_drk,list):
        
        sig = signal.dark_correction(sig = sig.copy(), 
                                     drk = sig_drk.copy())
            
        if external_info['debug']: pack_out['sig_drc'] = sig.copy()
        
        print('-- Dark signal structure succesfully removed!')
        
    
    return(sig, pack_out)

def quicklook(sig_raw, shots, meas_info, channel_info, time_info, 
              external_info, sig_drk):
    '''
    Perform the signal preprocessing for quicklooks in the following order
    
     -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels
    
     -- dead time correction: Performs the dead time correction onphoton channels
          
     -- background_calculation: Calculates the solar background per timeframe and channel 
     
     -- trigger_delay: Perform the trigger correction per channel
     
     -- trim_vertically: Trim channels up to a maximum altitude
     
     -- height_calculation: Calculates the height values per bin and channel
     
     -- range calculation: Calculates the range values per bin and channel
     
     -- background_correction: Performs the background correction on signals

     -- dark_correction: Removes the dark signal structure from analog channels
     
     -- range_correction: Performs the range correction on signals
     
     -- smoothing: Smooths the signals (sliding average)

    Returns:
        
    - sig:
        A 4D xarray with dimensions (iters, time, channel, bins) containing the
        signal profiles after preprocessing
    
    - pack_out:
        A dictionary that contains at least the following:
            
         -- ranges (range per bin and channel)
            
         -- heights(range per bin and channel)
         
         -- bgr(background per averaged timeframe and channel)
        
        and also the signal in all preprocessing stages if the debug argument 
        is set to True
            
    '''   
        
    print('-----------------------------------------')
    print('Start rayleigh signal preprocessing...')
    print('-----------------------------------------')
    
    idrk = external_info['dark_subtraction']
    idtc = external_info['dead_time_correction']
    iflt = external_info['signal_smoothing']

    
    alt_lim = external_info['trim_vertically']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
    smoothing_window = external_info['smoothing_window']
    smoothing_sbin = external_info['smoothing_sbin']
    smoothing_ebin = external_info['smoothing_ebin']
    trd_bins = channel_info.Trigger_Delay_Bins
    zenith_angle = meas_info.Laser_Pointing_Angle
    
    sig = sig_raw.copy()
    
    pack_out = dict()

    # --------------------------------------------------
    # Unit conversion - raw counts to MHz for the photon channels
    # --------------------------------------------------
    sig = signal.unit_conv_counts_to_MHz(sig = sig.copy(), 
                                         shots = shots.copy(), 
                                         resol = resol)
    
    if external_info['debug']: pack_out['sig_puc'] = sig.copy()
    
    
    print('-- Unit conversion from raw counts to MHz for the pc channels complete!')


    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if idtc:
        sig = signal.dead_time_correction(sig = sig.copy(), 
                                          dead_time = dead_time, 
                                          dead_time_cor_type = dead_time_cor_type)

        if external_info['debug']: pack_out['sig_dtc'] = sig.copy()

        
        print('-- Dead time correction for pc channels complete!')
        
    else:
        print('-- Warning: Dead time correction for pc channels deactivated!')


    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    bgr = signal.background_calculation(sig = sig.copy(), 
                                        lower_bin = bg_low,
                                        upper_bin = bg_high)

    pack_out['bgr'] = bgr.copy()

    print('-- Solar background succesfully calculated!')
    
    
    # --------------------------------------------------
    # Remove the pre-triggering region (or correct triger delays)
    # --------------------------------------------------
    sig = signal.trigger_delay(sig = sig.copy(), 
                               trigger_delay_bins = trd_bins)
    
    if external_info['debug']: pack_out['sig_trc'] = sig.copy()

    print('-- Triggering correction complete!')

    # --------------------------------------------------
    # Trim signal up to an upper altitude limit
    # --------------------------------------------------
    if alt_lim:
        sig = signal.trim_vertically(sig = sig.copy(), 
                                     ground_alt = ground_alt,
                                     zenith_angle = zenith_angle, 
                                     alt_lim = alt_lim,
                                     resol = resol)
    
        if external_info['debug']: pack_out['sig_trm'] = sig.copy()
        
        print(f'-- Range bins sucessfully trimmed above {alt_lim}m distance!')

    # --------------------------------------------------
    # Calculation of the signals ranges
    # --------------------------------------------------
    ranges = signal.range_calculation(bins = sig.copy().bins.values, 
                                      resol = resol)
    
    pack_out['ranges'] = ranges.copy()

        
    print('-- Ranges calculated per signal bin and channel!')


    # --------------------------------------------------
    # Calculation of the signals heights
    # --------------------------------------------------
    heights = signal.height_calculation(bins = sig.copy().bins.values, 
                                        ground_alt = ground_alt,
                                        zenith_angle = zenith_angle,                                       
                                        resol = resol)
    
    pack_out['heights'] = heights.copy()
        
    print('-- Heights calculated per signal bin and channel!')


    # # --------------------------------------------------
    # # Error simulation and xarray expansion to 4D
    # # --------------------------------------------------
    # if ierr:     
    #     if not isinstance(dark_pack,list): # normal signal with dark signal available
    #         sig, info_sig = error(sig.copy(), cfg, info_sig, frames, 
    #                               isdark, drk_inp = dark_avg.copy()) 

    #     if isinstance(dark_pack,list): # signal without dark signal availability
    #         sig, info_sig = error(sig.copy(), cfg, info_sig, frames, isdark) 

    #     sig_mc = sig.copy()
 
    #     print(f'-- Error simulations [{mc_sim}] and 4D array expansion complete!')


    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')


    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    if external_info['debug']: pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')


    # --------------------------------------------------
    # Smoothing
    # --------------------------------------------------
    if iflt == 1:
        sig = signal.smoothing(sig = sig, 
                               smoothing_window = smoothing_window,
                               smoothing_sbin = smoothing_sbin,
                               smoothing_ebin = smoothing_ebin)
            
        if external_info['debug']: pack_out['sig_flt'] = sig.copy()
        
        print('-- Signal smoothing complete!')

    # --------------------------------------------------
    # Dark correction
    # --------------------------------------------------
    if idrk and not isinstance(sig_drk,list):
        
        sig = signal.dark_correction(sig = sig.copy(), 
                                     drk = sig_drk.copy())
            
        if external_info['debug']: pack_out['sig_drc'] = sig.copy()
        
        print('-- Dark signal structure succesfully removed!')
        
    
    return(sig, pack_out)



def dark(sig_raw, shots, meas_info, channel_info, time_info, 
         external_info):
    '''
    Perform the dark signal preprocessing in the following order
     
     -- average_by_time: Average signals across the timeframes
     
     -- background_calculation: Calculates the solar background per timeframe and channel 
     
     -- trigger_delay: Perform the trigger correction per channel
     
     -- trim_vertically: Trim channels up to a maximum altitude
     
     -- height_calculation: Calculates the height values per bin and channel
     
     -- range calculation: Calculates the range values per bin and channel
     
     -- background_correction: Performs the background correction on signals
     
     -- range_correction: Performs the range correction on signals
     
     -- smoothing: Smooths the signals (sliding average)

    Returns:
        
    - sig:
        A 4D xarray with dimensions (iters, time, channel, bins) containing the
        signal profiles after preprocessing
    
    - pack_out:
        A dictionary that contains at least the following:
            
         -- ranges (range per bin and channel)
            
         -- heights(range per bin and channel)
         
         -- bgr(background per averaged timeframe and channel)
        
        and also the signal in all preprocessing stages if the debug argument 
        is set to True
            
    '''   
        
    print('-----------------------------------------')
    print('Start dark signal preprocessing...')
    print('-----------------------------------------')
    
    alt_lim = external_info['trim_vertically']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
    timescale = external_info['timescale']
    trd_bins = channel_info.Trigger_Delay_Bins
    zenith_angle = meas_info.Laser_Pointing_Angle
    
    sig = sig_raw.copy()
    
    pack_out = dict()

    # --------------------------------------------------
    # Temporal averaging 
    # --------------------------------------------------   
    sig, frames = signal.average_by_time(sig = sig.copy(), 
                                          timescale = timescale)
    
    if external_info['debug']: pack_out['sig_avg'] = sig.copy()
    
    print(f'-- Temporal averaging complete! [..{sig.time.size} averaged timeframe(s)]')            
   
  
    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    bgr = signal.background_calculation(sig = sig.copy(), 
                                        lower_bin = bg_low,
                                        upper_bin = bg_high)

    pack_out['bgr'] = bgr.copy()

    print('-- Solar background succesfully calculated!')
    
    
    # --------------------------------------------------
    # Remove the pre-triggering region (or correct triger delays)
    # --------------------------------------------------
    sig = signal.trigger_delay(sig = sig.copy(), 
                               trigger_delay_bins = trd_bins)
    
    if external_info['debug']: pack_out['sig_trc'] = sig.copy()

    print('-- Triggering correction complete!')

    # --------------------------------------------------
    # Trim signal up to an upper altitude limit
    # --------------------------------------------------
    if alt_lim:
        sig = signal.trim_vertically(sig = sig.copy(), 
                                     ground_alt = ground_alt,
                                     zenith_angle = zenith_angle, 
                                     alt_lim = alt_lim,
                                     resol = resol)
    
        if external_info['debug']: pack_out['sig_trm'] = sig.copy()
        
        print(f'-- Range bins sucessfully trimmed above {alt_lim}m distance!')

    # --------------------------------------------------
    # Calculation of the signals ranges
    # --------------------------------------------------
    ranges = signal.range_calculation(bins = sig.copy().bins.values, 
                                      resol = resol)
    
    pack_out['ranges'] = ranges.copy()

        
    print('-- Ranges calculated per signal bin and channel!')


    # --------------------------------------------------
    # Calculation of the signals heights
    # --------------------------------------------------
    heights = signal.height_calculation(bins = sig.copy().bins.values, 
                                        ground_alt = ground_alt,
                                        zenith_angle = zenith_angle,                                       
                                        resol = resol)
    
    pack_out['heights'] = heights.copy()
        
    print('-- Heights calculated per signal bin and channel!')


    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')


    # --------------------------------------------------
    # Smoothing
    # --------------------------------------------------
    sig = signal.smoothing(sig = sig, 
                           smoothing_window = 1000,
                           smoothing_sbin = 500,
                           smoothing_ebin = None)
        
    if external_info['debug']: pack_out['sig_flt'] = sig.copy()
    
    print('-- Signal smoothing complete!')
        

    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    if external_info['debug']: pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')
        
    
    return(sig, pack_out)

    # # --------------------------------------------------
    # # Error simulation and xarray expansion to 4D
    # # --------------------------------------------------
    # if ierr:     
    #     if not isinstance(dark_pack,list): # normal signal with dark signal available
    #         sig, info_sig = error(sig.copy(), cfg, info_sig, frames, 
    #                               isdark, drk_inp = dark_avg.copy()) 

    #     if isinstance(dark_pack,list): # signal without dark signal availability
    #         sig, info_sig = error(sig.copy(), cfg, info_sig, frames, isdark) 

    #     sig_mc = sig.copy()
 
    #     print(f'-- Error simulations [{mc_sim}] and 4D array expansion complete!')

