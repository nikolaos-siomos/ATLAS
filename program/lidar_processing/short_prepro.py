"""
@authors: N. Siomos & P. Paschou
"""

from lidar_processing import signal, diagnose
 
def rayleigh(sig_raw, shots, meas_info, channel_info, 
             external_info, time_info, sig_drk = []):
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
    
    isdk = external_info['skip_dark_subtraction']
    isdt = external_info['skip_dead_time_correction']
    itrv = external_info['vertical_trimming']
    
    alt_lim = external_info['vertical_limit']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    daq_range = channel_info.DAQ_Range
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
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
    # Detect saturation and potential clipping
    # --------------------------------------------------
    diagnose.detect_overflows(sig = sig.copy(), 
                              dead_time = dead_time, 
                              daq_range = daq_range)
    

    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if not isdt:
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
                                         timescale = -1)
    
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
    if itrv:
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
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')
    
    
    # --------------------------------------------------
    # Dark correction
    # --------------------------------------------------
    if not isdk and not isinstance(sig_drk,list):
        
        sig = signal.dark_correction(sig = sig.copy(), 
                                     drk = sig_drk.copy())
            
        if external_info['debug']: pack_out['sig_drc'] = sig.copy()
        
        print('-- Dark signal structure succesfully removed!')
        

    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')
    
    return(sig, pack_out, time_info)

def telecover(sig_raw, shots, meas_info, channel_info, external_info, 
              time_info, sig_drk = []):
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
    print('Start telecover signal preprocessing...')
    print('-----------------------------------------')
    
    isdk = external_info['skip_dark_subtraction']
    isdt = external_info['skip_dead_time_correction']
    itrv = external_info['vertical_trimming']
    
    alt_lim = external_info['vertical_limit']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    daq_range = channel_info.DAQ_Range
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
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
    # Detect saturation and potential clipping
    # --------------------------------------------------
    diagnose.detect_overflows(sig = sig.copy(), 
                              dead_time = dead_time, 
                              daq_range = daq_range)
    

    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if not isdt:
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
    sig, time_info = signal.average_by_group(sig = sig.copy(),
                                             time_info = time_info,
                                             start_time = 'Raw_Data_Start_Time',
                                             stop_time = 'Raw_Data_Stop_Time',
                                             grouper = 'sector')
    
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
    if itrv:
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
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')
    
    
    # --------------------------------------------------
    # Dark correction
    # --------------------------------------------------
    if not isdk and not isinstance(sig_drk,list):
        
        sig = signal.dark_correction(sig = sig.copy(), 
                                     drk = sig_drk.copy())
            
        if external_info['debug']: pack_out['sig_drc'] = sig.copy()
        
        print('-- Dark signal structure succesfully removed!')
        

    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')
        
    
    return(sig, pack_out, time_info)

def calibration(sig_raw, shots, meas_info, channel_info, external_info, 
                time_info, sig_drk = []):
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
    print('Start telecover signal preprocessing...')
    print('-----------------------------------------')
    
    isdk = external_info['skip_dark_subtraction']
    isdt = external_info['skip_dead_time_correction']
    itrv = external_info['vertical_trimming']
    
    alt_lim = external_info['vertical_limit']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    daq_range = channel_info.DAQ_Range
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
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
    # Detect saturation and potential clipping
    # --------------------------------------------------
    diagnose.detect_overflows(sig = sig.copy(), 
                              dead_time = dead_time, 
                              daq_range = daq_range)
    

    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if not isdt:
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
    sig, time_info = signal.average_by_group(sig = sig.copy(),
                                             time_info = time_info,
                                             start_time = 'Raw_Data_Start_Time',
                                             stop_time = 'Raw_Data_Stop_Time',
                                             grouper = 'calibrator_position')
    
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
    if itrv:
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
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')
    
    
    # --------------------------------------------------
    # Dark correction
    # --------------------------------------------------
    if not isdk and not isinstance(sig_drk,list):
        
        sig = signal.dark_correction(sig = sig.copy(), 
                                     drk = sig_drk.copy())
            
        if external_info['debug']: pack_out['sig_drc'] = sig.copy()
        
        print('-- Dark signal structure succesfully removed!')
        

    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')
        
    
    return(sig, pack_out, time_info)

def quicklook(sig_raw, shots, meas_info, channel_info, external_info, sig_drk):
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
    print('Start quicklook preprocessing...')
    print('-----------------------------------------')
    
    isdk = external_info['skip_dark_subtraction']
    isdt = external_info['skip_dead_time_correction']
    itrv = external_info['vertical_trimming']
    
    alt_lim = external_info['vertical_limit']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    daq_range = channel_info.DAQ_Range
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
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
    # Detect saturation and potential clipping
    # --------------------------------------------------
    diagnose.detect_overflows(sig = sig.copy(), 
                              dead_time = dead_time, 
                              daq_range = daq_range)
    

    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if not isdt:
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
    if itrv:
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
    # Solar backsground signal calculation
    # --------------------------------------------------
    sig = signal.background_correction(sig = sig.copy(), bgr = bgr.copy())
    
    if external_info['debug']: pack_out['sig_bgc'] = sig.copy()
    
    print('-- Solar background subtraction complete!')
    
    
    # --------------------------------------------------
    # Dark correction
    # --------------------------------------------------
    if not isdk and not isinstance(sig_drk,list):
        
        sig = signal.dark_correction(sig = sig.copy(), 
                                     drk = sig_drk.copy())
            
        if external_info['debug']: pack_out['sig_drc'] = sig.copy()
        
        print('-- Dark signal structure succesfully removed!')
        

    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')
        
    
    return(sig, pack_out)



def dark(sig_raw, shots, meas_info, channel_info, external_info, time_info):
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

    isdt = external_info['skip_dead_time_correction']
    itrv = external_info['vertical_trimming']
    
    alt_lim = external_info['vertical_limit']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    daq_range = channel_info.DAQ_Range
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = meas_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
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
    # Detect saturation and potential clipping
    # --------------------------------------------------
    diagnose.detect_overflows(sig = sig.copy(), 
                              dead_time = dead_time, 
                              daq_range = daq_range)
    

    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if not isdt:
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
                                          timescale = -1)
    
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
    if itrv:
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
                           smoothing_window = 500,
                           smoothing_sbin = 750,
                           smoothing_ebin = -1)
        
    if external_info['debug']: pack_out['sig_flt'] = sig.copy()
    
    print('-- Signal smoothing complete!')
        

    # --------------------------------------------------
    # Range correction
    # --------------------------------------------------
    sig = signal.range_correction(sig = sig, ranges = ranges)
    
    if external_info['debug']: pack_out['sig_rnc'] = sig.copy()
    
    print('-- Range correction complete!')
        
    
    return(sig, pack_out, time_info)

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

    # # --------------------------------------------------
    # # Calculation of the signals heights
    # # --------------------------------------------------
    # heights = signal.height_calculation(bins = sig.copy().bins.values, 
    #                                     ground_alt = ground_alt,
    #                                     zenith_angle = zenith_angle,                                       
    #                                     resol = resol)
    
    # pack_out['heights'] = heights.copy()
        
    # print('-- Heights calculated per signal bin and channel!')
