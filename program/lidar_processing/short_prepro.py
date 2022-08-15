"""
@authors: N. Siomos & P. Paschou
"""

from lidar_processing import signal
 
def rayleigh(sig_raw, shots, lidar_info, channel_info, time_info, 
              external_info, dark_pack = []):
    '''
    Perform the pre-processing in the signals
    - unit conversion (photon to MHz)
    - dead-time correction for pc channels
    - time averaging
    - simulate errors if specified
    - background subtraction
    - Trigger delay correction
    - Change of spatial resolution (optional)
    - Signal trimming (optional)
    - Removal of the dark signal structure (optional; only for normal signals)
    - Range correction
    - Smoothing

    
    Returns:
    - signals(xarrays) after each correction with 4-dimensions[iters, time, channel, range/bin]
    ''' 

    
    if isinstance(dark_pack,dict) == True:
        # smoothed rc sig
        dark = dark_pack['sig_fltr']
        # time averaged sig
        dark_avg = dark_pack['sig_avg']
    else:
        dark = []
        dark_avg = []    
        
    print('-----------------------------------------')
    print('Start signal preprocessing...')
    print('-----------------------------------------')
    
    ierr = external_info['error_simulation']
    idtc = external_info['dead_time_correction']
    idrk = external_info['dark_subtraction']
    iflt = external_info['signal_smoothing']
    iavg = external_info['signal_averaging']

    
    alt_lim = external_info['trim_vertically']
    bg_low = channel_info.Background_Low
    bg_high = channel_info.Background_High
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    ground_alt = lidar_info.Altitude_meter_asl
    resol = channel_info.Raw_Data_Range_Resolution
    timescale = external_info['timescale']
    trd_bins = channel_info.Trigger_Delay_Bins
    zenith_angle = lidar_info.Laser_Pointing_Angle
    
    sig = sig_raw.copy()

    # --------------------------------------------------
    # Unit conversion - raw counts to MHz for the photon channels
    # --------------------------------------------------
    sig = signal.unit_conv_counts_to_MHz(sig = sig.copy(), 
                                         shots = shots.copy(), 
                                         resol = resol)
    sig_puc = sig.copy()
    
    print('-- Unit conversion from raw counts to MHz for the pc channels complete!')


    # --------------------------------------------------
    # Dead time correction on photon counting channels 
    # --------------------------------------------------
    if idtc:
        sig = signal.dead_time_correction(sig = sig.copy(), 
                                          dead_time = dead_time, 
                                          dead_time_cor_type = dead_time_cor_type)
        sig_dtc = sig.copy()
        
        print('-- Dead time correction for pc channels complete!')
    else:
        print('-- Warning: Dead time correction for pc channels deactivated!')


    # --------------------------------------------------
    # Temporal averaging 
    # --------------------------------------------------   
    if idtc:     
        sig, frames = signal.average_by_time(sig = sig.copy(), 
                                             timescale = timescale)
        sig_avg = sig.copy()
        
        print(f'-- Temporal averaging complete! [..{sig_avg.time.size} averaged timeframe(s)]')            
   
    else:
        print('-- Warning: Temporal averaging deactivated! The processing might delay considerably!')


    # --------------------------------------------------
    # Solar backsground signal calculation
    # --------------------------------------------------
    bgr = signal.background(sig = sig.copy(), 
                            lower_bin = bg_low,
                            upper_bin = bg_high)
    
    print('-- Solar background succesfully calculated!')

    # --------------------------------------------------
    # Correct trigger delay bins (pre-trigger or delay)
    # --------------------------------------------------
    sig = signal.trigger_delay(sig = sig.copy(), 
                               trigger_delay_bins = trd_bins)
    
    sig_trc = sig.copy()
    
    print('-- Triggering correction complete!')

    # --------------------------------------------------
    # Trim signal up to an upper altitude limit
    # --------------------------------------------------
    sig = signal.trim_vertically(sig = sig.copy(), 
                                 ground_alt = ground_alt,
                                 zenith_angle = zenith_angle, 
                                 alt_lim = alt_lim,
                                 resol = resol)
    
    sig_trm = sig.copy()
    
    print(f'-- Range bins sucessfully trimmed above {alt_lim}m distance!')

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


    if external_info['debug']:
        pack_out = {'sig_raw' : sig_raw,
                    'sig_puc' : sig_puc,
                    'sig_dtc' : sig_dtc,
                    'sig_avg' : sig_avg,
                    'sig_trc' : sig_trc,
                    'sig_trm' : sig_trm,
                    'bgr' : bgr,
                    'frames' : frames}
    else:
        pack_out = []
    
    return(pack_out)


#     # Background subtraction
#         sig = signal.background_correction(sig, bc)
#         sig_bc = sig.copy()
#         print('-- Sky background subtraction complete!')
            

#     # Calculate and add the reference altitude/range index and the window bins around it
#         info_sig = height_scales.bin_conversions(info_sig, sig_interp, cfg)
        
#     #Remove dark structure from normal signal
#         if not isinstance(dark,list):
#             dark_bc = dark/np.power(dark.range.values, 2) # convert the smoothed rc drk sig to smoothed bc drk sig
#             sig = signal.dark_correction(sig.copy(), info_sig, dark_bc)
#             sig_dc = sig.copy()
#             print('-- Dark signal structure correction complete!')
#         else:
#             sig_dc = sig.copy()
            
#     # Range correction
#         sig = signal.range_correction(sig)
#         sig_rc = sig.copy()
#         print('-- Range correction complete!')
    

#     # Smoothing  
#         if ifltr == 1:
#             if isdark:
#                 lims = drk_lims
#                 wins = drk_win
    
#             if not isdark:                 
#                 lims = sig_lims
#                 wins = sig_win
            
#             # Assign the end_of_signal range in the last range lim for smoothing 
#             if len(lims) > 0 and sig.range.values[-1] > lims[-1]:
#                 lims[-1] = sig.range.values[-1]
                
#             # Calculate half smoothing window in bins (always an even integer)
#             ihwin = smooth.get_win(x = sig.range.values, lims = lims, 
#                                    wins = wins, step = resol)
                        
#             # Apply smoothing
#             # sig.copy().sel(iters=slice(0)) # to keep the 4-D dim
#             sig = smooth.smoothing(sig.copy(), ihwin, fltr_ord)
#             sig_fltr = sig.copy()
#             print('-- Signals smoothing complete!')
        
#         else:
#             sig_fltr = sig.copy()

        
#     else: # empty signal. All the processed xarrays are empty
#         sig_mc, sig_dt, sig_avg, sig_bc, sig_bc_prt = [],[],[],[],[] 
#         sig_interp, sig_bc_trm, sig_fltr, sig_rc = [],[],[],[]
#         sig_dc, sig_gl, sig_tot = [],[],[]
        

# # Group all signals in lists and name lists
#     sig_pack = {}
    
#     sigs = [sig_dt, sig_avg, sig_mc, sig_bc, sig_bc_prt, 
#             sig_interp, sig_bc_trm, sig_dc, sig_gl, sig_rc, 
#             sig_tot, sig_fltr]
    
#     keys = ['sig_dt', 'sig_avg', 'sig_mc', 'sig_bc', 'sig_bc_prt', 
#             'sig_interp', 'sig_bc_trm',  'sig_dc', 'sig_gl', 'sig_rc',
#             'sig_tot', 'sig_fltr']    
    
#     for i in range(len(keys)):
#         sig_pack[keys[i]] = sigs[i]
                
#     return(sig_pack, info_sig, timescale)   

# def select_timescale(cfg, isdark):
#     if isdark:
#         timescale = cfg.timescale_d
#     else:
#         timescale = cfg.timescale
#     return(timescale)