"""
@authors: N. Siomos & P. Paschou
"""

from lidar_processing import signal
 
def short_prepro(sig, shots, lidar_info, channel_info, time_info, dark_pack = []):
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
        
    resol = channel_info.Raw_Data_Range_Resolution 
    dead_time = channel_info.Dead_Time
    dead_time_cor_type = channel_info.Dead_Time_Correction_Type
    
    
    # Unit conversion - raw counts to MHz for the photon channels
    sig = signal.unit_conv_counts_to_MHz(sig = sig.copy(), 
                                         shots = shots.copy(), 
                                         resol = resol)
    sig_cu = sig.copy()

    # Dead time correction on photon counting channels (only for normal signal)
    sig = signal.dead_time_correction(sig = sig.copy(), 
                                      dead_time = dead_time, 
                                      dead_time_cor_type = dead_time_cor_type)
    sig_dt = sig.copy()
    
    return(sig_cu, sig_dt)

#     # Time averaging
#         sig, frames, timescale = \
#             signal.average_by_time(sig.copy(), timescale)
#         sig_avg = sig.copy()
#         if timescale != '':
#             print(f'-- Time averaging complete! [..{len(sig_avg.time)} averaged timeframe(s)]')            
#         else:
#             timescale='0min'

#     # Simulate errors, expands to 4D arrays even if error sim is deactivated       
#         if not isinstance(dark_avg,list): # normal signal with dark signal available
#             sig, info_sig = error(sig.copy(), cfg, info_sig, frames, 
#                                   isdark, drk_inp = dark_avg.copy()) 

#         if isinstance(dark_avg,list): # signal without dark signal availability
#             sig, info_sig = error(sig.copy(), cfg, info_sig, frames, isdark) 

#         sig_mc = sig.copy()
        
#         mc_sim = 0        
#         if cfg.ierror > 0:
#             mc_sim = cfg.iters
#         print(f'-- Error simulations [{mc_sim}] and 4D array expansion complete!')

#     # Background subtraction
#         bc = signal.background(sig.copy(), prt_lims.start, prt_lims.end)
#         sig = signal.background_correction(sig, bc)
#         sig_bc = sig.copy()
#         print('-- Sky background subtraction complete!')
            
#     # Correct trigger delay bins (pre-trigger or delay)
#         sig = signal.trigger_delay(sig.copy(), info_sig, tr_bins)
#         sig_bc_prt = sig.copy()
#         print('-- Pre-trigger or trigger delay correction complete!')
        
#     # Change spatial resolution and convert to range  
#         sig, info_sig, resol = \
#             signal.change_resol(sig.copy(), info_sig, new_resol)
#         sig_interp = sig.copy()
#         print(f'-- Common spatial resolution ({resol} m) '+\
#               'and bins converted to range!')
    
#     # Calculate and add the reference altitude/range index and the window bins around it
#         info_sig = height_scales.bin_conversions(info_sig, sig_interp, cfg)
        
#     # Signal Trimming
#         if alt_lim != '':
#             sig, alt_lim = signal.trim(sig.copy(), ground = ground_alt,
#                                        SZA = SZA, alt_lim = alt_lim)
#             sig_bc_trm = sig.copy()
#             print(f'-- Signals trimed at {alt_lim} m altitude!')  
#         else:
#             sig_bc_trm = []
        
#     #Remove dark structure from normal signal
#         if not isinstance(dark,list):
#             dark_bc = dark/np.power(dark.range.values, 2) # convert the smoothed rc drk sig to smoothed bc drk sig
#             sig = signal.dark_correction(sig.copy(), info_sig, dark_bc)
#             sig_dc = sig.copy()
#             print('-- Dark signal structure correction complete!')
#         else:
#             sig_dc = sig.copy()
            
#     # Glue of the smoothed bc corrected signals
#         if not isdark and iglue == 1:
#             gl = read_maps.gluing(maps_dir, 'gluing_map.ini')
            
#             if len(gl.map)>0:
#                 sig, info_sig = glue.gluing(sig.copy(), info_sig, gl, bc)
#                 print('-- Gluing complete!')
            
#         sig_gl = sig.copy()

#     # Range correction
#         sig = signal.range_correction(sig)
#         sig_rc = sig.copy()
#         print('-- Range correction complete!')
    
#     # Construction of total signal
#         if not isdark and itot == 1:
#             tot = read_maps.total_sig(maps_dir, 'total_map.ini')
            
#             if len(tot.map)>0:
#                 sig, info_sig = signal.construct_total(sig.copy(), info_sig, tot)
#                 print('-- Calibrated sum signals constructed!')
            
#         sig_tot = sig.copy()    

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