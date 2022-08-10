"""
@author: N. Siomos & P. Paschou

Main routine for pre-processing the raw signals
"""
from lidar_processing.short_prepro import short_prepro

def process_signals(meas_type, sig_raw, info_val, cfg, maps_dir, dark_pack = []):
    # Pre processing
    print('-----------------------------------------')
    print(f'Start processing {meas_type} signal measurements...')
    print('-----------------------------------------')
    
    if meas_type == 'dark':
        isdark=True
    else:
        isdark=False
    
    # No measurements found
    if len(sig_raw) == 0:
        stime = ''
        sig_pack = {}
        info_sig=[]
        if isdark and cfg.idark == 1:
            print(f'-- Warning! No files found in folder {meas_type} --> \n '+\
                  '---> Skipping processing of dark measurements') 
        else:
            raise ValueError(f'No files found in folder {meas_type}! Program Stopped')
    else:
        stime = str(sig_raw.time.values[0])

    # signal pre-processing procedures
    #** the short prepro returns empty lists in case of empty sig_raw
    if isdark:
        if cfg.idark == 1:
            sig_pack, info_sig, timescale = short_prepro(sig_raw.copy(), info_val, 
                                                         maps_dir = maps_dir, cfg = cfg,
                                                         isdark = True)
        if cfg.idark == 0:
            sig_pack = []
            info_sig = []
            timescale = []

    if isdark == False:
        sig_pack, info_sig, timescale = short_prepro(sig_raw.copy(), info_val, 
                                                     maps_dir = maps_dir, cfg = cfg,
                                                     isdark = False, dark_pack = dark_pack) 
        print(f'Processing of {meas_type} signals complete!')
        
    print('-----------------------------------------')
    print('')

    return(sig_pack, info_sig, stime, timescale)