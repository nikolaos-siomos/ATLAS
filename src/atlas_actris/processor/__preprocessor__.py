#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:48:11 2022

@author: nick
"""

import warnings, os, sys
from .readers.parse_args import call_parser, check_parser
from .readers import read_files
from .lidar_processing import short_prepro
from .arc import atmosphere
from .export import nc_dataset

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

def main(args, __version__):
              
    output_files = dict()
    
    print('-----------------------------------------')
    print('Initializing Preprocessor...')
    print('-----------------------------------------')
    print(' ')
    
    print('-----------------------------------------')
    print('Start reading the QA file(s)...')
    
    system_info, channel_info, time_info, time_info_d, \
        sig_raw, sig_raw_d, shots, shots_d = \
            read_files.short_reader(
                args['input_file'],
                exclude_telescope_type = args['exclude_telescope_type'], 
                exclude_channel_type = args['exclude_channel_type'], 
                exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                exclude_channel_subtype = args['exclude_channel_subtype'], 
                use_channels = args['channels'])
            

    print(f'Lidar: {system_info.Lidar_Name}')
    print('-----------------------------------------')
    print(' ')
    
    meas_type = system_info.Measurement_type
    meas_label = {'ray' : 'Rayleigh', 
                  'tlc' : 'Telecover', 
                  'pcb' : 'Polarization Calibration',
                  'drk' : 'Dark'}
    print(f'-- {meas_label[meas_type]} QA file succesfully parsed!')
    
    if 'Rayleigh_File_Name' in system_info.keys():
        ray_path = os.path.join(os.path.dirname(args['input_file']),  
                                system_info['Rayleigh_File_Name'])
        
        system_info_r, channel_info_r, time_info_r, time_info_dr, \
            sig_raw_r, sig_raw_dr, shots_r, shots_dr = \
                read_files.short_reader(
                    ray_path,
                    exclude_telescope_type = args['exclude_telescope_type'], 
                    exclude_channel_type = args['exclude_channel_type'], 
                    exclude_acquisition_mode = args['exclude_acquisition_mode'], 
                    exclude_channel_subtype = args['exclude_channel_subtype'], 
                    use_channels = args['channels'])
        
        if all([ch in channel_info_r.index.values 
                for ch in channel_info.index.values]):
            print(f"-- {meas_label['ray']} QA file used for calibration "+
                  "succesfully parsed!")
        else:
            raise Exception("-- Error! The Rayleigh measurement accompaning the " +
                            "calibration measurement {system_info['Rayleigh_File_Name']} " +
                            "lacks at least of the calibration measurement channels! " +
                            "Please make sure that the Rayleigh measuremen contains " +
                            "at least the channels channels necessary for the calibration")
    
    print('-----------------------------------------')
    print("")

    if meas_type != 'drk' and not isinstance(sig_raw_d, list):
        drk, drk_pack, time_info_d = \
            short_prepro.dark(sig_raw = sig_raw_d, 
                              shots = shots_d, 
                              system_info = system_info, 
                              channel_info = channel_info, 
                              time_info = time_info_d,
                              external_info = args)
        sig_drk = drk_pack['sig_bgc']
    else:
        sig_drk = []
    
    if meas_type != 'drk' and 'Rayleigh_File_Name' in system_info.keys() and not isinstance(sig_raw_d, list):
        drk_r, drk_pack_r, time_info_dr = \
            short_prepro.dark(sig_raw = sig_raw_dr, 
                              shots = shots_dr, 
                              system_info = system_info_r, 
                              channel_info = channel_info_r, 
                              time_info = time_info_dr,
                              external_info = args)
        sig_drk_r = drk_pack_r['sig_bgc']
    else:
        sig_drk_r = []
            
    if meas_type == 'ray':
                
        ray, ray_pack, time_info_ray = \
            short_prepro.standard(sig_raw = sig_raw, 
                                  shots = shots, 
                                  system_info = system_info, 
                                  channel_info = channel_info, 
                                  time_info = time_info,
                                  external_info = args,
                                  meas_type = meas_type,
                                  sig_drk = sig_drk)
        
        molec, molec_info, meteo = \
            atmosphere.short_molec(heights = ray_pack['heights'],
                                   ranges = ray_pack['ranges'],
                                   system_info = system_info, 
                                   channel_info = channel_info, 
                                   time_info = time_info_ray,
                                   external_info = args)
        
        output_files['ray'] = \
            nc_dataset.rayleigh(sig = ray, 
                                molec = molec,
                                meteo = meteo,
                                system_info = system_info, 
                                channel_info = channel_info, 
                                time_info = time_info_ray, 
                                molec_info = molec_info, 
                                heights = ray_pack['heights'], 
                                ranges = ray_pack['ranges'],
                                version = __version__,
                                dir_out = args['output_folder'])
    
    if meas_type == 'tlc':
        
        tlc, tlc_pack, time_info_tlc = \
            short_prepro.standard(sig_raw = sig_raw, 
                                  shots = shots, 
                                  system_info = system_info, 
                                  channel_info = channel_info, 
                                  time_info = time_info,
                                  external_info = args,
                                  meas_type = meas_type,
                                  sig_drk = sig_drk)
    
        output_files['tlc'] = \
            nc_dataset.telecover(sig = tlc, 
                                 system_info = system_info, 
                                 channel_info = channel_info, 
                                 time_info = time_info_tlc, 
                                 heights = tlc_pack['heights'], 
                                 ranges = tlc_pack['ranges'],
                                 version = __version__,
                                 dir_out = args['output_folder'])
    
    if meas_type == 'pcb':
                
        pcb, pcb_pack, time_info_pcb = \
            short_prepro.standard(sig_raw = sig_raw, 
                                  shots = shots, 
                                  system_info = system_info, 
                                  channel_info = channel_info, 
                                  time_info = time_info,
                                  external_info = args,
                                  meas_type = meas_type,
                                  sig_drk = sig_drk)
        
        ray, ray_pack, time_info_ray = \
            short_prepro.standard(sig_raw = sig_raw_r, 
                                  shots = shots_r, 
                                  system_info = system_info_r, 
                                  channel_info = channel_info_r, 
                                  time_info = time_info_r,
                                  external_info = args,
                                  meas_type = 'ray',
                                  sig_drk = sig_drk_r)
        
        molec, molec_info, meteo = \
            atmosphere.short_molec(heights = ray_pack['heights'],
                                   ranges = ray_pack['ranges'],
                                   system_info = system_info_r, 
                                   channel_info = channel_info_r, 
                                   time_info = time_info_ray,
                                   external_info = args)
            
        pcb_ch = pcb.channel.values
        
        output_files['pcb'] = \
            nc_dataset.calibration(sig = pcb, 
                                   sig_ray = ray.loc[dict(channel = pcb_ch)],
                                   molec = molec.loc[dict(channel = pcb_ch)],
                                   meteo = meteo.loc[dict(channel = pcb_ch)],
                                   system_info = system_info, 
                                   system_info_ray = system_info_r, 
                                   channel_info = channel_info, 
                                   channel_info_ray = channel_info.loc[pcb_ch], 
                                   time_info = time_info_pcb, 
                                   time_info_ray = time_info_ray, 
                                   molec_info = molec_info, 
                                   heights = pcb_pack['heights'], 
                                   ranges = pcb_pack['ranges'],
                                   ranges_ray = ray_pack['ranges'],
                                   heights_ray = ray_pack['heights'], 
                                   version = __version__,
                                   dir_out = args['output_folder'])
    
    if meas_type == 'drk':
        args['vertical_trimming'] = False
        drk, drk_pack, time_info_drk = \
            short_prepro.standard(sig_raw = sig_raw, 
                                  shots = shots, 
                                  system_info = system_info, 
                                  channel_info = channel_info, 
                                  time_info = time_info,
                                  external_info = args,
                                  meas_type = 'drk')
        sig_drk = []

        output_files['drk'] = \
            nc_dataset.dark(sig = drk, 
                            system_info = system_info, 
                            channel_info = channel_info, 
                            time_info = time_info, 
                            heights = drk_pack['heights'], 
                            ranges = drk_pack['ranges'],
                            version = __version__,
                            dir_out = args['output_folder'])

    if args['quicklook']:
                
        qck, qck_pack, time_info_qck = \
            short_prepro.standard(sig_raw = sig_raw, 
                                  shots = shots, 
                                  system_info = system_info, 
                                  channel_info = channel_info, 
                                  time_info = time_info,
                                  external_info = args,
                                  meas_type = 'qck',
                                  sig_drk = sig_drk)
            
        if meas_type == 'ray':
            output_files['qck'] = \
                nc_dataset.quicklook(sig = qck, 
                                     system_info = system_info, 
                                     channel_info = channel_info, 
                                     time_info = time_info_qck, 
                                     heights = qck_pack['heights'], 
                                     ranges = qck_pack['ranges'],
                                     version = __version__,
                                     meas_type = meas_type,
                                     dir_out = args['output_folder'],
                                     molec = molec,
                                     meteo = meteo,
                                     molec_info = molec_info)
        else:
            output_files['qck'] = \
                nc_dataset.quicklook(sig = qck, 
                                     system_info = system_info, 
                                     channel_info = channel_info, 
                                     time_info = time_info_qck, 
                                     heights = qck_pack['heights'], 
                                     ranges = qck_pack['ranges'],
                                     version = __version__,
                                     meas_type = meas_type,
                                     dir_out = args['output_folder'])
    
    else:
        output_files['qck'] = []

    return(output_files)
    
        
if __name__ == '__main__':
    sys.path.append('../')
    from version import __version__
    
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args, __version__)