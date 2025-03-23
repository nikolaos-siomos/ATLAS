#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob
from helper_functions import cleaner
from helper_functions import parse_master_args
from scc_converter.__scc_converter__ import main as __scc_converter__
from processor.__preprocessor__ import main as __preprocessor__
from visualizer.__quicklook__ import main as __quicklook__
from visualizer.__rayleigh_fit__ import main as __rayleigh_fit__
from visualizer.__telecover_quadrants__ import main as __telecover_quadrants__
from visualizer.__telecover_rings__ import main as __telecover_rings__
from visualizer.__polarization_calibration__ import main as __polarization_calibration__
from visualizer.__intercomparison__ import main as __intercomparison__
from scc_converter.readers.parse_args import check_parser as check_cnv
from processor.readers.parse_args import check_parser as check_prs
from visualizer.readers.parse_qck_args import check_parser as check_qck
from visualizer.readers.parse_ray_args import check_parser as check_ray
from visualizer.readers.parse_tlc_args import check_parser as check_tlc
from visualizer.readers.parse_pcb_args import check_parser as check_pcb
from visualizer.readers.parse_drk_args import check_parser as check_drk
from visualizer.readers.parse_cmp_args import check_parser as check_cmp
from scc_converter.readers.parse_args import view_parser as view_cnv
from processor.readers.parse_args import view_parser as view_prs
from visualizer.readers.parse_qck_args import view_parser as view_qck
from visualizer.readers.parse_ray_args import view_parser as view_ray
from visualizer.readers.parse_tlc_args import view_parser as view_tlc
from visualizer.readers.parse_pcb_args import view_parser as view_pcb
from visualizer.readers.parse_drk_args import view_parser as view_drk
from visualizer.readers.parse_cmp_args import view_parser as view_cmp
from visualizer.readers.parse_tlc_args import call_parser as parse_tlc
from scc_converter.readers.parse_args import call_parser as parse_cnv
from processor.readers.parse_args import call_parser as parse_prs
from visualizer.readers.parse_qck_args import call_parser as parse_qck
from visualizer.readers.parse_ray_args import call_parser as parse_ray
from visualizer.readers.parse_pcb_args import call_parser as parse_pcb
from visualizer.readers.parse_drk_args import call_parser as parse_drk

# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

def converter(mst_args, mst_cfg, cnv_out, processing, reprocess = True):

    cnv_ray_file = []
    cnv_tlc_file = []
    cnv_pcb_file = []
    cnv_drk_file = []
    
    if (processing['ray'] == True) or (processing['tlc'] == True) or \
        (processing['pcb'] == True) or (processing['drk'] == True):
        # Converter arguments
        cnv_args = parse_cnv()
    
        cnv_args = parse_master_args.substitute(org = cnv_args, rpl = mst_args)
        cnv_args = parse_master_args.substitute(org = cnv_args, rpl = mst_cfg.cnv)
        
    if processing['ray'] == True:
        
        cnv_ray_args = cnv_args.copy()   
        cnv_ray_args['mode'] = 'R'        
        cnv_ray_args['output_folder'] = cnv_out 
        cnv_ray_args = check_cnv(cnv_ray_args)
    
        os.makedirs(cnv_out, exist_ok = True)
    
        if reprocess == True:
            for file in glob.glob(os.path.join(cnv_out, '*_ray_ATLAS_*.nc')):
                os.remove(file)
                
        # QA files
        if os.path.exists(cnv_ray_args['rayleigh_folder']):                    
            cnv_ray_file = glob.glob(os.path.join(cnv_out, '*_ray_ATLAS_*.nc'))
            if len(cnv_ray_file) > 1:
                raise Exception(f'More than one rayleigh fit files detected in folder {cnv_out}. Please make sure that only one rayleigh file exists in that folder ')
           
            # Ececute scc_converter
            if len(cnv_ray_file) == 0 and processing['ray'] == True:
                view_cnv(cnv_ray_args)
                files = __scc_converter__(cnv_ray_args, __version__ = __version__)
                if files['rayleigh'] != None:
                    cnv_ray_file = [files['rayleigh']]
                else:
                    cnv_ray_file = []
        else:
            print(f"-- Warning: The {cnv_ray_args['rayleigh_folder']} folder does not exist. Please make sure to provide it if you want to process a Rayleigh fit test")
    else:
        # QA files
        cnv_ray_file = glob.glob(os.path.join(cnv_out, '*_ray_ATLAS_*.nc'))
        if len(cnv_ray_file) > 1:
            raise Exception(f'More than one rayleigh fit files detected in folder {cnv_out}. Please make sure that only one rayleigh file exists in that folder ')
       
        
    if processing['tlc'] == True:
        cnv_tlc_args = cnv_args.copy()  
        cnv_tlc_args['mode'] = 'T'  
        cnv_tlc_args['output_folder'] = cnv_out 
        cnv_tlc_args = check_cnv(cnv_tlc_args)
        
        os.makedirs(cnv_out, exist_ok = True)
    
        if reprocess == True:
            for file in glob.glob(os.path.join(cnv_out, '*_tlc_ATLAS_*.nc')):
                os.remove(file)
                
        # QA files  
        if os.path.exists(cnv_tlc_args['telecover_sectors_folder']) or \
                os.path.exists(cnv_tlc_args['telecover_rings_folder']):
                    
            cnv_tlc_file = glob.glob(os.path.join(cnv_out, '*_tlc_ATLAS_*.nc'))
            if len(cnv_tlc_file) > 1:
                raise Exception(f'More than one telecover files detected in folder {cnv_out}. Please make sure that only one telecover file exists in that folder ')
            
            # Ececute scc_converter
            if len(cnv_tlc_file) == 0 and processing['tlc'] == True:
                view_cnv(cnv_tlc_args)
                files = __scc_converter__(cnv_tlc_args, __version__ = __version__)
                if files['telecover'] != None:
                    cnv_tlc_file = [files['telecover']]                
                else:
                    cnv_ray_file = []
        else:
            print(f"-- Warning: Neither the {cnv_ray_args['telecover_sectors_folder']} nor the {cnv_ray_args['telecover_rings_folder']} folder does not exist. Please make sure to provide at least one of them if you want to process a telecover test")
        
    else:
        # QA files          
        cnv_tlc_file = glob.glob(os.path.join(cnv_out, '*_tlc_ATLAS_*.nc'))
        if len(cnv_tlc_file) > 1:
            raise Exception(f'More than one telecover files detected in folder {cnv_out}. Please make sure that only one telecover file exists in that folder ')


    if processing['pcb'] == True:

        cnv_pcb_args = cnv_args.copy()   
        cnv_pcb_args['mode'] = 'C' 
        cnv_pcb_args['output_folder'] = cnv_out 
        cnv_pcb_args = check_cnv(cnv_pcb_args)
        
        os.makedirs(cnv_out, exist_ok = True)
    
        if reprocess == True:
            for file in glob.glob(os.path.join(cnv_out, '*_pcb_ATLAS_*.nc')):
                os.remove(file)
                
        # QA files   
        if (os.path.exists(cnv_pcb_args['pcb_cal_p45_folder']) and \
            os.path.exists(cnv_pcb_args['pcb_cal_m45_folder'])) or \
                os.path.exists(cnv_pcb_args['pcb_cal_stc_folder']):
                    
            cnv_pcb_file = glob.glob(os.path.join(cnv_out, '*_pcb_ATLAS_*.nc'))        
            if len(cnv_pcb_file) > 1:
                raise Exception(f'More than one polarization calibration files detected in folder {cnv_out}. Please make sure that only one polarization calibration file exists in that folder ')         
    
            cnv_ray_file = glob.glob(os.path.join(cnv_out, '*_ray_ATLAS_*.nc'))
            if len(cnv_ray_file) > 1:
                raise Exception(f'More than one rayleigh fit files detected in folder {cnv_out}. Please make sure that only one polarization calibration file exists in that folder ')         
        
            if cnv_pcb_args['rayleigh_filename'] == None and len(cnv_ray_file) == 1:
                cnv_pcb_args['rayleigh_filename'] = os.path.basename(cnv_ray_file[0])
            else:
                raise Exception(f'Rayleigh file not found inside {cnv_out}. A Rayleigh measurement is mandatory for the pol. calibration test. Consider reprocessing with --quick_run deactivated ')
        
            # Ececute scc_converter
            if len(cnv_pcb_file) == 0 and processing['pcb'] == True:
                view_cnv(cnv_pcb_args)
                files = __scc_converter__(cnv_pcb_args, __version__ = __version__)
                if files['polarization_calibration'] != None:
                    cnv_pcb_file = [files['polarization_calibration']]
                else:
                    cnv_pcb_file = []
        else:
            print(f"-- Warning: None of the {cnv_ray_args['pcb_cal_p45_folder']}, {cnv_ray_args['pcb_cal_m45_folder']}, or {cnv_ray_args['pcb_cal_stc_folder']} was provided. Please make sure to provide at least the p45 and m45 or just the sta folder if you want to process a polarization calibration test")
            
    else:
        # QA files         
        cnv_pcb_file = glob.glob(os.path.join(cnv_out, '*_pcb_ATLAS_*.nc'))        
        if len(cnv_pcb_file) > 1:
            raise Exception(f'More than one polarization calibration files detected in folder {cnv_out}. Please make sure that only one polarization calibration file exists in that folder ')         

    if processing['drk'] == True:
        cnv_drk_args = cnv_args.copy()   
        cnv_drk_args['mode'] = 'D'        
        cnv_drk_args['output_folder'] = cnv_out 
        cnv_drk_args = check_cnv(cnv_drk_args)
    
        os.makedirs(cnv_out, exist_ok = True)
    
        if reprocess == True:
            for file in glob.glob(os.path.join(cnv_out, '*_drk_ATLAS_*.nc')):
                os.remove(file)
                
        # QA files
        if os.path.exists(cnv_drk_args['dark_folder']):                    
            cnv_drk_file = glob.glob(os.path.join(cnv_out, '*_drk_ATLAS_*.nc'))
            if len(cnv_drk_file) > 1:
                raise Exception(f'More than one rayleigh fit files detected in folder {cnv_out}. Please make sure that only one rayleigh file exists in that folder ')
           
            # Ececute scc_converter
            if len(cnv_drk_file) == 0 and processing['drk'] == True:
                view_cnv(cnv_drk_args)
                files = __scc_converter__(cnv_drk_args, __version__ = __version__)
                if files['dark'] != None:
                    cnv_drk_file = [files['dark']]
                else:
                    cnv_drk_file = []
        else:
            print(f"-- Warning: The {cnv_drk_args['dark_folder']} folder does not exist. Please make sure to provide it if you want to process a Dark test")
    else:
        # QA files
        cnv_drk_file = glob.glob(os.path.join(cnv_out, '*_drk_ATLAS_*.nc'))
        if len(cnv_drk_file) > 1:
            raise Exception(f'More than one dark files detected in folder {cnv_out}. Please make sure that only one dark exists in that folder ')
       

    files_out = {'ray' : cnv_ray_file, 
                 'tlc' : cnv_tlc_file, 
                 'pcb' : cnv_pcb_file,
                 'drk' : cnv_drk_file}
    
    return(files_out)

def preprocessor(mst_args, mst_cfg, cnv_files, prs_out, quicklook, processing, isday, reprocess = True):

    prs_ray_file = []
    prs_tlc_file = []
    prs_pcb_file = []
    prs_drk_file = []

    prs_qck_ray_file = []
    prs_qck_tlc_file = []
    prs_qck_pcb_file = []
    prs_qck_drk_file = []

    # Preprocessor arguments
    if (processing['ray'] == True and len(cnv_files['ray']) == 1) or \
        (processing['tlc'] == True and len(cnv_files['tlc']) == 1) or \
            (processing['pcb'] == True and len(cnv_files['pcb']) == 1) or \
                (processing['drk'] == True and len(cnv_files['drk']) == 1):
        prs_args = parse_prs()
    
        prs_args = parse_master_args.substitute(org = prs_args, rpl = mst_cfg.prs)
        
    if processing['ray'] == True and len(cnv_files['ray']) == 1:
        
        if isday == True:
            if prs_args['exclude_channel_type'] == None:
                prs_args['exclude_channel_type'] = ['v', 'f']
            else:
                prs_args['exclude_channel_type'] = \
                    prs_args['exclude_channel_type'] + ['v', 'f']
        
        prs_ray_args = prs_args.copy()  
        prs_ray_args['input_file'] = cnv_files['ray'][0]   
        prs_ray_args['quicklook'] = quicklook['ray']   
        prs_ray_args['output_folder'] = prs_out
        prs_ray_args = check_prs(prs_ray_args)
        
        os.makedirs(prs_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(prs_out, '*_ray_ATLAS_*.nc')):
                os.remove(file)
            for file in glob.glob(os.path.join(prs_out, '*_qck_ray_ATLAS_*.nc')):
                os.remove(file)
                
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_ray_file = glob.glob(os.path.join(prs_out, '*_ray_ATLAS_*_prepro.nc'))
        prs_ray_file = [item for item in prs_ray_file if 'qck_ray' not in item]

        if len(prs_ray_file) > 1:
            raise Exception(f'More than one rayleigh fit files detected in folder {prs_out}. Please make sure that only one rayleigh file exists in that folder ')
       
        if len(prs_ray_file) > 0:
            meas_ID = os.path.basename(prs_ray_file[0])[:15]
            prs_qck_ray_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_ray_ATLAS_*_prepro.nc'))
       
        # Excecute ATLAS preprocessor
        if len(prs_ray_file) == 0 and processing['ray'] == True:
            view_prs(prs_ray_args)
            files = __preprocessor__(prs_ray_args, __version__)
            prs_ray_file = files['ray']
            prs_qck_ray_file = files['qck'] 

    else:
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_ray_file = glob.glob(os.path.join(prs_out, '*_ray_ATLAS_*_prepro.nc'))
        prs_ray_file = [item for item in prs_ray_file if 'qck_ray' not in item]

        if len(prs_ray_file) > 1:
            raise Exception(f'More than one rayleigh fit files detected in folder {prs_out}. Please make sure that only one rayleigh file exists in that folder ')
       
        if len(prs_ray_file) > 0:
            meas_ID = os.path.basename(prs_ray_file[0])[:15]
            prs_qck_ray_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_ray_ATLAS_*_prepro.nc'))
   
    if processing['tlc'] == True and len(cnv_files['tlc']) == 1:
        
        prs_tlc_args = prs_args.copy()   
        prs_tlc_args['input_file'] = cnv_files['tlc'][0]        
        prs_tlc_args['quicklook'] = quicklook['tlc'] 
        prs_tlc_args['output_folder'] = prs_out
        prs_tlc_args = check_prs(prs_tlc_args)
    
        os.makedirs(prs_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(prs_out, '*_tlc_ATLAS_*.nc')):
                os.remove(file)
            for file in glob.glob(os.path.join(prs_out, '*_qck_tlc_ATLAS_*.nc')):
                os.remove(file)
                
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_tlc_file = glob.glob(os.path.join(prs_out, '*_tlc_ATLAS_*_prepro.nc'))
        prs_tlc_file = [item for item in prs_tlc_file if 'qck_tlc' not in item]

        if len(prs_tlc_file) > 1:
            raise Exception(f'More than one telecover files detected in folder {prs_out}. Please make sure that only one telecover file exists in that folder ')
    
        # Excecute ATLAS preprocessor
        if len(prs_tlc_file) > 0:
            meas_ID = os.path.basename(prs_tlc_file[0])[:15]
            prs_qck_tlc_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_tlc_ATLAS_*_prepro.nc'))
    
        if len(prs_tlc_file) == 0 and processing['tlc'] == True:
            view_prs(prs_tlc_args)
            files = __preprocessor__(prs_tlc_args, __version__)
            prs_tlc_file = files['tlc']
            prs_qck_tlc_file = files['qck']

    else:
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_tlc_file = glob.glob(os.path.join(prs_out, '*_tlc_ATLAS_*_prepro.nc'))
        prs_tlc_file = [item for item in prs_tlc_file if 'qck_tlc' not in item]

        if len(prs_tlc_file) > 1:
            raise Exception(f'More than one telecover files detected in folder {prs_out}. Please make sure that only one telecover file exists in that folder ')
    
        # Excecute ATLAS preprocessor
        if len(prs_tlc_file) > 0:
            meas_ID = os.path.basename(prs_tlc_file[0])[:15]
            prs_qck_tlc_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_tlc_ATLAS_*_prepro.nc'))
        else: prs_qck_tlc_file = []

    if processing['pcb'] == True and len(cnv_files['pcb']) == 1:

        prs_pcb_args = prs_args.copy()   
        prs_pcb_args['input_file'] = cnv_files['pcb'][0]   
        prs_pcb_args['quicklook'] = quicklook['pcb'] 
        prs_pcb_args['output_folder'] = prs_out
        prs_pcb_args = check_prs(prs_pcb_args)
        
        os.makedirs(prs_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(prs_out, '*_pcb_ATLAS_*.nc')):
                os.remove(file)
            for file in glob.glob(os.path.join(prs_out, '*_qck_pcb_ATLAS_*.nc')):
                os.remove(file)
            
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_pcb_file = glob.glob(os.path.join(prs_out, '*_pcb_ATLAS_*_prepro.nc'))
        prs_pcb_file = [item for item in prs_pcb_file if 'qck_pcb' not in item]

        if len(prs_pcb_file) > 1:
            raise Exception(f'More than one polarization calibration files detected in folder {prs_out}. Please make sure that only one polarization calibration file exists in that folder ')
        
        # Excecute ATLAS preprocessor            
        if len(prs_pcb_file) > 0:
            meas_ID = os.path.basename(prs_pcb_file[0])[:15]
            prs_qck_pcb_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_pcb_ATLAS_*_prepro.nc'))
    
        if len(prs_pcb_file) == 0 and processing['pcb'] == True:
            view_prs(prs_pcb_args)
            files = __preprocessor__(prs_pcb_args, __version__)
            prs_pcb_file = files['pcb']
            prs_qck_pcb_file = files['qck']

    else:
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_pcb_file = glob.glob(os.path.join(prs_out, '*_pcb_ATLAS_*_prepro.nc'))
        prs_pcb_file = [item for item in prs_pcb_file if 'qck_pcb' not in item]
   
        if len(prs_pcb_file) > 1:
            raise Exception(f'More than one polarization calibration files detected in folder {prs_out}. Please make sure that only one polarization calibration file exists in that folder ')
        
        # Excecute ATLAS preprocessor            
        if len(prs_pcb_file) > 0:
            meas_ID = os.path.basename(prs_pcb_file[0])[:15]
            prs_qck_pcb_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_pcb_ATLAS_*_prepro.nc'))
            
    if processing['drk'] == True and len(cnv_files['drk']) == 1:
        
        if isday == True:
            if prs_args['exclude_channel_type'] == None:
                prs_args['exclude_channel_type'] = ['v', 'f']
            else:
                prs_args['exclude_channel_type'] = \
                    prs_args['exclude_channel_type'] + ['v', 'f']
        
        prs_drk_args = prs_args.copy()  
        prs_drk_args['input_file'] = cnv_files['drk'][0]   
        prs_drk_args['quicklook'] = quicklook['drk']   
        prs_drk_args['output_folder'] = prs_out
        prs_drk_args = check_prs(prs_drk_args)
        
        os.makedirs(prs_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(prs_out, '*_drk_ATLAS_*.nc')):
                os.remove(file)
            for file in glob.glob(os.path.join(prs_out, '*_qck_drk_ATLAS_*.nc')):
                os.remove(file)
                
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_drk_file = glob.glob(os.path.join(prs_out, '*_drk_ATLAS_*_prepro.nc'))
        prs_drk_file = [item for item in prs_drk_file if 'qck_drk' not in item]

        if len(prs_drk_file) > 1:
            raise Exception(f'More than one dark files detected in folder {prs_out}. Please make sure that only one drk exists in that folder ')
       
        if len(prs_drk_file) > 0:
            meas_ID = os.path.basename(prs_drk_file[0])[:15]
            prs_qck_drk_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_drk_ATLAS_*_prepro.nc'))
       
        # Excecute ATLAS preprocessor
        if len(prs_drk_file) == 0 and processing['drk'] == True:
            view_prs(prs_drk_args)
            files = __preprocessor__(prs_drk_args, __version__)
            prs_drk_file = files['drk']
            prs_qck_drk_file = files['qck'] 

    else:
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_drk_file = glob.glob(os.path.join(prs_out, '*_drk_ATLAS_*_prepro.nc'))
        prs_drk_file = [item for item in prs_drk_file if 'qck_drk' not in item]

        if len(prs_drk_file) > 1:
            raise Exception(f'More than one dark files detected in folder {prs_out}. Please make sure that only one dark exists in that folder ')
       
        if len(prs_drk_file) > 0:
            meas_ID = os.path.basename(prs_drk_file[0])[:15]
            prs_qck_drk_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_qck_drk_ATLAS_*_prepro.nc'))
   
    files_out = {'ray' : prs_ray_file, 
                 'tlc' : prs_tlc_file, 
                 'pcb' : prs_pcb_file,
                 'drk' : prs_drk_file,
                 'qck_ray' : prs_qck_ray_file, 
                 'qck_tlc' : prs_qck_tlc_file, 
                 'qck_pcb' : prs_qck_pcb_file,
                 'qck_drk' : prs_qck_drk_file}
    
    return(files_out)

def ray_test(mst_cfg, prs_files, ray_out, processing, reprocess = True):

    if processing['ray'] == True and len(prs_files['ray']) == 1:
    
        # Rayleigh fit arguments
        ray_args = parse_ray()
        ray_args = parse_master_args.substitute(org = ray_args, rpl = mst_cfg.ray)

        ray_args = ray_args.copy()   
        
        ray_args['input_file'] = prs_files['ray'][0]   
        
        ray_args = check_ray(ray_args)
        
        ray_args['output_folder'] = ray_out
        os.makedirs(ray_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(ray_out, 'plots', '*_ray_*_ATLAS_*.png')):
                os.remove(file)
            for file in glob.glob(os.path.join(ray_out, 'ascii', '*_ray_*_ATLAS_*.txt')):
                os.remove(file)    
                
        # Excecute ATLAS visualizer
        view_ray(ray_args)
        __rayleigh_fit__(ray_args, __version__)

    return()

def intercomparison(cmp_args, input_files, cmp_out):
    
    cmp_args['input_files'] = input_files 
    cmp_args['output_folder'] = cmp_out
        
    os.makedirs(cmp_out, exist_ok = True)
    
    cmp_args = cmp_args.copy()   
    
    cmp_args = check_cmp(cmp_args)

    # Excecute ATLAS visualizer
    view_cmp(cmp_args)
    __intercomparison__(cmp_args, __version__)

    return()

def tlc_test(mst_cfg, prs_files, tlc_out, processing, reprocess = True):
    
    if processing['tlc'] == True  and len(prs_files['tlc']) == 1:

        # Telecover arguments
        tlc_args = parse_tlc()
        tlc_args = parse_master_args.substitute(org = tlc_args, rpl = mst_cfg.tlc)
                
        tlc_args = tlc_args.copy()   
        
        tlc_args['input_file'] = prs_files['tlc'][0]   
        
        tlc_args = check_tlc(tlc_args)
        
        tlc_args['output_folder'] = tlc_out
        os.makedirs(tlc_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(tlc_out, 'plots', '*_tlc_*_ATLAS_*.png')):
                os.remove(file)
            for file in glob.glob(os.path.join(tlc_out, 'ascii', '*_tlc_*_ATLAS_*.txt')):
                os.remove(file)    
       
        # Excecute ATLAS visualizer
        view_tlc(tlc_args)
        __telecover_quadrants__(tlc_args, __version__)

        __telecover_rings__(tlc_args, __version__)

    return()

def pcb_test(mst_cfg, prs_files, pcb_out, processing, reprocess = True):
    
    if processing['pcb'] == True and len(prs_files['pcb']) == 1:

        # Polarization Calibration arguments
        pcb_args = parse_pcb()
        pcb_args = parse_master_args.substitute(org = pcb_args, rpl = mst_cfg.pcb)

        pcb_args = pcb_args.copy()   
        
        pcb_args['input_file'] = prs_files['pcb'][0]   
        
        pcb_args = check_pcb(pcb_args)
        
        pcb_args['output_folder'] = pcb_out
        os.makedirs(pcb_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(pcb_out, 'plots', '*_pcb_*_ATLAS_*.png')):
                os.remove(file)
            for file in glob.glob(os.path.join(pcb_out, 'ascii', '*_pcb_*_ATLAS_*.txt')):
                os.remove(file)    
    
        # Excecute ATLAS visualizer
        view_pcb(pcb_args)
        __polarization_calibration__(pcb_args, __version__)

    return()

def drk_test(mst_cfg, prs_files, drk_out, processing, reprocess = True):
    
    if processing['drk'] == True and len(prs_files['drk']) == 1:

        # Polarization Calibration arguments
        drk_args = parse_drk()
        drk_args = parse_master_args.substitute(org = drk_args, rpl = mst_cfg.drk)

        drk_args = drk_args.copy()   
        
        drk_args['input_file'] = prs_files['drk'][0]   
        
        drk_args = check_drk(drk_args)
        
        drk_args['output_folder'] = drk_out
        os.makedirs(drk_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(drk_out, 'plots', '*_drk_*_ATLAS_*.png')):
                os.remove(file)
            for file in glob.glob(os.path.join(drk_out, 'ascii', '*_drk_*_ATLAS_*.txt')):
                os.remove(file)    
    
        # Excecute ATLAS visualizer
        # view_drk(drk_args)
        # __dark__(drk_args, __version__)

    return()
def quicklook(mst_cfg, prs_files, qck_out, quicklook, reprocess = True):

    if quicklook['ray'] == True and len(prs_files['qck_ray']) == 1:

        # Quicklook arguments
        qck_args = parse_qck()
        qck_args = parse_master_args.substitute(org = qck_args, rpl = mst_cfg.qck)

        qck_ray_args = qck_args.copy()   
        
        qck_ray_args['input_file'] = prs_files['qck_ray'][0]   
        
        qck_ray_args = check_qck(qck_ray_args)
       
        qck_ray_args['output_folder'] = qck_out
        os.makedirs(qck_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(qck_out, 'plots', '*_qck_ray_*_ATLAS_*.png')):
                os.remove(file)
    
    
        # Excecute ATLAS visualizer
        view_qck(qck_ray_args)
        __quicklook__(qck_ray_args, __version__, meas_type = 'ray')
    
    elif quicklook['ray'] == True and len(prs_files['qck_ray']) == 0:
        print('-- Warning: Rayleight quicklook visualization was enabled but no file from the preprocessing stage was found. Please make sure that --process_qck ray is at least provided. You might also need to deactivate --quick_run')
        
    if quicklook['tlc'] == True and len(prs_files['qck_tlc']) == 1:
        
        qck_args = parse_qck()
        qck_args = parse_master_args.substitute(org = qck_args, rpl = mst_cfg.qck)

        qck_tlc_args = qck_args.copy() 

        if len(prs_files['qck_tlc']) == 1 and len(prs_files['qck_tlc']) == 1:
            qck_tlc_args['input_file'] = prs_files['qck_tlc'][0]        
        
        qck_tlc_args = check_qck(qck_tlc_args)
    
        qck_tlc_args['output_folder'] = qck_out
        os.makedirs(qck_out, exist_ok = True)
    
        if reprocess == True:
            for file in glob.glob(os.path.join(qck_out, 'plots', '*_qck_tlc_*_ATLAS_*.nc')):
                os.remove(file)
        
        # Excecute ATLAS visualizer
        view_qck(qck_tlc_args)
        __quicklook__(qck_tlc_args, __version__, meas_type = 'tlc')

    elif quicklook['tlc'] == True and len(prs_files['qck_tlc']) == 0:
        print('-- Warning: Telecover quicklook visualization was enabled but no file from the preprocessing stage was found. Please make sure that --process_qck tlc is at least provided. You might also need to deactivate --quick_run')
        
    if quicklook['pcb'] == True and len(prs_files['qck_pcb']) == 1:
        
        qck_args = parse_qck()
        qck_args = parse_master_args.substitute(org = qck_args, rpl = mst_cfg.qck)

        qck_pcb_args = qck_args.copy()   

        if len(prs_files['qck_pcb']) == 1 and len(prs_files['qck_pcb']) == 1:
            qck_pcb_args['input_file'] = prs_files['qck_pcb'][0]   
        
        qck_pcb_args = check_qck(qck_pcb_args)    
        
        qck_pcb_args['output_folder'] = qck_out
        os.makedirs(qck_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(qck_out, 'plots', '*_qck_pcb_*_ATLAS_*.nc')):
                os.remove(file)
            
        # Excecute ATLAS visualizer
        view_qck(qck_pcb_args)
        __quicklook__(qck_pcb_args, __version__, meas_type = 'pcb')

    elif quicklook['pcb'] == True and len(prs_files['qck_pcb']) == 0:
        print('-- Warning: Polarization calibration quicklook visualization was enabled but no file from the preprocessing stage was found. Please make sure that --process_qck tlc is at least provided. You might also need to deactivate --quick_run')
            
    if quicklook['drk'] == True and len(prs_files['qck_drk']) == 1:

        # Quicklook arguments
        qck_args = parse_qck()
        qck_args = parse_master_args.substitute(org = qck_args, rpl = mst_cfg.qck)

        qck_drk_args = qck_args.copy()   
        
        qck_drk_args['input_file'] = prs_files['qck_drk'][0]   
        
        qck_drk_args = check_qck(qck_drk_args)
       
        qck_drk_args['output_folder'] = qck_out
        os.makedirs(qck_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(qck_out, 'plots', '*_qck_drk_*_ATLAS_*.png')):
                os.remove(file)
    
        # Excecute ATLAS visualizer
        view_qck(qck_drk_args)
        __quicklook__(qck_drk_args, __version__, meas_type = 'drk')
    
    elif quicklook['drk'] == True and len(prs_files['qck_drk']) == 0:
        print('-- Warning: Dark quicklook visualization was enabled but no file from the preprocessing stage was found. Please make sure that --process_qck ray is at least provided. You might also need to deactivate --quick_run')
                
    return()