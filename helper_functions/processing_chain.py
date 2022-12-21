#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:21:48 2022

@author: nick
"""

from version import __version__
import os, warnings, glob
from helper_functions import cleaner
from scc_converter.__scc_converter__ import main as __scc_converter__
from processor.__preprocessor__ import main as __preprocessor__
from visualizer.__quicklook__ import main as __quicklook__
from visualizer.__rayleigh_fit__ import main as __rayleigh_fit__
from visualizer.__telecover__ import main as __telecover__
from visualizer.__polarization_calibration__ import main as __polarization_calibration__
from visualizer.__intercomparison__ import main as __intercomparison__
from scc_converter.readers.parse_args import check_parser as check_cnv
from processor.readers.parse_args import check_parser as check_prs
from visualizer.readers.parse_qck_args import check_parser as check_qck
from visualizer.readers.parse_ray_args import check_parser as check_ray
from visualizer.readers.parse_tlc_args import check_parser as check_tlc
from visualizer.readers.parse_pcb_args import check_parser as check_pcb
from visualizer.readers.parse_cmp_args import check_parser as check_cmp
from scc_converter.readers.parse_args import view_parser as view_cnv
from processor.readers.parse_args import view_parser as view_prs
from visualizer.readers.parse_qck_args import view_parser as view_qck
from visualizer.readers.parse_ray_args import view_parser as view_ray
from visualizer.readers.parse_tlc_args import view_parser as view_tlc
from visualizer.readers.parse_pcb_args import view_parser as view_pcb
from visualizer.readers.parse_cmp_args import view_parser as view_cmp

# from visualizer.readers.parse_cmp_args import call_parser as parse_cmp

warnings.filterwarnings('ignore')

def converter(cnv_args, cnv_out, processing, reprocess = True):

    cnv_ray_file = []
    cnv_tlc_file = []
    cnv_pcb_file = []
    
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
        cnv_pcb_file = glob.glob(os.path.join(cnv_out, '*_pcb_ATLAS_*.nc'))        
        if len(cnv_pcb_file) > 1:
            raise Exception(f'More than one polarization calibration files detected in folder {cnv_out}. Please make sure that only one polarization calibration file exists in that folder ')         

        cnv_ray_file = glob.glob(os.path.join(cnv_out, '*_ray_ATLAS_*.nc'))
        if len(cnv_ray_file) > 1:
            raise Exception(f'More than one rayleigh fit files detected in folder {cnv_out}. Please make sure that only one polarization calibration file exists in that folder ')         
    
        if cnv_pcb_args['rayleigh_filename'] == None and len(cnv_ray_file) == 1:
            cnv_pcb_args['rayleigh_filename'] = os.path.basename(cnv_ray_file[0])
        else:
            raise Exception(f'Rayleigh file not found inside {cnv_out}. A Rayleigh measurement is mandatory for the pol. calibration test. Include the Rayleigh files if not already included and reprocess. If the files are already there consider setting newdata = True ')
    
        # Ececute scc_converter
        if len(cnv_pcb_file) == 0 and processing['pcb'] == True:
            view_cnv(cnv_pcb_args)
            files = __scc_converter__(cnv_pcb_args, __version__ = __version__)
            if files['polarization_calibration'] != None:
                cnv_pcb_file = [files['polarization_calibration']]
            else:
                cnv_pcb_file = []
    else:
        # QA files          
        cnv_pcb_file = glob.glob(os.path.join(cnv_out, '*_pcb_ATLAS_*.nc'))        
        if len(cnv_pcb_file) > 1:
            raise Exception(f'More than one polarization calibration files detected in folder {cnv_out}. Please make sure that only one polarization calibration file exists in that folder ')         


    files_out = {'ray' : cnv_ray_file, 
                 'tlc' : cnv_tlc_file, 
                 'pcb' : cnv_pcb_file}
    
    return(files_out)

def preprocessor(prs_args, cnv_files, prs_out, quicklook, processing, reprocess = True):

    prs_ray_file = []
    prs_tlc_file = []
    prs_pcb_file = []

    prs_ray_qck_file = []
    prs_tlc_qck_file = []
    prs_pcb_qck_file = []
    
    if processing['ray'] == True and len(cnv_files['ray']) == 1:
        
        prs_ray_args = prs_args.copy()  
        prs_ray_args['input_file'] = cnv_files['ray'][0]   
        prs_ray_args['quicklook'] = quicklook['ray']   
        prs_ray_args['output_folder'] = prs_out
        prs_ray_args = check_prs(prs_ray_args)
        
        os.makedirs(prs_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(prs_out, '*_ray_ATLAS_*.nc')):
                os.remove(file)
            for file in glob.glob(os.path.join(prs_out, '*_ray_qck_ATLAS_*.nc')):
                os.remove(file)
                
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_ray_file = glob.glob(os.path.join(prs_out, '*_ray_ATLAS_*_prepro.nc'))
         
        if len(prs_ray_file) > 1:
            raise Exception(f'More than one rayleigh fit files detected in folder {prs_out}. Please make sure that only one rayleigh file exists in that folder ')
       
        if len(prs_ray_file) > 0:
            meas_ID = os.path.basename(prs_ray_file[0])[:15]
            prs_ray_qck_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_ray_qck_ATLAS_*_prepro.nc'))
       
        # Excecute ATLAS preprocessor
        if len(prs_ray_file) == 0 and processing['ray'] == True:
            view_prs(prs_ray_args)
            files = __preprocessor__(prs_ray_args, __version__)
            prs_ray_file = files['ray']
            prs_ray_qck_file = files['qck'] 

    else:
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_ray_file = glob.glob(os.path.join(prs_out, '*_ray_ATLAS_*_prepro.nc'))
 
        if len(prs_ray_file) > 1:
            raise Exception(f'More than one rayleigh fit files detected in folder {prs_out}. Please make sure that only one rayleigh file exists in that folder ')
       
        if len(prs_ray_file) > 0:
            meas_ID = os.path.basename(prs_ray_file[0])[:15]
            prs_ray_qck_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_ray_qck_ATLAS_*_prepro.nc'))
   
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
            for file in glob.glob(os.path.join(prs_out, '*_tlc_qck_ATLAS_*.nc')):
                os.remove(file)
                
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_tlc_file = glob.glob(os.path.join(prs_out, '*_tlc_ATLAS_*_prepro.nc'))
    
        if len(prs_tlc_file) > 1:
            raise Exception(f'More than one telecover files detected in folder {prs_out}. Please make sure that only one telecover file exists in that folder ')
    
        # Excecute ATLAS preprocessor
        if len(prs_tlc_file) > 0:
            meas_ID = os.path.basename(prs_tlc_file[0])[:15]
            prs_tlc_qck_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_tlc_qck_ATLAS_*_prepro.nc'))
    
        if len(prs_tlc_file) == 0 and processing['tlc'] == True:
            view_prs(prs_tlc_args)
            files = __preprocessor__(prs_tlc_args, __version__)
            prs_tlc_file = files['tlc']
            prs_tlc_qck_file = files['qck']

    else:
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_tlc_file = glob.glob(os.path.join(prs_out, '*_tlc_ATLAS_*_prepro.nc'))
    
        if len(prs_tlc_file) > 1:
            raise Exception(f'More than one telecover files detected in folder {prs_out}. Please make sure that only one telecover file exists in that folder ')
    
        # Excecute ATLAS preprocessor
        if len(prs_tlc_file) > 0:
            meas_ID = os.path.basename(prs_tlc_file[0])[:15]
            prs_tlc_qck_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_tlc_qck_ATLAS_*_prepro.nc'))
        else: prs_tlc_qck_file = []

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
            for file in glob.glob(os.path.join(prs_out, '*_pcb_qck_ATLAS_*.nc')):
                os.remove(file)
            
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_pcb_file = glob.glob(os.path.join(prs_out, '*_pcb_ATLAS_*_prepro.nc'))
        
        if len(prs_pcb_file) > 1:
            raise Exception(f'More than one polarization calibration files detected in folder {prs_out}. Please make sure that only one polarization calibration file exists in that folder ')
        
        # Excecute ATLAS preprocessor            
        if len(prs_pcb_file) > 0:
            meas_ID = os.path.basename(prs_pcb_file[0])[:15]
            prs_pcb_qck_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_pcb_qck_ATLAS_*_prepro.nc'))
    
        if len(prs_pcb_file) == 0 and processing['pcb'] == True:
            view_prs(prs_pcb_args)
            files = __preprocessor__(prs_pcb_args, __version__)
            prs_pcb_file = files['pcb']
            prs_pcb_qck_file = files['qck']

    else:
        # Preprocessed files - shoulde be zero if reprocess is True
        prs_pcb_file = glob.glob(os.path.join(prs_out, '*_pcb_ATLAS_*_prepro.nc'))
        
        if len(prs_pcb_file) > 1:
            raise Exception(f'More than one polarization calibration files detected in folder {prs_out}. Please make sure that only one polarization calibration file exists in that folder ')
        
        # Excecute ATLAS preprocessor            
        if len(prs_pcb_file) > 0:
            meas_ID = os.path.basename(prs_pcb_file[0])[:15]
            prs_pcb_qck_file = \
                glob.glob(os.path.join(prs_out, f'{meas_ID}*_pcb_qck_ATLAS_*_prepro.nc'))
            
    files_out = {'ray' : prs_ray_file, 
                 'tlc' : prs_tlc_file, 
                 'pcb' : prs_pcb_file,
                 'ray_qck' : prs_ray_qck_file, 
                 'tlc_qck' : prs_tlc_qck_file, 
                 'pcb_qck' : prs_pcb_qck_file}
    
    return(files_out)

def ray_test(ray_args, prs_files, ray_out, processing, reprocess = True):

    if processing['ray'] == True and len(prs_files['ray']) == 1:
    
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

def tlc_test(tlc_args, prs_files, tlc_out, processing, reprocess = True):
    
    if processing['tlc'] == True  and len(prs_files['tlc']) == 1:
        
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
        __telecover__(tlc_args, __version__)

    return()

def pcb_test(pcb_args, prs_files, pcb_out, processing, reprocess = True):
    
    if processing['pcb'] == True and len(prs_files['pcb']) == 1:
    
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

def quicklook(qck_args, prs_files, qck_out, quicklook, reprocess = True):

    if quicklook['ray'] == True and len(prs_files['ray_qck']) == 1:

        qck_ray_args = qck_args.copy()   
        
        qck_ray_args['input_file'] = prs_files['ray_qck'][0]   
        
        qck_ray_args = check_qck(qck_ray_args)
       
        qck_ray_args['output_folder'] = qck_out
        os.makedirs(qck_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(qck_out, 'plots', '*_ray_qck_*_ATLAS_*.png')):
                os.remove(file)
    
    
        # Excecute ATLAS visualizer
        view_qck(qck_ray_args)
        __quicklook__(qck_ray_args, __version__, meas_type = 'ray')
    
    elif quicklook['ray'] == True and len(prs_files['ray_qck']) == 0:
        print('-- Warning: Rayleight quicklook visualization was enabled but no file from the preprocessing stage was found. Please make sure that the quicklook option is enable in the preprocessor. If not, the measurement might need to be reprocessed (set newdata to True)')
        
    if quicklook['tlc'] == True and len(prs_files['tlc_qck']) == 1:
        
        qck_tlc_args = qck_args.copy() 

        if len(prs_files['tlc_qck']) == 1 and len(prs_files['tlc_qck']) == 1:
            qck_tlc_args['input_file'] = prs_files['tlc_qck'][0]        
        
        qck_tlc_args = check_qck(qck_tlc_args)
    
        qck_tlc_args['output_folder'] = qck_out
        os.makedirs(qck_out, exist_ok = True)
    
        if reprocess == True:
            for file in glob.glob(os.path.join(qck_out, 'plots', '*_tlc_qck_*_ATLAS_*.nc')):
                os.remove(file)
        
        # Excecute ATLAS visualizer
        view_qck(qck_tlc_args)
        __quicklook__(qck_tlc_args, __version__, meas_type = 'tlc')

    elif quicklook['tlc'] == True and len(prs_files['tlc_qck']) == 0:
        print('-- Warning: Telecover quicklook visualization was enabled but no file from the preprocessing stage was found. Please make sure that the quicklook option is enable in the preprocessor. If not, the measurement might need to be reprocessed (set newdata to True)')
        
    if quicklook['pcb'] == True and len(prs_files['pcb_qck']) == 1:
        
        qck_pcb_args = qck_args.copy()   

        if len(prs_files['pcb_qck']) == 1 and len(prs_files['pcb_qck']) == 1:
            qck_pcb_args['input_file'] = prs_files['pcb_qck'][0]   
        
        qck_pcb_args = check_qck(qck_pcb_args)    
        
        qck_pcb_args['output_folder'] = qck_out
        os.makedirs(qck_out, exist_ok = True)
        
        if reprocess == True:
            for file in glob.glob(os.path.join(qck_out, 'plots', '*_pcb_qck_*_ATLAS_*.nc')):
                os.remove(file)
            
        # Excecute ATLAS visualizer
        view_qck(qck_pcb_args)
        __quicklook__(qck_pcb_args, __version__, meas_type = 'pcb')

    elif quicklook['pcb'] == True and len(prs_files['pcb_qck']) == 0:
        print('-- Warning: Polarization calibration quicklook visualization was enabled but no file from the preprocessing stage was found. Please make sure that the quicklook option is enable in the preprocessor. If not, the measurement might need to be reprocessed (set newdata to True)')
                    
    return()