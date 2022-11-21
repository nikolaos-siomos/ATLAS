"""
@authors: Nikolaos Siomos (nikolaos.siomos@lmu.de)

================
General:
    Parses the command line arguments provided when calling __main_ 
    The information is stored in a python dictionary
    
Returns:
    
    args:
        A dictionary with all the information provided as command line arguments
        
"""

import argparse
import os

def call_parser():
        
    """Collects the information included as commandline arguments. 
    """
    
    print('Parsing master arguments...')    
    
    parser = argparse.ArgumentParser(
    	description='arguments ')
    

    parser.add_argument('-f', '--parent_folder', metavar='parent_folder', 
                        type=str, nargs='?',  default = None,
                        help='The path to the parent folder that contains the normal folder and all other optional input folders (dark, atmosphere, overlap). If no results folder is provided, it will be exported here by default')

    parser.add_argument('--dark_folder', metavar='dark_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the dark folder. Defaults to a drk folder inside the parent folder')

    parser.add_argument('--rayleigh_folder', metavar='rayleigh_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the rayleigh fit measurement folder. Defaults to a nrm folder inside the parent folder')
    
    parser.add_argument('--telecover_sectors_folder', metavar='telecover_sectors_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the telecover folder that contains the sector files. Defaults to a tlc/sectors folder inside the parent folder')

    parser.add_argument('--telecover_rings_folder', metavar='telecover_rings_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the telecover folder that contains the ring (inner/outer) files. Defaults to a tlc/rings folder inside the parent folder')
        
    parser.add_argument('--pcb_cal_p45_folder', metavar='pcb_cal_p45_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the polarization calibration +45 folder. Defaults to a pcb/+45 folder inside the parent folder')

    parser.add_argument('--pcb_cal_m45_folder', metavar='pcb_cal_m45_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the polarization calibration -45 folder. Defaults to a pcb/-45 folder inside the parent folder')

    parser.add_argument('--pcb_cal_stc_folder', metavar='pcb_cal_stc_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the polarization calibration folder for a calibration with a single calibrator position. Defaults to a pcb/stc folder inside the parent folder')

    parser.add_argument('--radiosonde_folder', metavar='radiosonde_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the radiosonde folder. The radiosonde file that is closest to the measurement within 12h will be selected if more than 1 files are provided inside')    
    
    parser.add_argument('--converter_out', metavar='converter_out', 
                        type=str, nargs='?', default = None,
                        help='The path to the folder where the converted netcdf files will be placed. This optional argument can be used if the folder must be placed out of the parent_folder. Defaults to parent_folder/netcdf ')

    parser.add_argument('--preprocessor_out', metavar='preprocessor_out', 
                        type=str, nargs='?', default = None,
                        help='The path to the folder where the preprocessed netcdf files will be placed. This optional argument can be used if the folder must be placed out of the parent_folder. Defaults to parent_folder/netcdf ')

    parser.add_argument('--visualizer_out', metavar='visualizer_out', 
                        type=str, nargs='?', default = None,
                        help='The path to the folder where the plots will be placed. This optional argument can be used if the folder must be placed out of the parent_folder. Defaults to parent_folder/plots')

    parser.add_argument('-c', '--config_file', metavar='config_file', 
                        type=str, nargs='?', default = None,
                        help='The path to the configuration file that contains the necessary metadata. This optional argument can be used if the settings folder must be placed out of the parent_folder (default)')            

    parser.add_argument('-s', '--settings_file', metavar='settings_file', 
                        type=str, nargs='?', default = None,
                        help='The absolute path to the settings file that contains options to run ATLAS. ')            

    args = vars(parser.parse_args())
    
    return(args)

def check_parser(args):

    if args['parent_folder'] != None:
        if not os.path.exists(args['parent_folder']) == True:
            raise Exception("The provided parent_folder {args['parent_folder']} does not exists. Please provide a valid path ")   

    if args['config_file'] == None:
        if args['parent_folder'] != None:
            args['settings_file'] = os.path.join(args['parent_folder'],'settings_file.ini')  
        else:
            raise Exception("Neither the parent folder nor the path to the settings file where provided. Please provide at least one of the two. The settings file path, unless specified defaults to ./<parent_folder>/settings_file.ini ")

    if not os.path.exists(args['settings_file']):
        raise Exception("The path to the settings file does not exist. Please provide a valid path. ")

    cnv_path = args['converter_out']
    prs_path = args['preprocessor_out']
    vis_path = args['visualizer_out']
    
    if cnv_path== None:
        cnv_path = os.path.join(args['parent_folder'], 'netcdf', 'converter')
        args['converter_out'] = cnv_path
        os.makedirs(cnv_path, exist_ok = True)
    elif os.path.isdir(os.path.split(os.path.abspath(cnv_path))[0]) or \
        os.path.isdir(os.path.split(os.path.split(os.path.abspath(cnv_path))[0])[0]) or \
            os.path.isdir(os.path.split(os.path.split(os.path.split(os.path.abspath(cnv_path))[0])[0])[0]):
        os.makedirs(cnv_path, exist_ok = True)        

    if prs_path== None:
        prs_path = os.path.join(args['parent_folder'], 'netcdf', 'preprocessor')
        args['preprocessor_out'] = prs_path
        os.makedirs(prs_path, exist_ok = True)
    elif os.path.isdir(os.path.split(os.path.abspath(prs_path))[0]) or \
        os.path.isdir(os.path.split(os.path.split(os.path.abspath(prs_path))[0])[0]) or \
            os.path.isdir(os.path.split(os.path.split(os.path.split(os.path.abspath(prs_path))[0])[0])[0]):
        os.makedirs(prs_path, exist_ok = True)          

    if vis_path== None:
        vis_path = os.path.join(args['parent_folder'], 'plots')
        args['visualizer_out'] = vis_path
        os.makedirs(vis_path, exist_ok = True)
    elif os.path.isdir(os.path.split(os.path.abspath(vis_path))[0]) or \
        os.path.isdir(os.path.split(os.path.split(os.path.abspath(vis_path))[0])[0]) or \
            os.path.isdir(os.path.split(os.path.split(os.path.split(os.path.abspath(vis_path))[0])[0])[0]):
        os.makedirs(vis_path, exist_ok = True)       

    return(args)

def substitute(org, rpl):
    
    for key in org.keys():
        if key in rpl.keys():
            if rpl[key] != None:
                if key != 'output_folder':
                    org[key] = rpl[key]
                else:
                    org[key] = os.path.join(rpl[key],'scc_converter')
                    
                
    return(org)