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

import argparse, os, sys

def call_parser(get_defaults = False):
        
    """Collects the information included as commandline arguments. 
    """
    
    print('Parsing master arguments...')    
    
    if get_defaults:
        sys.argv = [sys.argv[0]]
        
    parser = argparse.ArgumentParser(
    	description='arguments ')
    

    parser.add_argument('-f', '--parent_folder', metavar='parent_folder', 
                        type=str, nargs='?',  default = None,
                        help='The path to the parent folder that contains the normal folder and all other optional input folders (dark, atmosphere, overlap). If no results folder is provided, it will be created here by default')

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

    parser.add_argument('--radiosonde', metavar='radiosonde', 
                        type=str, nargs='?', default = None,
                        help='The path to the radiosonde file or the radiosonde folder. Provide either a folder path containing ascii radiosonde files or a netcdf file with the radiosonde data in the expected ATLAS (SCC) format. If a folder path is provided, the radiosonde file that is closest to the measurement within 12h will be selected if more than 1 files are provided inside')    

    parser.add_argument('--converter_out', metavar='converter_out', 
                        type=str, nargs='?', default = None,
                        help='The path to the folder where the converted netcdf files will be placed. This optional argument can be used if the folder must be placed out of the parent_folder. Defaults to parent_folder/netcdf ')

    parser.add_argument('--preprocessor_out', metavar='preprocessor_out', 
                        type=str, nargs='?', default = None,
                        help='The path to the folder where the preprocessed netcdf files will be placed. This optional argument can be used if the folder must be placed out of the parent_folder. Defaults to parent_folder/netcdf ')

    parser.add_argument('--visualizer_out', metavar='visualizer_out', 
                        type=str, nargs='?', default = None,
                        help='The path to the folder where the final plots and ascii files will be placed inside the ./visualizer_out/plots and ./visualizer_out/ascii folders. This optional argument can be used if the folder must be placed out of the parent_folder. Defaults to parent_folder')

    parser.add_argument('--ascii_folder', metavar='ascii_folder', 
                        type=str, nargs='?', default = None,
                        help='The path to the folder where the plots will be placed. This optional argument can be used if the folder must be placed out of the parent_folder. Defaults to parent_folder/plots')

    parser.add_argument('-c', '--config_file', metavar='config_file', 
                        type=str, nargs='?', default = None,
                        help='The path to the configuration file that contains the necessary metadata. This optional argument can be used if the settings folder must be placed out of the parent_folder (default)')            

    parser.add_argument('-s', '--settings_file', metavar='settings_file', 
                        type=str, nargs='?', default = None,
                        help='The absolute path to the settings file that contains options to run ATLAS. ')            

    parser.add_argument('--file_format', metavar='file_format', 
                        type=str, nargs='?', default = None,
                        help='The format of the raw lidar files. Currently licel and polly_xt are supported')            

    parser.add_argument('--operation_mode', metavar='operation_mode', 
                        type=str, nargs='?', default = 'labeling',
                        help='Choose one of: a) labeling: Use when submitting a measurement to CARS. Makes some SCC related field in the configuration file mandatory. b) testing: Use for causal testing of a measurement when there is no dedicated SCC configuration Defaults to: labeling')            

    parser.add_argument('--isday', metavar='isday', 
                        type=bool, default = False,
                        action = argparse.BooleanOptionalAction,
                        help='The format of the raw lidar files. Currently licel and polly_xt are supported')            

    parser.add_argument('--quick_run', metavar='quick_run', 
                        type=bool, default = False,
                        action = argparse.BooleanOptionalAction,
                        help='Defaults to: False. If set to True the converter and the preprocessing modules will not be called if the algorithm detects output files already produced by them for a specific measurement. This mainly saves time during execution.')            

    parser.add_argument('--process', metavar='process', 
                        type=str, nargs='?', default = ['ray','tlc','pcb','drk'],
                        help='The user can choose specific QA test(s) to process. It is only taken into account if visualize is set to True. Use any of: ray (Rayleigh Fit), tlc (Telecover Test), pcb (Polarization Calibration), off (Nothing will be processed). Defaults to: ray, tlc, pcb')            

    parser.add_argument('--process_qck', metavar='process_qck', 
                        type=str, nargs='*', default = ['ray','tlc','pcb','drk'],
                        help='Choose which quicklooks to process.  It is only taken into account if visualize is set to True. It must be a subset of process. Note that each visualizer module (ray, tlc, pcb) can create their own quicklooks. Use any of: ray (Rayleigh Fit), tlc (Telecover Test), pcb (Polarization Calibration), off (Nothing will be processed). Defaults to: ray tlc pcb')            

    parser.add_argument('--slice_rayleigh', metavar='slice_rayleigh', 
                        type=str, nargs='*', default = [None, None],
                        help='Provide temporal limits for the processing of the normal measurement. Use the following format: HHMM for the limits. For example use: --slice_rayleigh 2300, 0130 to slice between 23:00 UTC and 01:30 UTC (the next day) Defaults to: None, None ')            

    # parser.add_argument('-e', '--export_legacy', metavar = 'export_legacy',
    #                     type = bool, default = False, 
    #                     action = argparse.BooleanOptionalAction,
    #                     help = 'If called then data will be exported to the legacy QA file format of EARLINET. Defaults to False ')

    args = vars(parser.parse_args())
    
    return(args)

def check_parser(args):

    raw_fmts = ['licel','polly_xt','licel_matlab','polly_xt_first','licel_old2rack']
    if args['file_format'] == None:
        raise Exception("The mandatory argument file_format was not provided. Please add it with the --file_format option ")   
    elif args['file_format'] not in raw_fmts:
        raise Exception(f"The provided file_format: {args['file_format']} is not supported. Please choose among: {raw_fmts}")   
        
    if args['parent_folder'] != None:
        if not os.path.exists(args['parent_folder']):
            raise Exception(f"The provided parent_folder {args['parent_folder']} does not exists. Please provide a valid path ")   

    if args['settings_file'] == None:
        if os.path.exists(args['parent_folder']):
            args['settings_file'] = os.path.join(args['parent_folder'],'settings_file.ini')  
            if not os.path.exists(args['settings_file']):
                print('-----------------------------------------')
                print(f"-- Warning: The default path for the settings file inside the parent_folder {args['settings_file']} does not exist. Default settings will be applied! Please make sure this is what you want, otherwise provide a settings_file.ini ")
                print('-----------------------------------------')
        else:
            raise Exception("Neither the parent folder nor the path to the settings file where provided. Default settings will be applied! Please make sure this is what you want, otherwise provide the settings_file.ini path explicitely ")
    else:
        if not os.path.exists(args['settings_file']):
            raise Exception(f"The provided path to the settings file {args['settings_file']} does not exist. Please provide a valid path or skip it entirely to assume default settings ")
        
    # if args['settings_file'] == None:
    #     if args['parent_folder'] != None:
    #         args['settings_file'] = os.path.join(args['parent_folder'],'settings_file.ini')  
    #     else:
    #         raise Exception("Neither the parent folder nor the path to the settings file where provided. Please provide at least one of the two. The settings file path, unless specified defaults to ./<parent_folder>/settings_file.ini ")

    # if not os.path.exists(args['config_file']):
    #     raise Exception(f"The path to the config file {args['config_file']} does not exist. Please provide a valid path. ")


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
        vis_path = os.path.join(args['parent_folder'])
        args['visualizer_out'] = vis_path
        os.makedirs(vis_path, exist_ok = True)
    elif os.path.isdir(os.path.split(os.path.abspath(vis_path))[0]) or \
        os.path.isdir(os.path.split(os.path.split(os.path.abspath(vis_path))[0])[0]) or \
            os.path.isdir(os.path.split(os.path.split(os.path.split(os.path.abspath(vis_path))[0])[0])[0]):
        os.makedirs(vis_path, exist_ok = True)       
        
    if args['radiosonde'] == None:
        raise Exception("The radiosonde argument is mandatory! Please provide it. ")
        
    if os.path.isdir(args['radiosonde']) == False and os.path.isfile(args['radiosonde']) == False:
        raise Exception("The radiosonde folder/file path does not exist. Please provide an existing path. ")

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