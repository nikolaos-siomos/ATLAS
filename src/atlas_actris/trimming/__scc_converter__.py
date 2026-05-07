"""
@author: N. Siomos & P. Paschou

Main algorithm for pre-processing the raw signals and retrieve the optical products
"""
import warnings, os, sys
from .readers.parse_args import call_parser, check_parser
from .tools import process, automate

warnings.filterwarnings('ignore')

def main(args, __version__, dry_run = False):
        
    print('-----------------------------------------')
    print('Initializing Converter...')
    print('-----------------------------------------')
    print(' ')
    
    # Identify the measurement type (rayleigh , telecover, or polarization_calibration)    
    meas_type = automate.get_meas_type(args)

    # In[1]
    #------------------------------------------------------------
    # A) Read and pre-process the signals
    #------------------------------------------------------------
    allowed_types = ['rayleigh', 'telecover', 
                     'polarization_calibration', 'dark'] 
    
    processors = {'rayleigh' : process.rayleigh,
                  'telecover' : process.telecover,
                  'polarization_calibration' : process.polarization_calibration,
                  'dark' : process.dark}
    
    modes = {'rayleigh' : 'R',
             'telecover' : 'T',
             'polarization_calibration' : 'C',
             'dark' : 'D'}
    
    output_files = {'rayleigh' : None,
                    'telecover' : None,
                    'polarization_calibration' : None,
                    'dark' : None,
                    'radiosonde' : None} 
    
    # Call all the processors sequentially
    for mtype in allowed_types: 
        if mtype in meas_type and (args['mode'] == 'A' or args['mode'] == modes[mtype]):
            
            print(f"Processing {mtype} measurement")
            
            nc_fname = processors[mtype](args, version = __version__)
            
            if nc_fname[0] != None:
                if mtype == 'rayleigh':
                    args['rayleigh_filename'] = os.path.basename(nc_fname[0])
                    args['radiosonde_filename'] = os.path.basename(nc_fname[1])
                    output_files['rayleigh'] = nc_fname[0]
                    output_files['radiosonde'] = nc_fname[0]
                else:
                    output_files[mtype] = nc_fname[0]
        
        elif mtype not in meas_type and args['mode'] == 'A':
            print(f"--Warning: No {mtype} files were processed!")
            print("")
            
        elif mtype not in meas_type and args['mode'] == modes[mtype]:
            print(f"--Warning: No {mtype} folder was detected despite the provided mode ({modes[mtype]})!")
    
    return(output_files)

if __name__ == '__main__':
    
    sys.path.append('../')
    from version import __version__
    
    # Get the command line argument information
    args = call_parser()
    
    # Call main
    main(args, __version__)
