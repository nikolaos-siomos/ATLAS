"""
@author: P. Paschou & N. Siomos

Readers for the map files in settings folder
"""
import os, re
import configparser
import pandas as pd
from .read_config import comma_split, read_var, check_len

class gluing():
    def __init__(self, folder, filename):
        
        path = check_path(os.path.join(folder, filename))
        
        parser = configparser.ConfigParser()
        parser.read(path)
        
        self.map = read_section(parser['Signals_ID'], object)
        
        # Signal Threshold
        self.l_sig_min = read_var(parser['Thresholds']['low_sig_min'], float)
        self.u_sig_max = read_var(parser['Thresholds']['up_sig_max'], float)
        # Range Threshold
        self.l_range_min = read_var(parser['Thresholds']['low_range_min'], float)         
        self.u_range_max = read_var(parser['Thresholds']['up_range_max'], float)         
        # Gluing Metrics
        self.ref_sig = parser['Metrics']['sig_ref']
        self.hwin_len = read_var(parser['Metrics']['searching_hwindow'], float)
        self.cor_coef = read_var(parser['Metrics']['correlation_coef'], float)
        
class molecular():
    def __init__(self, folder, filename):
        
        path = check_path(os.path.join(folder, filename))
        
        parser = configparser.ConfigParser()
        parser.read(path)
        
        self.elastic = read_section(parser['Elastic'], object)
        self.depol = read_section(parser['Depol'], object)
        self.raman = read_section(parser['Raman'], object)

#read map in .ini format and export class instances
class total_sig():

    def __init__(self, folder, filename):
        path = check_path(os.path.join(folder, filename))
        
        parser = configparser.ConfigParser()
        parser.read(path)
        #put the signals ID in a DataFrame
        self.map = read_section(parser['Signals_ID'], object)
        #put the Parameters in a DataFrame
        self.params = read_section(parser['Parameters'], float)
            

#read map in .ini format and export class instances
class products():
    
    def __init__(self, folder, filename):
        path = check_path(os.path.join(folder, filename))    
        
        parser = configparser.ConfigParser()
        parser.read(path)
        
        #maps in dataframews
        
        # klett backscatter  
        if parser.has_section('Klett_bsc'):
            self.kl_bsc = read_section(parser['Klett_bsc'], object)
        else:
            self.kl_bsc = []
        
        # raman extinction
        if parser.has_section('Raman_ext'):
            self.rm_ext = read_section(parser['Raman_ext'], object, skip_vars=['deriv_lims','deriv_win','deriv_order'])
            lims = comma_split(parser['Raman_ext']['deriv_lims'], float)
            win = comma_split(parser['Raman_ext']['deriv_win'], float)        
            order = read_var(parser['Raman_ext']['deriv_order'], int) 
            self.deriv = {'lims':lims, 'win':win, 'order':order}
        else:
            self.rm_ext = []
        # raman backscatter
        if parser.has_section('Raman_bsc'):
            self.rm_bsc = read_section(parser['Raman_bsc'], object)
        else:
            self.rm_bsc = []
            
        #scattering ratio
        if parser.has_section('Scattering_ratio'):
            self.sc_r = read_section(parser['Scattering_ratio'], object)
        else:
            self.sc_r = []

        #attenuated backscatter ratio
        if parser.has_section('Normalized_Attenuated_Backscatter'):
            self.atten_nrm = \
                read_section(parser['Normalized_Attenuated_Backscatter'], object)
        else:
            self.atten_nrm = []
                    
        #volume depol ratio
        if parser.has_section('VDR'):
            self.vdr = read_section(parser['VDR'], object)
        else:
            self.vdr = []

        #volume depol ratio from Polly (cross and total sig)
        if parser.has_section('VDR_XT'):
            self.vdr_xt = read_section(parser['VDR_XT'], object)
        else:
            self.vdr_xt = []
        
        #particle depol ratio
        if parser.has_section('PDR'):
            self.pdr = read_section(parser['PDR'], object)
        else:
            self.pdr = []
        
        # Lidar ratio    
        if parser.has_section('Lidar_ratio'):
            self.lr = read_section(parser['Lidar_ratio'], object)
        else:
            self.lr = []

        # Angstrom exponent (bae or eae)
        if parser.has_section('Angstrom_exponent'):
            self.angstr = read_section(parser['Angstrom_exponent'], object)
        else:
            self.angstr = []
            
        # Depolarization Ratio Conversion
        if parser.has_section('DR_Conversions'):
            self.dr_conv = read_section(parser['DR_Conversions'], object)
        else:
            self.dr_conv = []
            
        # Backscatter total/parallel/cross Conversion
        if parser.has_section('BSC_Conversions'):
            self.bsc_conv = read_section(parser['BSC_Conversions'], object)
        else:
            self.bsc_conv = []

        # Parallel/Cross corrected Signal Conversion
        if parser.has_section('SIG_Conversions'):
            self.sig_conv = read_section(parser['SIG_Conversions'], object)
        else:
            self.sig_conv = []
            
        # Exports
        if parser.has_section('Export'):
            self.export = comma_split(parser['Export']['list'], str)
        else:
            self.export = []
            
#read map in .ini format and export class instances
class overlap():
    
    def __init__(self, folder, filename, iovl):

        calc_map, cor_map = [], []
        ovl_file, change_lr = [], []
        iovl_plots, ovl_plots = [], []

        #if iovl in [1,2]:
        if iovl > 0:
            path = check_path(os.path.join(folder, filename))    
            parser = configparser.ConfigParser()
            parser.read(path)
            
            #Read the correction module related parameters 
            cor_map = read_section(parser['Correction'], object, skip_vars=['ovl_file'])
            if iovl == 1:
                ovl_file = str(parser['Correction']['ovl_file'])
            
            if iovl == 2: #Read the map inside the calculation module
                change_lr = comma_split(parser['Calculation']['lidar_ratio_region'], float)
                calc_map = read_section(parser['Calculation'], object, skip_vars=['lidar_ratio_region'])
                
            iovl_plots = int(parser['Plotting']['iovl_plots'])
            if iovl_plots == 1:
                ovl_plots = read_section(parser['Plotting'], object, skip_vars=['iovl_plots'],return_empty=0)
            
        self.filepath = ovl_file
        self.cor_map = cor_map
        self.calc_map = calc_map
        self.lr_region = change_lr
        self.iplots = iovl_plots
        self.plt_map = ovl_plots

#read map in .ini format and export class instances
class export():
    
    def __init__(self, folder, filename):
        path = check_path(os.path.join(folder, filename))    

        parser = configparser.ConfigParser()
        parser.read(path)
        
        # ASSA 
        if parser.has_section('ASSA'):
            self.assa_prof = comma_split(parser['ASSA']['profiles'], str)
            self.assa_mol = comma_split(parser['ASSA']['molecular'], str)
        else:
            self.assa_prof = []
            self.assa_mol = []
        
#read map in .ini format and export class instances        
class plt_options():
    def __init__(self,folder, filename, module):
        path = check_path(os.path.join(folder, filename))    
        
        parser = configparser.ConfigParser()
        parser.read(path)

        # rfit, prod = [], []
        # telc, qckl = [], []

        if parser.has_section(module):

            if module == 'Rayleigh_Fit':   
                self.rfit = read_section(parser['Rayleigh_Fit'], object, return_empty=0)
                
                
            if module == 'Products':
                self.dpi = int(parser['Products']['dpi'])
                self.vis_prod = comma_split(parser['Products']['vis_prod'], object)
                self.plt_type = int(parser['Products']['plt_type'])
                
                self.wl = comma_split(parser['Products']['wl'], object)
                self.wl_color = comma_split(parser['Products']['wl_color'], object)
                check_len(self.wl, 'wl', self.wl_color, 'wl_color', 'plot_options.ini')
                
                self.min_alt = float(parser['Products']['min_alt'])
                self.max_alt = float(parser['Products']['max_alt'])
                self.step_alt = float(parser['Products']['step_alt'])
                
                self.prod_type = comma_split(parser['Products']['prod_type'], object)
                self.min_lim = comma_split(parser['Products']['min_lim'], float)
                check_len(self.prod_type, 'prod_type', self.min_lim, 'min_lim', 'plot_options.ini')
                
                self.max_lim = comma_split(parser['Products']['max_lim'], float)
                check_len(self.prod_type, 'prod_type', self.max_lim, 'max_lim', 'plot_options.ini')

                self.step = comma_split(parser['Products']['step'], float)
                check_len(self.prod_type, 'prod_type', self.step, 'step', 'plot_options.ini')                

                
            if module == 'Signals':
                self.dpi = int(parser['Products']['dpi'])

                self.vis_sig = comma_split(parser['Signals']['vis_sig'], object)                                
                self.wl = comma_split(parser['Signals']['wl'], object)
                self.wl_color = comma_split(parser['Signals']['wl_color'], object)
                check_len(self.wl, 'wl', self.wl_color, 'wl_color', 'plot_options.ini')

                self.min_alt = float(parser['Signals']['min_alt'])
                self.max_alt = float(parser['Signals']['max_alt'])
                self.step_alt = float(parser['Signals']['step_alt'])

                
            if module == 'Telecover':
                self.telc = read_section(parser['Telecover'], object, return_empty=0)


            if module == 'Calibration':
                self.clb = read_section(parser['Calibration'], object, return_empty=0)

            
            if module == 'Quicklooks':
                skip_prms = ['vis_sig', 'min_sig', 'max_sig', 'step_sig',
                             'vis_prod', 'min_prod', 'max_prod', 'step_prod']
                
                self.qcklk = read_section(parser['Quicklooks'], object, 
                                          skip_vars=skip_prms, return_empty=0)
                
                #read the rest params
                self.qlk_vis_sig = comma_split(parser['Quicklooks']['vis_sig'], object)
                self.qlk_min_sig = comma_split(parser['Quicklooks']['min_sig'], object)
                self.qlk_max_sig = comma_split(parser['Quicklooks']['max_sig'], object)
                self.qlk_step_sig = comma_split(parser['Quicklooks']['step_sig'], object)

                self.qlk_vis_prod = comma_split(parser['Quicklooks']['vis_prod'], object)
                self.qlk_min_prod = comma_split(parser['Quicklooks']['min_prod'], object)
                self.qlk_max_prod = comma_split(parser['Quicklooks']['max_prod'], object)
                self.qlk_step_prod = comma_split(parser['Quicklooks']['step_prod'], object)
                            
        else:
            raise ModuleNotFoundError(f'Module {module} not found in Plotting_options.ini settings file. Programm stopped')
        

    
#-------------------------- HELPER FUNCTIONS ----------------------------------
def read_section(section, dtype=object, skip_vars=[], return_empty=1):
    # Reads the whole or part of the section and returns a Pandas Series
    first = 1
    for key in section:
        if key not in skip_vars:
            temp = pd.Series([i.strip() for i in re.split(',', section[key]) if i !=''], 
                             dtype = dtype, name = key)
            #temp = pd.Series([i for i in section[key].split(sep=', ')], dtype = dtype, name = key)
            if first == 1:
                map_info = temp
                first = 0
            else:
                map_info = pd.concat([map_info, temp], axis = 1)
    
    # if return empty is 1 and the map has at least one NaN value then return an empty map
    if return_empty==1 and map_info.isnull().values.any():
    #if return_empty==1 and not map_info.all(axis = None, skipna = False): 
        map_info = []

    return(map_info)


def check_path(path):
    if not os.path.exists(path):
        folder, filename = os.path.split(path)
        raise OSError(f'The file {filename} does not exist in {folder}. Programm stopped')
    return(path)
