"""
@author: P. Paschou

Plotting the processed signals and the retrieved products 
"""
import os, re
import matplotlib.pyplot as plt
from helper_functions.helper_functions import find_xarray_index 
#from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def plot_ovl(dir_plt, ovl, plt_pack, 
                     info_sig=[], info_prod=[], dpi_v=150):

    iplots = ovl.iplots 
    ovl_plots = ovl.plt_map
    
    sig_4_plot = plt_pack['signal']
    ovl_4_plot = plt_pack['overlap']
    prod_4_plot = plt_pack['product']
    
    if iplots==1:    
        dir_plt_ovl=os.path.join(dir_plt, 'overlap')
        os.makedirs(dir_plt_ovl, exist_ok=True)
                
        #signal plots
        if len(sig_4_plot)>0:
            
            #range_ulim_idx = find_val_index(range_ulim, float, sig_4_plot.range)
            
            for i_d in sig_4_plot.id.values:
                
                fig,ax= plt.subplots(figsize=(7,8))
                plt_name = f'{i_d}.png'
                plt_title = f'Signal correction during overlap calculation iterations \n\n'+\
                            f'Channel: {i_d} ({int(info_sig.wave.loc[i_d])} nm)'
                
                full_ovl_idx = int(info_sig.loc[i_d,'full_ovl_idx'])

                for n in sig_4_plot.iteration.values:
                    ax.plot(sig_4_plot.loc[n,i_d,:full_ovl_idx].values,
                            sig_4_plot.range.loc[:full_ovl_idx].values,
                            label=f'n:{n}')
                
                ax.set(xlabel='Signal (A.U.)', ylabel = 'Range (m)', title=plt_title)
                ax.legend(loc=0)
                ax.minorticks_on()
                ax.ticklabel_format(axis='x', style = 'sci', scilimits=(0,0))
                
                fig.savefig(os.path.join(dir_plt_ovl, plt_name), 
                            bbox_inches = 'tight', dpi = dpi_v)
                
                plt.close(fig)
        
        #product plots
        if len(prod_4_plot)>0:
            #range_ulim_idx = find_val_index(range_ulim, float, prod_4_plot.range)
            
            for i_d in prod_4_plot.id.values:
                
                fig,ax= plt.subplots(figsize=(7,8))
                plt_name = f'{i_d}.png'
                plt_title = f'Signal correction during overlap calculation iterations \n\n'+\
                            f'Product: {info_prod.method.loc[i_d]} '+\
                            f'{info_prod.type_short.loc[i_d]} {info_prod.wave.loc[i_d]} nm'
                
                product_name = re.sub('_',' ',info_prod.type.loc[i_d].title())
                
                full_ovl_idx = int(info_sig.loc[i_d,'full_ovl_idx'])

                if info_prod.units.loc[i_d] == '':
                    units = ''
                else:
                    units =f'({info_prod.units.loc[i_d]})'
                
                for n in prod_4_plot.iteration.values:
                    ax.plot(prod_4_plot.loc[n,i_d,:full_ovl_idx].values,
                            prod_4_plot.range.loc[:full_ovl_idx].values,
                            label=f'n:{n}')
                
                ax.set(xlabel=f'{product_name} {units}', ylabel = 'Range (m)', title=plt_title)
                ax.legend(loc=0)
                ax.minorticks_on()
                ax.ticklabel_format(axis='x', style = 'sci', scilimits=(0,0))
                
                fig.savefig(os.path.join(dir_plt_ovl, plt_name), 
                            bbox_inches = 'tight', dpi = dpi_v)
                
                plt.close(fig)
        
        #overlap plots
        if len(ovl_4_plot)>0:
            
            #range_ulim_idx = find_val_index(range_ulim, float, ovl_4_plot.range)
            
            for i_d in ovl_4_plot.id.values:
                
                fig,ax= plt.subplots(figsize=(7,8))
                plt_name = f'{i_d}.png'
                plt_title = f'Overlap profile through the iterations \n\n'+\
                            f'overlap: {i_d} '
             
                full_ovl_idx = int(info_sig.loc[i_d,'full_ovl_idx'])

                for n in ovl_4_plot.iteration.values:
                    ax.plot(ovl_4_plot.loc[n,i_d,:full_ovl_idx].values,
                            ovl_4_plot.range.loc[:full_ovl_idx].values,
                            label=f'n:{n}')
                
                ax.set(xlabel=f'{i_d}', ylabel = 'Range (m)', title=plt_title)
                ax.legend(loc=0)
                ax.minorticks_on()
                
                fig.savefig(os.path.join(dir_plt_ovl, plt_name), 
                            bbox_inches = 'tight', dpi = dpi_v)
                
                plt.close(fig)
    
    return()

# def make_plot(fig, ax, xarray_plot, i_d, plt_name, plt_title, full_ovl):
    
#     for n in xarray_plot.iteration.values:
#         ax.plot(xarray_plot.loc[n,i_d,:full_ovl].values,
#                 xarray_plot.range.loc[:full_ovl].values,
#                 label=f'n:{n}')
    
#     ax.set(xlabel='Signal (A.U.)', ylabel = 'Range (m)', title=plt_title)
#     ax.legend(loc=0)
#     ax.minorticks_on()
#     ax.ticklabel_format(axis='x', style = 'sci', scilimits=(0,0))
        
#     return(fig, ax)