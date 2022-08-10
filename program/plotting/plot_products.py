"""
@author: P. Paschou & N. Siomos

Plot the retrieved products 
"""
import os, re, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from helper_functions.helper_functions import datetime_unit
from readers.read_maps import plt_options
import matplotlib as mpl

ln_wdth = 2.
ax_lblsz = 18#15
ax_tlsz = 19#15
xtck_lblsz = ytck_lblsz = 18#15
ax_lblpd = 18#15.0
lgd_fsz = 18#18

# Plot gridlines in plots?
gd = True
gd_c = 'grey'
gd_l = 'dashed'


def plot_prod(cfg, maps_dir, prod_map, prod, info_prod):

    # Extract the desired info from config file
    iplt = cfg.iplt
    
    if iplt > 1 and len(prod)>0:
        #create dir
        dir_plt = cfg.directories.plots
        output_dir = os.path.join(dir_plt, 'products')            
        os.makedirs(output_dir, exist_ok=True)        
        
        #read the plotting options file
        options = plt_options(maps_dir, 'plot_options.ini','Products')        
        
        #meke plots
        plot_prod_per_wl(output_dir, cfg, prod_map, prod, info_prod, options)

        plot_prod_all_wl(output_dir, cfg, prod_map, prod, info_prod, options)        
        
    return()

#Plotting the products in figures per product per time
def plot_prod_per_wl(output_dir, cfg, prod_map, prod, info_prod, options):

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         
    
        # Extract the desired info from config file
        lidar_info = cfg.lidar
        timescale = cfg.timescale
        iprod_hscale = cfg.iprod_hscale
            
        if options.plt_type in [0,2]:
    
            # Extract the plotting options
            dpi_v = options.dpi        
            vis_prod = options.vis_prod
            # if plot from prod_map empty, plot all products
            if len(vis_prod)==0:
                vis_prod = prod.product.values
            
            lcolor = pd.Series(data=options.wl_color, index=options.wl)
            lim = pd.DataFrame(data=[options.min_lim,options.max_lim],
                               index=['low', 'up'],
                               columns = options.prod_type,
                               dtype = float)
            step = pd.Series(data = options.step, index = options.prod_type, dtype=float)
            
            y_min = options.min_alt/1e3 # m to km
            y_max = options.max_alt/1e3 # m to km
            ystep = options.step_alt/1e3 # m to km
            y_ticks = np.arange(y_min, y_max+ystep, ystep)
            
            # Define the height scale in km
            if iprod_hscale =='altitude':
                ytype = 'altitude a.s.l.'
                prod_h = prod.altitude.values/1e3                    
            else:
                ytype = 'range'
                prod_h = prod.range.values/1e3        
            
            for t in range(prod.time.size):
                ind_t = dict(time = t)  
    
                date = prod[ind_t].time.dt.strftime('%Y-%m-%d').values
                stime = prod[ind_t].time.dt.strftime('%H:%M').values        
                dt_str = prod[ind_t].time.dt.strftime('%Y%m%dT%H%M%SZ').values        
        
                # Calculate the ending time of the timeframe        
                t_unit = datetime_unit(timescale).lower()
                number=int(''.join(filter(str.isdigit, timescale)))        
        
                t2 = prod[ind_t].time+np.timedelta64(number,t_unit)        
                stime2 = t2.dt.strftime('%H:%M').values        
                
                for c in vis_prod:
      
                    if not all(np.isnan(prod[ind_t].loc[dict(product = c)].values)):
    
                        #product_name = re.sub('_',' ',info_prod.type[c].capitalize()) #for capitalize only the first letter 
                        product_name = re.sub('_',' ',info_prod.type[c].title()) # title() : capitalize first letter of each word
                        method_name = info_prod.method[c].title()
                        full_product_name = ' '.join([method_name, product_name])
            
                        plt_name = f'{c}_{dt_str}.png'
                        plt_title = f'{lidar_info.lr_name} @ {lidar_info.location} \n'+\
                                    f'Time frame (UTC): {date}, {stime} - {stime2} \n' +\
                                    f'{full_product_name.lstrip()}'
        
                        xstep = step[info_prod.type_short[c]]
                        llim = lim[info_prod.type_short[c]].low
                        ulim = lim[info_prod.type_short[c]].up
                        x_ticks = np.arange(np.round(llim, decimals=3), ulim + xstep, xstep)
                        
                        if info_prod.units[c] == '':
                            units = ''
                        else:
                            units =f'[{info_prod.units[c]}]'
                        
                        fig, ax = plt.subplots(figsize=(10,10))
                
                        ax.plot(prod[ind_t].loc[dict(product = c)].values, prod_h, linewidth=1.2,
                                color = lcolor[str(info_prod.wave[c])], label = f'{info_prod.wave[c]} nm')
                        ax.set(xlabel = f'{product_name} {units}', ylabel = f'{ytype.title()} [km]', 
                               title = plt_title,
                               xlim = [llim, ulim],
                               ylim = [y_min, y_max],
                               xticks = x_ticks, yticks = y_ticks)
                        
                        ax.axvline(x = 0., ymin = y_min, ymax = y_max,
                                   color='silver', linestyle='--',alpha=0.5)
                        if  gd:
                            ax.grid(which='major',axis='both', alpha=0.2, linestyle=gd_l, color=gd_c)

                        ax.minorticks_on()
                        ax.legend(loc=1)
                        
                        if ulim <=1e-3 or ulim >=1e3: 
                            ax.ticklabel_format(axis='x', style = 'sci', scilimits=(0,0))
                        
                        fig.savefig(os.path.join(output_dir, plt_name), 
                                    bbox_inches = 'tight', dpi = dpi_v)
                        
                        plt.close(fig)
    
        print('-- Plotting of all products (*one plot per product per wavelength) complete!')

    return()

# #Plotting the products in figures product type (i.e. bsc, ext, vldr), per method(i.e. klett, raman) and per time
def plot_prod_all_wl(output_dir, cfg, prod_map, prod, info_prod, options):
    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        # Extract the desired info from config file
        lidar_info = cfg.lidar
        timescale = cfg.timescale
        iprod_hscale = cfg.iprod_hscale
    
        if options.plt_type > 0:
            output_dir = os.path.join(output_dir, 'all')            
            os.makedirs(output_dir, exist_ok=True)
    
            # Extract the plotting options
            dpi_v = options.dpi
            lcolor = pd.Series(data=options.wl_color, index=options.wl)
            lim = pd.DataFrame(data=[options.min_lim,options.max_lim],
                               index=['low', 'up'],
                               columns = options.prod_type,
                               dtype = float)
            step = pd.Series(data = options.step, index = options.prod_type, dtype=float)
            
            y_min = options.min_alt/1e3 # m to km
            y_max = options.max_alt/1e3 # m to km
            ystep = options.step_alt/1e3 # m to km
            y_ticks = np.arange(y_min, y_max+ystep, ystep)
            
            # Define the height scale in km
            if iprod_hscale =='altitude':
                ytype = 'altitude a.s.l.'
                prod_h = prod.altitude.values/1e3                    
            else:
                ytype = 'range'
                prod_h = prod.range.values/1e3        
    
            # Select only desired products from product_map
            vis_prod = options.vis_prod
            # if plot from prod_map empty, plot all products
            if len(vis_prod)==0:
                vis_prod = prod.product.values
            prod = prod.copy().sel(product = vis_prod)
            info_prod = info_prod.loc[vis_prod,:]
    
            group_by_type = info_prod.groupby('type_short').groups # dict with keys:type and values:indexes
            #group_by_type = pd.DataFrame(group_by_type.values(), index = group_by_type.keys())
    
            for prod_type in group_by_type.keys():
                #if group_by_type[prod_type].shape[0]>1:
                
                group_by_method = info_prod.loc[group_by_type[prod_type]].groupby('method').groups
                
                for m in group_by_method.keys():
                    
                    for t in range(prod.time.size):
                        ind_t = dict(time = t)  
    
                        date = prod[ind_t].time.dt.strftime('%Y-%m-%d').values
                        stime = prod[ind_t].time.dt.strftime('%H:%M').values        
                        dt_str = prod[ind_t].time.dt.strftime('%Y%m%dT%H%M%SZ').values        
                
                        # Calculate the ending time of the timeframe        
                        t_unit = datetime_unit(timescale).lower()
                        number=int(''.join(filter(str.isdigit, timescale)))        
                
                        t2 = prod[ind_t].time+np.timedelta64(number,t_unit)        
                        stime2 = t2.dt.strftime('%H:%M').values        
                        
                        # Create the figure per product type, per method and per time
                        fig, ax = plt.subplots(figsize=(6,8)) 
                        skip_fig = 1 # if no product is ploted then the fig is not saved
                        
                        for p in group_by_method[m]:
                            
                            if not all(np.isnan(prod[ind_t].loc[dict(product = p)].values)):
                                ax.plot(prod[ind_t].loc[dict(product = p)].values, prod_h, 
                                        color = lcolor[str(info_prod.wave[p])], 
                                        label = f'{info_prod.wave[p]} nm')
                                ax.axvline(x = 0., ymin = y_min,ymax = y_min,
                                           color='silver', linestyle='--',alpha=0.5)
        
                                skip_fig = 0
        
                        if  gd:
                            ax.grid(which='major',axis='both', alpha=0.2, linestyle=gd_l, color=gd_c)
                        
                        #set the plot axes and title
                        product_name = re.sub('_',' ',info_prod.type[p].title()) # title() : capitalize first letter of each word
                        full_product_name = ' '.join([m.title(), product_name])
            
                        plt_name = f'{prod_type}_{m}_{dt_str}.png'
                        plt_title = f'{lidar_info.lr_name} @ {lidar_info.location}\n'+\
                                    f'Time frame (UTC): {date}, {stime} - {stime2}\n' +\
                                    f'{full_product_name.lstrip()}'
        
                        xstep = step[info_prod.type_short[p]]
                        llim = lim[info_prod.type_short[p]].low
                        ulim = lim[info_prod.type_short[p]].up
                        x_ticks = np.arange(np.round(llim, decimals=3), ulim + xstep, xstep)
                        
                        if info_prod.units[p] == '':
                            units = ''
                        else:
                            units = f'[{info_prod.units[p]}]'
                        
                        
                        if skip_fig == 0: # save the fig if there are plots of the products
                            ax.set(xlabel = f'{product_name} {units}', ylabel = f'{ytype.title()} [km]', 
                                   title = plt_title,
                                   xlim = [llim, ulim],
                                   ylim = [y_min, y_max],
                                   xticks = x_ticks, yticks = y_ticks)
                            ax.minorticks_on()
                            ax.legend(loc=1)
                            if ulim <=1e-3 or ulim >=1e3: 
                                ax.ticklabel_format(axis='x', style = 'sci', scilimits=(0,0))
                                            
                            
                            fig.savefig(os.path.join(output_dir, plt_name), 
                                        bbox_inches = 'tight', dpi = dpi_v)
                        
                        plt.close(fig)
    
                    
        print('-- Plotting of all products (*one plot per product type) complete!')

    return()