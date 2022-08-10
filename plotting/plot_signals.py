"""
@author: P. Paschou & N. Siomos

Plot the processed signals
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from helper_functions.helper_functions import datetime_unit
from readers.read_maps import plt_options
import matplotlib as mpl
import math

# mpl.rcParams.update(mpl.rcParamsDefault)
ln_wdth = 2.
ax_lblsz= 18
ax_tlsz = 20
xtck_lblsz = ytck_lblsz = 18
ax_lblpd = 15
lgd_fsz = 20


lcolor=pd.Series(data=['dodgerblue','aqua','limegreen','yellowgreen','red', 'magenta', 'salmon', 'royalblue'], 
                index=['355', '387', '532', '607', '1064', '355/532', '355/1064', '532/1064'])


#Plotting the range corrected signals in figures per channel per time
def plot_sig(cfg, maps_dir, sig,  sig_alt, info_sig, isdark=0):  

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        # Extract the desired info from config file
        lidar_info = cfg.lidar
        iplt = cfg.iplt
        dir_plt = cfg.directories.plots
        isig_hscale = cfg.isig_hscale
        
        if iplt in [1,3] and len(sig)>0:
    
            #read the plotting options file
            options = plt_options(maps_dir, 'plot_options.ini','Signals')        
    
            # Extract the plotting options
            dpi_v = options.dpi
            lcolor = pd.Series(data=options.wl_color, index=options.wl)
            ymin = options.min_alt/1e3 # m to km
            ymax = options.max_alt/1e3 # m to km
            ystep = options.step_alt/1e3 # m to km
            y_ticks = np.arange(ymin, ymax+ystep, ystep)
            
    
            # Define the height scale in km
            if isig_hscale =='altitude':
                ytype = 'altitude a.s.l.'
                sig_h = sig.altitude.values/1e3                    
            else:
                ytype = 'range'
                sig_h = sig.range.values/1e3        
    
            vis_sig = options.vis_sig
            if len(vis_sig) == 0:
                vis_sig = sig.channel.values
    
            for t in range(sig.time.size):
                ind_t = dict(time = t)            
    
                date = sig[ind_t].time.dt.strftime('%Y-%m-%d').values        
                stime = sig[ind_t].time.dt.strftime('%H:%M').values        
                dt_str = sig[ind_t].time.dt.strftime('%Y%m%dT%H%M%SZ').values        
        
                # Calculate the ending time of the timeframe        
                t_unit = datetime_unit(cfg.timescale).lower()
                number=int(''.join(filter(str.isdigit, cfg.timescale)))        
        
                t2 = sig[ind_t].time+np.timedelta64(number,t_unit)        
                stime2 = t2.dt.strftime('%H:%M').values        
        
                for c in vis_sig:
                    
                    if not all(np.isnan(sig[ind_t].loc[dict(channel = c)].values)):
        
                        magnitude = math.floor(math.log(np.nanmax(sig[ind_t].loc[dict(channel = c)].values),10))
                        max_lim = int(np.nanmax(sig[ind_t].loc[dict(channel = c)].values)/10**magnitude +1) * (10**magnitude)
                        xstep = max_lim/10 # to create 10 xticks
                        x_ticks = np.arange(0, max_lim + xstep, xstep)
            
                        
                        output_dir = os.path.join(dir_plt, 'pre_processing')            
                        os.makedirs(output_dir, exist_ok=True)
            
                        plt_name = f'rc_{c}_{dt_str}.png'
                        plt_title = f'{lidar_info.lr_name} @ {lidar_info.location} \n'+\
                                    f'Time frame (UTC): {date}, {stime} - {stime2} \n' +\
                                    f'Channel: {c}'
                        
                        fig, ax = plt.subplots(figsize=(9,12))
                
                        ax.plot(sig[ind_t].loc[dict(channel = c)].values, sig_h, 
                                color = lcolor[str(int(info_sig.wave[c]))], 
                                label = f'{int(info_sig.wave[c])} nm')
                        ax.set(xlabel = 'Signal [A.U.]', ylabel = f'{ytype.capitalize()} [km]', 
                               title = plt_title,
                               xlim = [0, max_lim], xticks = x_ticks,
                               ylim = [ymin,ymax], yticks = y_ticks)
                        ax.ticklabel_format(axis='x', style = 'sci', scilimits=(0,0))
                        ax.minorticks_on()                    
                        ax.legend(loc=1)
                        
                        fig.savefig(os.path.join(output_dir, plt_name), 
                                    bbox_inches = 'tight', dpi = dpi_v)
                        
                        plt.close(fig)
            if isdark == 0:
                printer = 'normal'
            if isdark == 1:
                printer = 'dark'
            
            print(f'-- Plotting of range-corrected {printer} signals complete!')
    
    return()
