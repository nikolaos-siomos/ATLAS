"""
@author: P. Paschou

Plotting routines for Quicklooks
"""
import os, re
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from plotting import make_colormap
from helper_functions.helper_functions import check_class_attr


ln_wdth = 2.
ax_lblsz = 15
ax_tlsz = 15
xtck_lblsz = ytck_lblsz =15
ax_lblpd = 15
lgd_fsz = 20


# figure size
x_inches = 12
y_inches = 6

def sig_quicklook(cfg, sig, info, options, ver=''):

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        colormap_style = select_colormap(options.qcklk) # make_colormap.cmap # 'Spectral_r' # 
        
        ql_dir = cfg.directories.plots
        isig = int(options.qcklk.isig[0])
    
        if isig == 1 and len(sig) > 0:
            # Create dir if not exist
            # out_dir = os.path.join(ql_dir,'pre_processing')
            out_dir = os.path.join(ql_dir,'quicklooks')            
            os.makedirs(out_dir, exist_ok=True)
    
            # Extract info from config file
            sza = np.round(float(cfg.angles.zenith), decimals=1)
            location = cfg.lidar.location
            lat = np.round(float(cfg.lidar.latitude), decimals=1)
            lon = np.round(float(cfg.lidar.longitude), decimals=1)
            elev = int(np.round(float(cfg.lidar.altitude), decimals=0))
            lidar = cfg.lidar.lr_name
            isig_hscale = cfg.isig_hscale
            
            loc_info = f'{location} (lat: {lat}, lon: {lon}, elev: {elev} m)'
    
            # Extract plotting options
            dpi_val = int(options.qcklk.dpi[0])    
            inorm = int(options.qcklk.inorm[0])
            itime = int(options.qcklk.itime[0])
            check_class_attr(options.qcklk,'scale_sig','linear')
            scale_sig = options.qcklk.scale_sig[0]
            # check_class_attr(options.qcklk,'scale_sig','linear')
            # color_map = 
            
            vis_sig = options.qlk_vis_sig
            min_sig = options.qlk_min_sig
            max_sig = options.qlk_max_sig
            step_sig = options.qlk_step_sig
        
            # datetime title settings
            if itime == 1:
                tstart = mdates.date2num(np.datetime64(options.qcklk.tstart[0]))
                tstop = mdates.date2num(np.datetime64(options.qcklk.tstop[0]))
            
            date = sig.time[0].dt.strftime('%Y/%m/%d').values        
            date_out = sig.time[0].dt.strftime('%Y%m%d').values
            t1 = sig.time[0].dt.strftime('%H:%M').values
            t1_out = sig.time[0].dt.strftime('%H%M').values
            t2 = sig.time[-1].dt.strftime('%H:%M').values
            t2_out = sig.time[-1].dt.strftime('%H%M').values        
            
            # convert meas timeframe to numbers (easier to plot)
            date_num = mdates.date2num(sig.time.values)
            
            # formating the time in x label
            tstep = int(options.qcklk.step_time[0])
            minutes = mdates.MinuteLocator(interval = tstep)#4 #20
            t_fmt = mdates.DateFormatter('%H:%M')
            
    
            # Select the desired products for plotting
            ch_list = vis_sig
            sig_plt = sig.loc[{'channel':ch_list}].copy()#*1e5      
            
            if inorm == 1:
                rc_att = 'att_bsc'
                plt_info = 'Atten. Volume Backscatter Signal' #'Attenuatted Volume Backscatter Signal'
                # cbar_unit = '$sr^{-1} m^{-1}$'    
            else:
                rc_att = 'rc'
                plt_info = 'RC Signal'
                # cbar_unit = '[A.U.]'    
            
            # Define the height scale
            if isig_hscale =='altitude':
                ytype = 'altitude a.s.l.'
                alt = sig.altitude.values/1e3                    
            else:
                ytype = 'range'
                alt = sig.range.values/1e3        
                
            cbar_unit = '[A.U.]'
                            
            # altitude/range axis settings
            ymin = float(options.qcklk.min_alt[0])/1e3# m to km
            ymax = float(options.qcklk.max_alt[0])/1e3 # m to km
            ystep = float(options.qcklk.step_alt[0])/1e3 # m to km
            alt_tick = np.arange(ymin, ymax+ystep, ystep)
            h_lim = f'{int(ymax)}km'
    
            
            for i in range(len(ch_list)): 
            
                ch = ch_list[i]
                ind = dict(channel=ch)
                #colorbar settings
                vllim = float(min_sig[i])
                vulim = float(max_sig[i])
                vstep = float(step_sig[i])
                vlevel = np.arange(vllim, vulim+vstep, vstep)
                
                
                # signal info (total/co-pol/cross-pol for lin/cir)
                #laser polarization
                if info.laser_pol.loc[ch] == 1.0:
                    sig_info = 'Linear'
                else:
                    sig_info = 'Circular'
    
                # channel type (total, co-/cross-polar)
                if info.ch_pol.loc[ch] in ['p','r']:
                    sig_info = ' '.join((sig_info,'co-pol'))
                if info.ch_pol.loc[ch] in ['s','l']:
                    sig_info = ' '.join((sig_info,'cross-pol'))
    
                plt_info_c = ' '.join((sig_info, plt_info))
    
         
                # Title
                plt_title = make_title(date, t1, t2, sza, lidar, loc_info)
           
                #Create the figure
                fig, ax = plt.subplots(figsize=(x_inches,y_inches))
            
                # Set the axes limits and xaxis format
                ax.axis([date_num.min(), date_num.max(), ymin, ymax])
                plt.xticks(rotation=30)
                ax.xaxis.set_major_locator(minutes)
                ax.xaxis.set_major_formatter(t_fmt)

                # 2-D plot                
                if scale_sig == 'log': # sig colorbar logarithmic scale 
                    qlook = ax.pcolormesh(date_num, alt, sig_plt.loc[ind].transpose().values,
                                          cmap=colormap_style,
                                          norm=mpl.colors.LogNorm(vmin=vllim, vmax=vulim))
                    cbar = fig.colorbar(qlook, pad=0.01)
                    scale_type = '_log'
                
                else: # sig colorbar linear scale 
                    qlook = ax.pcolormesh(date_num, alt, sig_plt.loc[:,ch,:].transpose().values,
                                          cmap=colormap_style, vmin=vllim, vmax=vulim)
                
                    cbar = fig.colorbar(qlook, pad=0.01, ticks=vlevel) 
                    cbar.ax.ticklabel_format(style='sci', scilimits=(0, 0))    
                    scale_type=''
                    
                # apply colorbar settings
                cbar.set_label(f'{plt_info_c} {cbar_unit}')
    
                # axes settings
                ax.set(ylabel=f'{ytype.capitalize()} [km]', yticks= alt_tick, 
                       xlabel='Time [UTC]', xlim=[date_num[0],date_num[-1]])
                ax.set_title(plt_title, pad=15.)
                ax.minorticks_on()
            
                if itime == 1: # plot vertical lines to indicate the averaged timeframe for retrieval
                    ax.axvline(x=tstart, ymin=0, ymax=1, ls='--', lw='2.', c='k')
                    ax.axvline(x=tstop, ymin=0, ymax=1, ls='--', lw='2.', c='k')
    
                fig_name = f'{lidar.lower()}_quicklook{scale_type}_{date_out}_{t1_out}_{t2_out}_{ytype[:8]}_{rc_att}_{ch}_{h_lim}{ver}.png'
                fig.savefig(os.path.join(out_dir, fig_name), 
                            bbox_inches = 'tight', dpi = dpi_val)
                plt.close()    
            
            print('-- Quicklook plots of signals complete!')
            
    return()

def prod_quicklook(cfg, prod, info_prod, options, ver=''):

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        colormap_style = select_colormap(options.qcklk) # make_colormap.cmap # 'Spectral_r' # 

        ql_dir = cfg.directories.plots
        iprod = int(options.qcklk.iprod[0])    
    
        if iprod == 1 and len(prod) > 0:
            # Create dir if not exist
            # out_dir = os.path.join(ql_dir, 'products')        
            out_dir = os.path.join(ql_dir,'quicklooks')            
            os.makedirs(out_dir, exist_ok=True)
    
            # Extract info from config file
            sza = np.round(float(cfg.angles.zenith), decimals=1)
            location = cfg.lidar.location
            lat = np.round(float(cfg.lidar.latitude), decimals=1)
            lon = np.round(float(cfg.lidar.longitude), decimals=1)
            elev = int(np.round(float(cfg.lidar.altitude), decimals=0))
            lidar = cfg.lidar.lr_name
            iprod_hscale = cfg.iprod_hscale
            
            loc_info = f'{location} (lat: {lat}, lon: {lon}, elev: {elev} m)'
            # Extract plotting options
            dpi_val = int(options.qcklk.dpi[0])
            itime = int(options.qcklk.itime)
                
            # prod ids for plotting
            vis_prod = options.qlk_vis_prod
            # settings for product colorbar        
            min_prod = options.qlk_min_prod
            max_prod = options.qlk_max_prod
            step_prod = options.qlk_step_prod
                
            # datetime title settings
            if itime == 1:
                tstart = mdates.date2num(np.datetime64(options.qcklk.tstart[0]))
                tstop = mdates.date2num(np.datetime64(options.qcklk.tstop[0]))
    
            date = prod.time[0].dt.strftime('%Y/%m/%d').values 
            date_out = prod.time[0].dt.strftime('%Y%m%d').values
            t1 = prod.time[0].dt.strftime('%H:%M').values
            t2 = prod.time[-1].dt.strftime('%H:%M').values
            t1_out = prod.time[0].dt.strftime('%H%M').values
            t2_out = prod.time[-1].dt.strftime('%H%M').values
            
            # convert meas timeframe to numbers (easier to plot)
            date_num = mdates.date2num(prod.time.values)
            
            # formating the time in x label
            tstep = int(options.qcklk.step_time[0])
            minutes = mdates.MinuteLocator(interval = tstep)#4 #20
            t_fmt = mdates.DateFormatter('%H:%M')
    
            # Define the height scale
            if iprod_hscale =='altitude':
                ytype = 'altitude a.s.l.'
                alt = prod.altitude.values/1e3                    
            else:
                ytype = 'range'
                alt = prod.range.values/1e3        
            
            # altitude/range axis settings
            ymin = float(options.qcklk.min_alt[0])/1e3# m to km
            ymax = float(options.qcklk.max_alt[0])/1e3 # m to km
            ystep = float(options.qcklk.step_alt[0])/1e3 # m to km
            alt_tick = np.arange(ymin, ymax+ystep, ystep)
            h_lim = f'{int(ymax)}km'
            
            # Select the desired products for plotting
            p_list = vis_prod
            prod_plt = prod.loc[{'product':p_list}].copy()
       
            for j in range(len(p_list)): 
            
                p = p_list[j]
                ind = dict(product = p)
                p_type = re.sub('_',' ',info_prod.type[p].title())
        
                #colorbar settings
                vllim = float(min_prod[j])
                vulim = float(max_prod[j])
                vstep = float(step_prod[j]) 
                vlevel = np.arange(vllim, vulim+vstep, vstep)
                
                # Title
                plt_title = make_title(date, t1, t2, sza, lidar, loc_info)
    
                #Create the figure            
                fig, ax = plt.subplots(figsize=(x_inches,y_inches))    
            
                # Set the axes limits and xaxis format
                ax.axis([date_num.min(), date_num.max(), ymin, ymax])
                plt.xticks(rotation=30)
                ax.xaxis.set_major_locator(minutes)
                ax.xaxis.set_major_formatter(t_fmt)
                                
                # 2-D plot
                qlook = ax.pcolormesh(date_num, alt, prod_plt.loc[ind].transpose().values,
                                      cmap=colormap_style, vmin=vllim, vmax=vulim)
                
                # apply colorbar settings
                cbar = fig.colorbar(qlook, pad=0.01, ticks=vlevel) 
                cbar.set_label(f'{p_type} {info_prod.units[p]}')
                if vulim <= 1e-3 or vulim >= 1e3:
                    cbar.ax.ticklabel_format(style='sci', scilimits=(0, 0))
            
                # axes settings
                ax.set(ylabel=f'{ytype.capitalize()} [km]', yticks= alt_tick, 
                       xlabel='Time [UTC]', xlim=[date_num[0],date_num[-1]])
                ax.set_title(plt_title, pad=15.)
    
                ax.minorticks_on()
            
                if itime ==1:# plot vertical lines to indicate the averaged timeframe for retrieval
                    ax.axvline(x=tstart, ymin=0, ymax=1, ls='--', lw='2.', c='k')
                    ax.axvline(x=tstop, ymin=0, ymax=1, ls='--', lw='2.', c='k')
                    
                fig_name = f'{lidar.lower()}_quicklook_{date_out}_{t1_out}_{t2_out}_{ytype[:8]}_{info_prod.type_short[p]}_{h_lim}{ver}.png'
                fig.savefig(os.path.join(out_dir, fig_name), 
                            bbox_inches = 'tight', dpi = dpi_val)
                plt.close()
        
            print('-- Quicklook plots of products complete!')
    
    return()

def qclk_screenout(sig_r, num):
    
    plt.close()
    fig, ax = plt.subplots(figsize=(10,8))
    sig_r[0,num-10:num,4,:].plot(hue='time', x='range',xlim=[0,3e3],ylim=[0,2e6],
                                ax=ax)    

    return()

def make_title(date, t1, t2, sza, lidar, loc_info):
    
    plot_title = 'Time-Height cross sections\n'+\
                f'On {date} from {t1} to {t2} UTC, '+r'$\nearrow$'+\
                    f' {sza}'+r'$^{o}$ off-zenith' +\
                        f'\n{lidar} at {loc_info}'
    
    return(plot_title)

def select_colormap(qcklc):
    
    # check if attr exists in plot_options.ini
    qcklc = check_class_attr(qcklc,'colormap_style','Spectral_r')
    
    if 'Spectral_r' in qcklc.colormap_style[0]:
        color_map_style = 'Spectral_r'
    elif 'cmap'== qcklc.colormap_style[0]:
        color_map_style = make_colormap.cmap    
    elif 'cmap2'== qcklc.colormap_style[0]:
        color_map_style = make_colormap.cmap2
    elif 'cmap3'== qcklc.colormap_style[0]:
        color_map_style = make_colormap.cmap3
    elif 'cmap4'== qcklc.colormap_style[0]:
        color_map_style = make_colormap.cmap4

    return(color_map_style)