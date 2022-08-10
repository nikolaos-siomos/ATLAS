"""
@author: N. Siomos & P. Paschou 

Plooting routine for the telecover test
"""
from scipy import stats
import math
import matplotlib as mpl
from matplotlib import pyplot as plt
from pathlib import Path
from matplotlib import gridspec
from PIL import Image
import numpy as np
from helper_functions.helper_functions import datetime_unit
import os
import sys

ln_wdth = 2.
ax_lblsz = 35
xtck_lblsz = ytck_lblsz =28
ax_lblpd = 30
lgd_fsz = 20

xtick_major_size = 12
xtick_minor_size = 8
text_size = 35


def plot_rfit_one(dir_plt, sig_nrm, sig_mol, alt_lims, cfg, plt_opt):

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         
    
        lidar      = cfg.lidar
        timescale  = cfg.timescale
        sza        = np.round(float(cfg.angles.zenith), decimals=1)
        #info from plt options
        dpi_val    = int(plt_opt.dpi[0])
        subplt_num = plt_opt.subplt_num[0]
        scale      = str(plt_opt.scale[0])
        min_alt    = float(plt_opt.min_alt_ray[0])
        max_alt    = float(plt_opt.max_alt_ray[0])
        step_alt   = float(plt_opt.step_alt_ray[0])
    
        # channels, time, altitude arrays    
        num_ch = sig_nrm.channel.size
        time = sig_nrm.time.values
        alt = sig_nrm.altitude.values
        channel = sig_nrm.channel.values
        
        # For empty subplt_num, create one figure with all subplots
        if subplt_num != subplt_num:
            subplt_num = num_ch
        else:
            subplt_num = int(subplt_num)
        
        if max_alt > alt[-1]:
            max_alt = np.round(alt[-1],decimals=-3)
            
        # alt_ticks = np.arange(0, alt[-1]+1e3, 1e3)
        alt_ticks = np.arange(min_alt, max_alt+step_alt, step_alt)    
    
        for t in range(time.size):
            ind_t = dict(time = t)        
            # Create the datetime strings
            dt_str = sig_nrm[ind_t].time.dt.strftime('%Y%m%dT%H%M%SZ').values        
            date = sig_nrm[ind_t].time.dt.strftime('%Y-%m-%d').values
            #start time
            s_time = sig_nrm[ind_t].time.dt.strftime('%H:%M').values
            #end_time
            t_unit = datetime_unit(timescale).lower() # time unit e.g. s, m, h, d 
            number=int(''.join(filter(str.isdigit, timescale)))        
            e_time = (sig_nrm[ind_t].time + np.timedelta64(number, t_unit)).dt.strftime('%H:%M').values
            
            # Create title name of figure
            title_name = f'Time frame: {date}, {s_time} - {e_time} UTC\n'+\
                f'Measurement angle: {sza}'+'$^{o}$ off-zenith\n'+\
                    f'{lidar.lr_name} @ {lidar.location}'
            
            # Create the full path with the plot name
            export_path = os.path.join(dir_plt, f'rayleigh_fit_{lidar.lr_name}_{dt_str}_{scale}_(1).png')
            
            # Create temp. directory for saving the subplots plots per timeframe
            temp_dir = Path('./tmp')
            temp_dir.mkdir(exist_ok=True)
            
            # Number of sub-figures with -subplt_num- channels plotted per sub-figure                 
            figures_num = int(np.ceil(num_ch / subplt_num))
            modulo = num_ch % subplt_num
            
            # Loop to create the sub-figures and distribute accordingly the...
            # ... vis_map (channels) in the sub-figures
            for sub_fig in range(figures_num):
    
                start = sub_fig * subplt_num
                stop = (sub_fig+1) * subplt_num
                
                # if it's in the last sub figure, end loop in the last channel
                if (sub_fig == figures_num-1) and (modulo != 0):
                    stop = start + modulo
                
                if figures_num != 1:
                    # Modify the full path with the (sub)figure name
                    export_path = os.path.join(dir_plt, f'rayleigh_fit_{dt_str}_{scale}_(1)_part{sub_fig}.png')
                
                #Loop to make the channel subplots in the (sub)figure
                for j in range(start, stop):
                    
                    # Create the figure and the axis
                    fig_size = (30,6)
                    
                    if j==start: # make slightly bigger the first plot to fit the suptitle
                        fig_size = (30,7.5)
                    fig = plt.figure(figsize=fig_size)
                    spec = gridspec.GridSpec(ncols=1, nrows=1)
                    
                    ax = fig.add_subplot(spec[0, 0])
                    fig.subplots_adjust(left=0.1, right=0.95)
                    
                    # Select the normal and mol signal
                    nrm = sig_nrm.sel(altitude=slice(0,max_alt)).isel(channel = j, time = t).values
                    mol = sig_mol.sel(altitude=slice(0,max_alt)).isel(channel = j).values
                    alt=alt[:mol.size]
                    # Select the corresponding lims indicating the normalization region
                    llim = alt_lims.isel(channel = j).sel(lims = 'llim')
                    ulim = alt_lims.isel(channel = j).sel(lims = 'ulim')
                    
                    #rfit plot
                    ax.plot(alt, nrm, color='k')
                    ax.plot(alt, mol, color='r')
                    ax.axvline(llim, color='b', linestyle='--')
                    ax.axvline(ulim, color='b', linestyle='--')
                    ax.axvspan(llim, ulim, facecolor = 'dodgerblue', alpha = 0.5)
    
                                                    # '\n Att. Bsc [$sr^{-1} m^{-1}$]',        
                    ax.set(ylabel = channel[j]+'\n Att. Bsc Signal [A.U]', 
                           yscale = scale, 
                           ylim = (min(mol)/5., max(mol)*20.),
                           xlim = (min_alt, max_alt), xticks = alt_ticks,
                           xticklabels = (alt_ticks/1e3).astype(int).astype(str))
    
                    ax.tick_params(axis='both', which='major', size=xtick_major_size)
                    ax.tick_params(axis='both', which='minor', size=xtick_minor_size)
                    
                    ax.minorticks_on()
                    
                    # Set the x label only in the last plot
                    if j == stop - 1:
                        ax.set_xlabel('Altitude [km]')
        
                    # set the figure title in the first subplot
                    if j == start:
                        # create space for the suptitle
                        # ax.set_title(' \n \n \n', fontsize = 35) 
                        fig.suptitle(title_name, fontsize = 35)
    
                    # Save the fig
                    spec.tight_layout(fig)            
                    fig.savefig(temp_dir / f"{j}ch{sub_fig}.png", 
                                bbox_inches = 'tight', dpi = dpi_val)
                    plt.close(fig)
                
                # stitch all figures (per channel) to one image
                stitch_subplots_to_image(temp_dir, export_path,f'*{sub_fig}.png')
    
            temp_dir.rmdir()
            
    return()


def plot_rfit_two(dir_plt, sig_nrm, sig_mol, alt_lims, cfg, plt_opt):
    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         
    
        lidar      = cfg.lidar
        timescale  = cfg.timescale 
        sza        = np.round(float(cfg.angles.zenith), decimals=1)
    
        dpi_val    = int(plt_opt.dpi[0])
        subplt_num = plt_opt.subplt_num[0]
        scale      = plt_opt.scale[0]
        #axes info for rfit
        min_alt_ray    = float(plt_opt.min_alt_ray[0])
        max_alt_ray    = float(plt_opt.max_alt_ray[0])
        step_alt_ray   = float(plt_opt.step_alt_ray[0])
        #axes info for rel dif
        min_alt_dif    = float(plt_opt.min_alt_dif[0])
        max_alt_dif    = float(plt_opt.max_alt_dif[0])
        step_alt_dif   = float(plt_opt.step_alt_dif[0])
        min_rel_dif    = float(plt_opt.min_rel_dif[0])
        max_rel_dif    = float(plt_opt.max_rel_dif[0])
        step_rel_dif   = float(plt_opt.step_rel_dif[0])
        
        # channels, time, altitude arrays
        num_ch = sig_nrm.channel.size
        time = sig_nrm.time.values
        alt = sig_nrm.altitude.values
        channel = sig_nrm.channel.values
    
        # For empty subplt_num, create one figure with all subplots
        if subplt_num != subplt_num:
            subplt_num = num_ch
        else:
            subplt_num = int(subplt_num)
        
        max_alt = max(max_alt_ray, max_alt_dif)
        if max_alt > alt[-1]:
            max_alt = np.round(alt[-1],decimals=-3)
            # max_alt_ray = np.round(alt[-1],decimals=-3)
            # max_alt_dif = np.round(alt[-1],decimals=-3)
            
        alt_ticks_ray = np.arange(min_alt_ray, max_alt_ray+step_alt_ray, step_alt_ray)
        alt_ticks_dif = np.arange(min_alt_dif, max_alt_dif+step_alt_dif, step_alt_dif)
        dif_ticks = np.arange(min_rel_dif, max_rel_dif+step_rel_dif, step_rel_dif)
        
        for t in range(time.size):
            ind_t = dict(time = t)
            # Create the datetime strings
            dt_str = sig_nrm[ind_t].time.dt.strftime('%Y%m%dT%H%M%SZ').values        
            date = sig_nrm[ind_t].time.dt.strftime('%Y-%m-%d').values
            #start time
            s_time = sig_nrm[ind_t].time.dt.strftime('%H:%M').values
            #end_time
            t_unit = datetime_unit(timescale).lower() # time unit e.g. s, m, h, d 
            number=int(''.join(filter(str.isdigit, timescale)))        
            e_time = (sig_nrm[ind_t].time + np.timedelta64(number, t_unit)).dt.strftime('%H:%M').values
            
            # Create title name of figure
            title_name = f'Time frame: {date}, {s_time} - {e_time} UTC \n'+\
                f'Measurement angle: {sza}'+'$^{o}$ off-zenith \n'+\
                    f'{lidar.lr_name} @ {lidar.location}'
            
            # Create the full path with the plot name
            export_path = os.path.join(dir_plt, f'rayleigh_fit_{lidar.lr_name}_{dt_str}_{scale}_(2).png')
            
            # Create temp. directory for saving the subplots plots per timeframe
            temp_dir = Path('./tmp')
            temp_dir.mkdir(exist_ok=True)
            
            # Number of sub-figures with -subplt_num- channels plotted per sub-figure                 
            figures_num = int(np.ceil(num_ch / subplt_num))
            modulo = num_ch % subplt_num
            
            # Loop to create the sub-figures and distribute accordingly the...
            # ... vis_map (channels) in the sub-figures
            for sub_fig in range(figures_num):
    
                start = sub_fig * subplt_num
                stop = (sub_fig+1) * subplt_num
                
                # if it's in the last sub figure, end loop in the last channel
                if (sub_fig == figures_num-1) and (modulo != 0):
                    stop = start + modulo
                
                if figures_num != 1:
                    # Modify the full path with the (sub)figure name
                    export_path = os.path.join(dir_plt, f'rayleigh_fit_{dt_str}_{scale}_(2)_part{sub_fig}.png')
    
                #Loop to make the channel subplots in the (sub)figure              
                for j in range(start, stop):
                    
                    # Create the figure and the axis
                    fig_size = (30,6)
                    
                    if j == start: # make slightly bigger the first plot to fit the suptitle
                        fig_size = (30,7.5) # 7 an exw 2 seires title
                    fig = plt.figure(figsize=fig_size)
                    spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[8, 4])
                    
                    # for 2 subplots (rayleigh + relative_diff)
                    ax_ray = fig.add_subplot(spec[0, 0])
                    ax_dif = fig.add_subplot(spec[0, 1])
                    fig.subplots_adjust(left=0.1, right=0.95)
    
                    # Select the corresponding lims indicating the normalization region
                    llim = alt_lims.isel(channel = j).sel(lims = 'llim').values
                    ulim = alt_lims.isel(channel = j).sel(lims = 'ulim').values
                    
                    # Select the normalised and molec signal
                    nrm = sig_nrm.sel(altitude=slice(0,max_alt)).isel(channel = j,time = t).values
                    mol = sig_mol.sel(altitude=slice(0,max_alt)).isel(channel = j).values
                    alt=alt[:mol.size]
                    alt_refmsk = (alt >= llim) & (alt <= ulim)
    
                    #Relative difference 
                    dif = (nrm - mol)
                    rel_dif = dif/mol
                    magnitude = math.floor(math.log(np.nanmax(mol[alt_refmsk]),10))
    
                    # standard error of the mean (sem) of differences inside the ref_region
                    # stdev = np.std(rel_dif[alt_refmsk]) # swsto me ddof =1
                    # sem = stdev/np.sqrt(rel_dif[alt_refmsk].size)
                    sem = stats.sem(dif[alt_refmsk], nan_policy='omit')/(10**magnitude)
                    
                    #rfit plot
                    ax_ray.plot(alt, nrm, color='k')
                    ax_ray.plot(alt, mol, color='r')
                    ax_ray.axvline(llim, color='b', linestyle='--')
                    ax_ray.axvline(ulim, color='b', linestyle='--')
                    ax_ray.axvspan(llim, ulim, facecolor = 'dodgerblue', alpha = 0.5)
                    
                                                    # '\n Att. Bsc [$sr^{-1} m^{-1}$]',    
                    ax_ray.set(ylabel = channel[j]+'\n Att. Bsc Signal [A.U.]', 
                              yscale = scale, 
                              ylim = (min(mol)/5., max(mol)*20.),
                              xlim = (min_alt_ray, max_alt_ray), xticks = alt_ticks_ray,
                              xticklabels = (alt_ticks_ray/1e3).astype(int).astype(str))
                    
                    ax_ray.tick_params(axis='both', which='major', size=xtick_major_size)
                    ax_ray.tick_params(axis='both', which='minor', size=xtick_minor_size)
                    
                    ax_ray.minorticks_on()
                    
                    # rel dif plot
                    ax_dif.plot(alt, rel_dif, color='k')
                    ax_dif.axhline(y=0,xmin=0, xmax=1, color='lime')
                    ax_dif.axvline(llim, color='b', linestyle='--')
                    ax_dif.axvline(ulim, color='b', linestyle='--')
                    ax_dif.axvspan(llim, ulim, facecolor = 'dodgerblue', alpha = 0.5)
                    # max_rel_dif-step_rel_dif/2
                    # ax_dif.annotate(f'sem = {np.round(sem, 4)}', (ulim+50.,max_rel_dif-step_rel_dif/2),
                    #                 fontsize = text_size, facecolor = 'b')
                    ax_dif.text(ulim*1.02,max_rel_dif-step_rel_dif/2,f'sem = {np.round(sem, 4)}', 
                                fontsize = text_size, color = 'blue',# backgroundcolor='white',
                                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
                    ax_dif.set(ylabel = 'Relative diff.',
                               ylim = (min_rel_dif, max_rel_dif), yticks = dif_ticks,
                               xlim = (min_alt_dif, max_alt_dif), xticks = alt_ticks_dif,
                               xticklabels = (alt_ticks_dif/1e3).astype(int).astype(str))
    
                    ax_dif.tick_params(axis='both', which='major', size=xtick_major_size)
                    ax_dif.tick_params(axis='both', which='minor', size=xtick_minor_size)
    
                    ax_dif.minorticks_on()
                                
        
                    # Set the x label(s) only in the last subplot
                    if j == stop - 1:
                        ax_ray.set_xlabel('Altitude [km]')
                        ax_dif.set_xlabel('Altitude [km]')                
        
                    # set the figure title in the first subplot
                    if j == start:
                        # create space for the suptitle
                        # ax_ray.set_title(' \n \n \n', fontsize = 35) 
                        fig.suptitle(title_name, fontsize = 35)             
    
                    # Save the fig
                    spec.tight_layout(fig)            
                    fig.savefig(temp_dir / f"{j}ch{sub_fig}.png", 
                                bbox_inches = 'tight', dpi = dpi_val)
                    plt.close(fig)
    
                # stitch all figures (per channel) to one image
                stitch_subplots_to_image(temp_dir, export_path,f'*{sub_fig}.png')
    
            temp_dir.rmdir()
    return()
            

def stitch_subplots_to_image(temp_dir, export_path, search_type='*.png'):

    # path list of the created subfigs in temp dir sorted by time
    f_list = sorted(temp_dir.glob(search_type), key=os.path.getctime)

    # Stitch PNGs together
    plots = [Image.open(name) for name in f_list]
    widths, heights = zip(*(img.size for img in plots))
    
    # Target dimension
    dims = (max(widths), sum(heights))
    stitched_plot = Image.new('RGB', dims)
    
    # Paste plots into new image
    offset = 0
    for plot in plots:
        stitched_plot.paste(plot, (0, offset))
        offset += plot.size[1]
    
    # Write final plot, delete temp files
    stitched_plot.save(export_path)
    
    for temp_image in temp_dir.glob(search_type):
        temp_image.unlink()

    return()