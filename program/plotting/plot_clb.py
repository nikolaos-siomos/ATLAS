"""
Plotting for calibration routine and calibration tests

@author: N. Siomos & P. Paschou
"""
import numpy as np
import os, re
from os.path import join as join
import xarray as xr
from molecular.height_scales import range_to_alt
from lidar_processing.construct import extract_error
from readers.read_maps import plt_options
from matplotlib import pyplot as plt
import matplotlib as mpl

ln_wdth = 2.
ax_lblsz = 15
ax_tlsz = 15
xtck_lblsz = 15
ytck_lblsz = 15
ax_lblpd = 15.0
lgd_fsz= 15

text_font = 14

def plot_d90(rt_p, rt_m, rt_pm, cal_fac, angle, range_d90, c_mode, cfg, maps_dir):
# def plot_d90(rt_p, rt_m, rt_pm, cal_fac, angle, cfg, l_r, u_r, dpi_val, ishow, label):
# plots the rt_p45, rt_m45 and the eta_d90{geometric mean of rt_p45*rt_m45}

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        iplt = cfg.iplt
        dirs = cfg.directories
        sza  = np.round(float(cfg.angles.zenith), decimals=1)
    
        # Extract plotting options from plot_options.ini
        options = plt_options(maps_dir, 'plot_options.ini','Calibration').clb
        dpi_val = int(options.dpi[0])
        label = options.d90_type[0]
        ishow = int(options.ishow[0])
        l_alt = float(options.min_alt[0])
        u_alt = float(options.max_alt[0])
        step = float(options.step_alt[0])
        
        alt_ticks = np.arange(l_alt, u_alt+step, step) 
        # rt_p, rt_m, rt_pm, cal_fac, angle = delta_90(cfg, clb, sig_p, sig_m)
    
        # range limits where the calibration factor was calculated 
        #-> convert from m to km
        # llim_cf = range_to_alt(range_d90[0], 0., sza) # range_d90[0] 
        # ulim_cf = range_to_alt(range_d90[1], 0., sza) # range_d90[1]

        llim_cf = range_to_alt(range_d90[0], float(cfg.lidar.altitude), sza) # range_d90[0] 
        ulim_cf = range_to_alt(range_d90[1], float(cfg.lidar.altitude), sza) # range_d90[1]

        
        if iplt in [2,3]:
    
            print('-----------------------------------------')
            print('Generating calibration factor plots!')
            print('-----------------------------------------')        
    
            path = join(dirs.plots, 'calibration', 'd90', label)
            os.makedirs(path, exist_ok=True)
            
            # Extract error Remove the iters dimension
            rt_p, _ = extract_error(rt_p)
            rt_m, _ = extract_error(rt_m)
            rt_pm, _ = extract_error(rt_pm)
            cal_fac, cal_fac_err  = extract_error(cal_fac)
            angle, angle_err = extract_error(angle)
        
            for k in range(rt_p.time.size):
                t = str(rt_p.time.values[k])[0:19]
        
                for j in range(rt_p.channel.size):
                    
                    sel_t_ch = dict(time=k, channel=j)
                    
                    if not all(np.isnan(rt_p[sel_t_ch])):
                        
                        cal_f = np.round(cal_fac[sel_t_ch].values, decimals = 4)
                        cal_f_err = np.round(cal_fac_err[sel_t_ch].values, decimals = 6)
                        
                        offset = np.round(angle[sel_t_ch].values, decimals = 2)
                        offset_err = np.round(angle_err[sel_t_ch].values, decimals = 4)
                        
                        fig, ax = plt.subplots(figsize=(10,8))
                        rt = rt_p.channel.values[j]
            
                        ax.plot(rt_p[sel_t_ch].values, rt_p.alt.values)
                        ax.plot(rt_m[sel_t_ch].values, rt_m.alt.values)
                        ax.plot(rt_pm[sel_t_ch].values, rt_pm.alt.values)
                        
                        # plot a faded region where the calibr. factor was calculated
                        ax.axhspan(ymin=llim_cf, ymax=ulim_cf,
                                   facecolor='grey', alpha=0.5, zorder=3)
                        
                        rt_lim = 2. * cal_f #2.5
    
                        title_name = f'Setup: {rt}, Frame: {t}\n'+\
                            f'Channel mode: {c_mode}\n'+\
                            f'Measurement angle: {sza}'+'$^{o}$ off-zenith'
                    
                        ax.set(title = title_name,
                               xlim = [0, rt_lim], xlabel='Calibration Factor',
                               ylim = [l_alt, u_alt], ylabel='Altitude [km]',
                               yticks = alt_ticks,
                               yticklabels = (alt_ticks/1e3).astype(str))
                        
                        ax.annotate(f'cal. factor: {cal_f} '+r'$\pm$'+f' {cal_f_err}', 
                                    (0.03*rt_lim, 0.1*(u_alt-l_alt)),
                                    fontsize=text_font, backgroundcolor='white')
                        ax.annotate(f'offs. angle: {offset} '+r'$\pm$'+f' {offset_err}', 
                                    (0.03*rt_lim, 0.05*(u_alt-l_alt)),
                                    fontsize=text_font, backgroundcolor='white')
                        
                        ax.minorticks_on()
                        
                        ax.legend(['plus','minus','gmean'], loc=1)
    
                        fig.savefig(join(path, f'{rt}_Frame_{k}_{c_mode}.png'), dpi = dpi_val)
                        
                        if ishow == 1:
                            plt.show()
                        else:
                            plt.close()
    
    return()


def plot_d90_test_QWP(rt_p, rt_m, rt_pm, cal_fac, angle, c_mode, cfg, maps_dir):
    
# plots the rt_p45, rt_m45 and the eta_d90{geometric mean of rt_p45*rt_m45}
# Plots also the offset angle of the QWP_LB 
    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        iplt = cfg.iplt
        dirs = cfg.directories
        sza  = np.round(float(cfg.angles.zenith), decimals=1)
    
    
        # Extract plotting options from plot_options.ini
        options = plt_options(maps_dir, 'plot_options.ini','Calibration').clb
        dpi_val = int(options.dpi[0])
        label = options.d90_type[0]
        ishow = int(options.ishow[0])
        l_alt = float(options.min_alt[0])
        u_alt = float(options.max_alt[0])
        step = float(options.step_alt[0])
        
        alt_ticks = np.arange(l_alt, u_alt+step, step) 
    
        
        if iplt in [2,3]:
    
            print('-----------------------------------------')
            print('Generating calibration factor plots!')
            print('-----------------------------------------')        
    
            path = join(dirs.plots, 'calibration', 'd90', label)
            os.makedirs(path, exist_ok=True)
            
            # Extract error Remove the iters dimension
            rt_p, _ = extract_error(rt_p)
            rt_m, _ = extract_error(rt_m)
            rt_pm, _ = extract_error(rt_pm)
            cal_fac, cal_fac_err  = extract_error(cal_fac)
            angle, angle_err = extract_error(angle)
    
            for k in range(rt_p.time.size):
                t = str(rt_p.time.values[k])[0:19]
        
                for j in range(rt_p.channel.size):
                    
                    sel_t_ch = dict(time=k, channel=j)
                    
                    if not all(np.isnan(rt_p[sel_t_ch])):
                        
                        cal_f = np.round(cal_fac[sel_t_ch].values, decimals = 4)
                        cal_f_err = np.round(cal_fac_err[sel_t_ch].values, decimals = 6)
                        
                        offset = np.round(angle[sel_t_ch].values, decimals = 2)
                        offset_err = np.round(angle_err[sel_t_ch].values, decimals = 4)
                        
                        fig, ax = plt.subplots(figsize=(10,8))
                        rt = rt_p.channel.values[j]
            
                        ax.plot(rt_p[sel_t_ch].values, rt_p.alt.values)
                        ax.plot(rt_m[sel_t_ch].values, rt_m.alt.values)
                        ax.plot(rt_pm[sel_t_ch].values, rt_pm.alt.values)
    
                        
                        # rt_lim = 1.5*np.nanmax(rt_pm[sel_t_ch].loc[dict(alt=slice(l_alt+0.5e3,u_alt-0.5e3))].values)
                        rt_lim = 2. * cal_f #2.5
    
    
                        title_name = f'Setup: {rt}, Frame: {t}\n'+\
                            f'Channel mode: {c_mode}\n'+\
                            f'Measurement angle: {sza}'+'$^{o}$ off-zenith'
                    
                        ax.set(title = title_name,
                                xlim = [0, rt_lim], xlabel='Calibration Factor',
                                ylim = [l_alt, u_alt], ylabel='Altitude [km]',
                                yticks = alt_ticks,
                                yticklabels = (alt_ticks/1e3).astype(str))
                        
                        ax.annotate(f'cal. factor: {cal_f} '+r'$\pm$'+f' {cal_f_err}', 
                                      (0.03*rt_lim, 0.1*(u_alt-l_alt)), #(0.05*rt_lim, 0.95*(u_alt-l_alt))
                                      fontsize=text_font, backgroundcolor='white')
                        ax.annotate(f'offs. angle QWP_LB: {offset} '+r'$\pm$'+f' {offset_err}', 
                                      (0.03*rt_lim, 0.05*(u_alt-l_alt)), #(0.05*rt_lim, 0.85*(u_alt-l_alt))
                                      fontsize=text_font, backgroundcolor='white')
                        
                        ax.legend(['plus','minus','gmean'], loc=1)
    
                        ax.minorticks_on()                                                                    
    
                        fig.savefig(join(path, f'{rt}_Frame_{k}_{c_mode}.png'), dpi = dpi_val)                    
                        if ishow == 1:
                            plt.show()
                        else:
                            plt.close()

    return()


def plot_rt(sig_meas, cfg, clb, l_r, u_r, dpi_val, ishow, label):
# calculates the ratio R/T, corrects it with the calib_f and plots it

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        iplt = cfg.iplt    
        rt_id     = clb.rt_id_rt
        ch_up     = clb.ch_up_rt
        ch_dn     = clb.ch_dn_rt
        cal_fact  = clb.cal_fact
        
        eta = xr.DataArray(cal_fact, coords=[rt_id], dims=['channel'])
        
        alt = range_to_alt(sig_meas.range.values, 
                           float(cfg.lidar.altitude), cfg.angles.zenith)
        
        ratio = sig_meas.loc[:,ch_up,:].values/sig_meas.loc[:,ch_dn,:].values
        rt_meas = xr.DataArray(ratio, 
                               coords=[sig_meas.time, rt_id, alt],
                               dims=['time', 'channel', 'alt'])
        
        # corrects the signals ratio (R/T) with the calibration factor
        rt_meas = rt_meas/eta
        
        print(rt_meas.loc[:,:,500:1500].mean(axis=2,skipna=True))
    
        for j in range(rt_meas.channel.size):   
            for k in range(rt_meas.time.size):
                ch = rt_meas.channel.values[j]
                t = str(rt_meas.time.values[k])[0:19]
                
                plt.title(f'Channel: {ch}, Frame: {t}' )
                plt.plot(rt_meas[k,j,:], rt_meas.alt)
                
                rt_lim = 1.3*np.nanmax(rt_meas[k,j,:].loc[l_r:u_r].values)
                # rt_lim = 0.4
                
                plt.ylim([l_r, u_r])
                plt.xlim([0, rt_lim])
                
                if iplt in [2,3]:
                    path = join(cfg.directories.plots, 'calibration', 'ratios', label)
                    if not os.path.exists(path):
                        os.makedirs(path)
     
                    #stime = str(t)[0:19]
                    dt_str = re.split('[- : .]', t)
                    dt_str = ''.join(dt_str[:-2] + [dt_str[-1][:2]])
            
                    plt.savefig(join(path,f'{ch}_{dt_str}.png'), dpi = dpi_val)
                    if ishow == 1:
                        plt.show()
                    else:
                        plt.close()
    
    return()


def plot_d90_tri(sig, sig_p, sig_m, cfg, clb, l_r, u_r, dpi_val, ishow, label):
    # Calculates the ratio R/T of 3 meas (norm, +45, -45) and plots the ratios

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         
    
        # rt_id = clb.rt_id_d90
        # ch_up = clb.ch_up_d90
        # ch_dn = clb.ch_dn_d90
        iplt = cfg.iplt
        rt_id = clb.rt_id_d90
        ch_up = clb.ch_up_d90
        ch_dn = clb.ch_dn_d90
        # clb_r = clb.range_d90
        
        alt = range_to_alt(sig_p.range.values, 
                           float(cfg.lidar.altitude), cfg.angles.zenith)
    
        ratio_p = sig_p.loc[:,ch_up,:].values/sig_p.loc[:,ch_dn,:].values    
        rt_p = xr.DataArray(ratio_p, 
                            coords=[sig_p.time, rt_id, alt],
                            dims=['time', 'channel', 'alt'])    
        
        ratio_m = sig_m.loc[:,ch_up,:].values/sig_m.loc[:,ch_dn,:].values    
        
        rt_m = xr.DataArray(ratio_m, 
                               coords=[sig_m.time, rt_id, alt],
                               dims=['time', 'channel', 'alt'])    
        
        ratio = sig.loc[:,ch_up,:].values/sig.loc[:,ch_dn,:].values
        
        rt = xr.DataArray(ratio, 
                          coords=[sig.time, rt_id, alt],
                          dims=['time', 'channel', 'alt'])        

        for j in range(rt_p.channel.size):   
            for k in range(rt_p.time.size):
                ch = rt_p.channel.values[j]
                t = str(rt_p.time.values[k])[0:19]
    
                plt.title(f'Channel: {ch}, Frame: {t}' )
                plt.plot(rt_p[k,j,:], rt_p.alt)
                plt.plot(rt_m[k,j,:], rt_m.alt)
                plt.plot(rt[k,j,:], rt.alt)
                
                # rt_lim = 1.3*np.nanmax([rt_p[k,j,:].loc[l_r:u_r].values,
                #                         rt_m[k,j,:].loc[l_r:u_r].values,
                #                         rt[k,j,:].loc[l_r:u_r].values])
                rt_lim = 0.25            
                # print(rt[k,j,1000:2500].mean(skipna = True))
                
                plt.ylim([l_r, u_r])
                plt.xlim([0, rt_lim])
                # plt.xticks(np.arange(0,1.1,0.1)) 
                
                plt.legend(['+45','-45','0'])
                
                if iplt in [2,3]:
                    path = join(cfg.directories.plots, 'calibration', 'd90_tri', label)
                    if not os.path.exists(path):
                        os.makedirs(path)
                    plt.savefig(join(path, f'{ch}_Frame_{k}.png'), dpi = dpi_val)
                    if ishow == 1:
                        plt.show()
                    else:
                        plt.close()
        
    return()