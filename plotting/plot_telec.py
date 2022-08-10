"""
@author: N. Siomos & P. Paschou 

Plooting routine for the telecover test
"""

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import os

# mpl.rcParams.update(mpl.rcParamsDefault)
ln_wdth = 2.5
ax_lblsz = 22 
ax_tlsz = 24
xtck_lblsz = ytck_lblsz = 20
lgd_fsz = 22
ax_lblpd = 18.0
xtck_mj_sz = ytck_mj_sz = 6
xtck_mn_sz = ytck_mn_sz = 4
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.minor.size'] = 4


def plot_telec(telc_l, sig_telec_bc, sig_norm_rc, dev_all, ndev_all, 
               dev_N2N, ndev_N2N, dev_EW, ndev_EW, 
               sectors, info, datetime, cfg, options):

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz,
                        'xtick.major.size':xtck_mj_sz, 'ytick.major.size':ytck_mj_sz,
                        'xtick.minor.size':xtck_mn_sz, 'ytick.minor.size': ytck_mn_sz}):         

        sza = np.round(float(cfg.angles.zenith), decimals=1)
    
        
        # Extract plotting options
        dpi_val = int(options.dpi[0])
        
        dev_lim = float(options.crit_dev[0])#0.05
        
        # Normalization options
        norm_range = float(options.norm_range[0]) # range in meters
        norm_hwin   = float(options.norm_hwin[0]) # range in meters
        # limits of range axis (meters) and ticks step
        lrange = float(options.min_range[0])
        urange = float(options.max_range[0]) 
        step_range = float(options.step_range[0])
        # limits of raw signal axis and ticks step    
        min_dev_sig = float(options.min_dev_sig[0])
        max_dev_sig = float(options.max_dev_sig[0])
        step_dev_sig = float(options.step_dev_sig[0])
        # limits of normalized signal axis and ticks step
        min_dev_nrm = float(options.min_dev_nrm[0])
        max_dev_nrm = float(options.max_dev_nrm[0])
        step_dev_nrm = float(options.step_dev_nrm[0])
    
        rg_ticks = np.arange(lrange, urange+step_range, step_range)
        sg_dev_ticks = np.arange(min_dev_sig, max_dev_sig+step_dev_sig, step_dev_sig)
        nrm_dev_ticks = np.arange(min_dev_nrm, max_dev_nrm+step_dev_nrm, step_dev_nrm)
    
        
        stime_t = datetime.dt.strftime('%Y-%m-%d, %H:%M:%S UTC').values #%Y-%m-%dT%H:%M:%SZ
        dt_str = datetime.dt.strftime('%Y%m%dT%H%M').values
        
        print('-----------------------------------------')
        print('Generating telecover plots...')
        print('-----------------------------------------')
    
        for ch in sig_telec_bc.channel.values:
            
            # Make the plots if the signal values are not NaN in one of the sectors (i.e. empty for the rest sectors too)
            if info.loc[ch,'ch_mode']== 0. and np.isnan(sig_telec_bc.loc['north',ch,:].values).all() == False:
            
                y_max_bc  = np.round(sig_telec_bc.loc[:,ch,:].loc[:,lrange:urange].max()*1.5)
                y_max_nrc = np.round(sig_norm_rc.loc[:,ch,:].loc[:,lrange:urange].max()*1.5)
                # y_max_nrc = 1.5
                
                fig,((ax11,ax12),(ax21,ax22),(ax31,ax32)) =\
                    plt.subplots(nrows=3,ncols=2,
                                 figsize = (20,16)) #16, 12
                # Fig row with signals
                #subplot 11
                sig_telec_bc.loc[:,ch,:].plot.line(hue='sector', ylim = [0, y_max_bc], ax=ax11)
                ax11.legend(sectors['sector'], loc = 'upper right')
                ax11.set(xlabel='Range [m]', xlim = [lrange, urange],
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         title='Raw Signal')
                ax11.minorticks_on()
    
                #subplot 12
                sig_norm_rc.loc[:,ch,:].plot.line(hue='sector', ylim = [0, y_max_nrc], ax=ax12)
                ax12.legend(sectors['sector'], loc = 'upper right')
                ax12.set(xlabel='Range [m]', xlim = [lrange, urange], 
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         title=f'Normalized RC Signal (@ {int(norm_range)}'+r' $\pm$ '+f'{int(norm_hwin)} m)')
                ax12.minorticks_on()
                
                # Fig row with signal deviations
                #subplot 21
                dev_all.loc[:,ch,:].plot.line(hue='sector', ax = ax21)
                ax21.legend(sectors['sector'], loc = 'upper right')
                ax21.set(xlabel='Range [m]', xlim = [lrange, urange],
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         ylim = [min_dev_sig,max_dev_sig], yticks = sg_dev_ticks,
                         title='Signal Deviations')
                ax21.plot([lrange,urange],[dev_lim, dev_lim],'--',color='black')
                ax21.plot([lrange,urange],[-dev_lim,-dev_lim],'--',color='black')
                ax21.minorticks_on()
                
                #subplot 22
                ndev_all.loc[:,ch,:].plot.line(hue='sector', ax = ax22)
                ax22.legend(sectors['sector'], loc = 'upper right')
                ax22.set(xlabel='Range [m]', xlim = [lrange, urange],
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         ylim = [min_dev_nrm,max_dev_nrm], yticks = nrm_dev_ticks,
                         title='Normalized Signal Deviations')
                ax22.plot([lrange,urange],[dev_lim,dev_lim],'--',color='black')
                ax22.plot([lrange,urange],[-dev_lim,-dev_lim],'--',color='black')
                ax22.minorticks_on()                
                
                # Fig row with sector signal deviations
                #subplot 31
                dev_N2N.loc[ch,:].plot(color = 'magenta', ax = ax31)
                dev_EW.loc[ch,:].plot(color = 'cyan', ax = ax31)
                ax31.set(xlabel='Range [m]', xlim = [lrange, urange],
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         ylim = [min_dev_sig,max_dev_sig], yticks = sg_dev_ticks,
                         title='Sector signal deviations')
                ax31.legend(['N2 - N', 'E - W'], loc = 'upper right')
                ax31.plot([lrange,urange],[dev_lim,dev_lim],'--',color='black')
                ax31.plot([lrange,urange],[-dev_lim,-dev_lim],'--',color='black')      
                ax31.minorticks_on()
                
                #subplot 32
                ndev_N2N.loc[ch,:].plot(color = 'magenta', ax=ax32)
                ndev_EW.loc[ch,:].plot(color = 'cyan', ax=ax32)
                ax32.set(xlabel='Range [m]', xlim = [lrange, urange], 
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         ylim=[min_dev_nrm,max_dev_nrm], yticks = nrm_dev_ticks,
                         title='Sector normalized signal deviations')
                ax32.legend(['N - N2', 'E - W'], loc = 'upper right')
                ax32.plot([lrange,urange],[dev_lim,dev_lim],'--',color='black')
                ax32.plot([lrange,urange],[-dev_lim,-dev_lim],'--',color='black')            
                ax32.minorticks_on()
                
                fig.tight_layout()
                fig.subplots_adjust(top=0.86, hspace = 0.55, wspace = 0.15)
                                        #0.90
                # figure title
                title_name = f'Telecover test on channel {ch}\n'+\
                                f'Time frame: {stime_t}\n'+\
                                    f'Measurement angle: {sza}'+'$^{o}$ off-zenith'
                fig.suptitle(title_name, fontsize = 26)
                
                # dir for saving
                dir_plt = os.path.join(cfg.directories.plots, telc_l)
                
                if not os.path.exists(dir_plt):
                    os.makedirs(dir_plt)
                
                fpath = os.path.join(dir_plt, f'{dt_str}_{ch}.png')
                fig.savefig(fpath, dpi = dpi_val)
                plt.close('all')
    return()

def plot_telec_lite(telc_l, sig_telec_bc, sig_norm_rc, dev_all, ndev_all, 
                     dev_N2N, ndev_N2N, dev_EW, ndev_EW, 
                     sectors, info, datetime, cfg, options):

    # set the mpl rc params for plotting
    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':23, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz,
                        'xtick.major.size':xtck_mj_sz, 'ytick.major.size':ytck_mj_sz,
                        'xtick.minor.size':xtck_mn_sz, 'ytick.minor.size': ytck_mn_sz}):         
    
        sza = np.round(float(cfg.angles.zenith), decimals=1)
    
        
        # Extract plotting options
        dpi_val = 300 #int(options.dpi[0])
        
        dev_lim = float(options.crit_dev[0])#0.05
        
        # Normalization options
        norm_range = float(options.norm_range[0]) # range in meters
        norm_hwin   = float(options.norm_hwin[0]) # range in meters
        # limits of range axis (meters) and ticks step
        lrange = float(options.min_range[0])
        urange = float(options.max_range[0]) 
        step_range = float(options.step_range[0])
        # limits of raw signal axis and ticks step    
        min_dev_sig = float(options.min_dev_sig[0])
        max_dev_sig = float(options.max_dev_sig[0])
        step_dev_sig = float(options.step_dev_sig[0])
        # limits of normalized signal axis and ticks step
        min_dev_nrm = float(options.min_dev_nrm[0])
        max_dev_nrm = float(options.max_dev_nrm[0])
        step_dev_nrm = float(options.step_dev_nrm[0])
    
        rg_ticks = np.arange(lrange, urange+step_range, step_range)
        sg_dev_ticks = np.arange(min_dev_sig, max_dev_sig+step_dev_sig, step_dev_sig)
        nrm_dev_ticks = np.arange(min_dev_nrm, max_dev_nrm+step_dev_nrm, step_dev_nrm)
    
        
        stime_t = datetime.dt.strftime('%Y-%m-%d, %H:%M:%S UTC').values #%Y-%m-%dT%H:%M:%SZ
        dt_str = datetime.dt.strftime('%Y%m%dT%H%M').values
        
        print('-----------------------------------------')
        print('Generating telecover plots...')
        print('-----------------------------------------')
    
        for ch in sig_telec_bc.channel.values:
            
            # Make the plots if the signal values are not NaN in one of the sectors (i.e. empty for the rest sectors too)
            if info.loc[ch,'ch_mode']== 0. and np.isnan(sig_telec_bc.loc['north',ch,:].values).all() == False:
            
                y_max_bc  = np.round(sig_telec_bc.loc[:,ch,:].loc[:,lrange:urange].max()*1.5)
                y_max_nrc = np.round(sig_norm_rc.loc[:,ch,:].loc[:,lrange:urange].max()*1.5)
                # y_max_nrc = 1.5
                
                fig,(ax1,ax2) =\
                    plt.subplots(nrows=2,ncols=1,
                                 figsize = (14,10)) #20,16
                # Fig row with signals
    
                #subplot 1
                sig_norm_rc.loc[:,ch,:].plot.line(hue='sector', ylim = [0, y_max_nrc], ax=ax1)
                ax1.legend(sectors['sector'], loc = 'upper right')
                ax1.set(xlabel='Range [m]', xlim = [lrange, urange], 
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         ylabel='Signal [A.U.]',
                         title=f'Normalized RC Signal (@ {int(norm_range)}'+r' $\pm$ '+f'{int(norm_hwin)} m)')
                ax1.minorticks_on()
                
                # Fig row with signal deviations
                
                #subplot 2
                ndev_all.loc[:,ch,:].plot.line(hue='sector', ax = ax2)
                ax2.legend(sectors['sector'], loc = 'upper right')
                ax2.set(xlabel='Range [m]', xlim = [lrange, urange],
                         xticks = rg_ticks, #xticklabels = (rg_ticks/1e3).astype(str),
                         ylim = [min_dev_nrm,max_dev_nrm], yticks = nrm_dev_ticks,
                         ylabel='Signal Dev.',title='Normalized Signal Deviations')
                ax2.plot([lrange,urange],[dev_lim,dev_lim],'--',color='black')
                ax2.plot([lrange,urange],[-dev_lim,-dev_lim],'--',color='black')
                ax2.minorticks_on()                
                
                fig.tight_layout()
                fig.subplots_adjust(top=0.80, hspace = 0.55)
                #                         #0.86          0.55           0.15
                # figure title
                title_name = f'Telecover test on channel {ch}\n'+\
                                f'Time frame: {stime_t}\n'+\
                                    f'Measurement angle: {sza}'+'$^{o}$ off-zenith'
                fig.suptitle(title_name, fontsize = 24)
    
                # dir for saving
                dir_plt = os.path.join(cfg.directories.plots, telc_l)
                
                if not os.path.exists(dir_plt):
                    os.makedirs(dir_plt)
                
                fpath = os.path.join(dir_plt, f'{dt_str}_{ch}_lite.png')
                fig.savefig(fpath, dpi = dpi_val)
                plt.close('all')
    
    return()