"""
@authors: N. Siomos & P. Paschou
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime as dt

ln_wdth = 2.5 # 2
ax_lblsz = 18
ax_tlsz = 15.5#15
xtck_lblsz = 18
ytck_lblsz = 18
ax_lblpd = 15.
lgd_fsz = 15#18

# linecolor plotting
eve_c = 'tab:blue'
lid_c = 'tab:red'
lid_nf_c = 'tab:purple'
# err_c = 'grey'

# linestyle
eve_ls = 'solid'
lid_ls = 'dashed'
lid_nf_ls = 'dashed'

#fig size
x_inch = 7 #7 # 4.5
y_inch = 10.1 #7 # 6.5

# Plot gridlines in plots?
gd = True
gd_c = 'grey'
gd_l = 'dashed'
# Plot minorticks in plots?
minor_ticks = True

mLDR = 0.00586 # 0.00592 # paper
mLDR_err = 0.00004
mCDR = 0.011904 # 0.012 # paper
mCDR_err = 0.00009

def lineplots(title, fpath, dpiv, axeslims, size, yscale, 
              lstyle, legend, colors, x_label, y_label,
              y, x, xerr, nsigma, llim, ulim):

    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

    
        # plt.figure(figsize = size)
        fig, ax = plt.subplots(figsize=size)    
        # ax.set_title(title)
    
        for i in range(len(legend)):
            # mask data for the plotted lims
            mask_lim = (y[i] > llim[i]) & (y[i] < ulim[i])
    
            ax.plot(x[i][mask_lim], y[i][mask_lim]*yscale, linestyle=lstyle[i],
                     color = colors[i])
                
        if len(legend)>1:
            ax.legend(legend,loc=1)
    
        profllim = axeslims[0]
        if profllim != 0:
            ax.axvline(ls='-', lw='1.5', c='k', alpha = 0.6)
        
        for i in range(len(legend)):
            mask_lim = (y[i] > llim[i]) & (y[i] < ulim[i])
    
            ax.fill_betweenx(y[i][mask_lim]*yscale,
                              x[i][mask_lim] - nsigma*xerr[i][mask_lim], 
                              x[i][mask_lim] + nsigma*xerr[i][mask_lim],
                              alpha = 0.2, color = colors[i]) #err_c) #
    
        # if 'VLDR' in legend:
        #     ax.axvline(x=mLDR, ymin=0, ymax=1, ls='--', lw='2.', c='k', alpha = 0.2)
        #     if 'VCDR' in legend:
        #         ax.axvline(x=mCDR, ymin=0, ymax=1, ls='--', lw='2.', c='k', alpha = 0.2)
    
        ax.set_title(title, fontsize=20)
        ax.axis(axeslims)
        if gd:
            ax.grid(which='major',axis='both', alpha=0.2, linestyle=gd_l, color=gd_c)
        if minor_ticks:
            ax.minorticks_on()
        ax.set(ylabel = y_label, xlabel = x_label)
        fig.tight_layout()
        fig.savefig(fpath, dpi=dpiv)
        plt.close()
    
    return()


def lineplots_all_in_one(gen_title, fpath, y, profs, profs_err, xlabels, legends, 
                         lstyles, lcolors, axslims, llims, ulims, yscale, ylabel, nsigma,
                         dpiv, fig_size=[21, 8.]):

    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

    
        sbplt_num = len(profs)
        fig, axs = plt.subplots(nrows=1, ncols=sbplt_num, sharey=False, figsize=fig_size)
        
        for sb in range(sbplt_num):
        #     axs = 0
            lgd = legends[sb]
            x = profs[sb]
            xerr = profs_err[sb]
            lstyle = lstyles[sb]
            axeslims = axslims[sb]
            lcolor = lcolors[sb]
            llim = llims[sb]
            ulim = ulims[sb]
            
            x_label = xlabels[sb]
            
            for i in range(len(lgd)):
                # mask data for the plotted lims
                mask_lim = (y > llim) & (y < ulim)
        
                axs[sb].plot(x[i][mask_lim], y[mask_lim]*yscale, linestyle=lstyle[i],
                            color = lcolor[i])
    
            if len(lgd)>1:
                axs[sb].legend(lgd,loc=1, fontsize=17)
    
            profllim = axeslims[0]
            if profllim != 0:
                axs[sb].axvline(ls='-', lw=1.5, c='k', alpha = 0.5)
                
            for i in range(len(lgd)):
                mask_lim = (y > llim) & (y < ulim)
        
                axs[sb].fill_betweenx(y[mask_lim]*yscale,
                                      x[i][mask_lim] - nsigma*xerr[i][mask_lim], 
                                      x[i][mask_lim] + nsigma*xerr[i][mask_lim],
                                      alpha = 0.15, color = lcolor[i]) # alpha = 0.25, color = err_c) # 
    
            # if 'VLDR' in lgd:
            #     axs[sb].axvline(x=mLDR, ymin=0, ymax=1, ls='--', lw='2.', c='k', alpha = 0.2)
            #     if 'VCDR' in lgd:
            #         axs[sb].axvline(x=mCDR, ymin=0, ymax=1, ls='--', lw='2.', c='k', alpha = 0.3)
    
            if gd:
                axs[sb].grid(which='major',axis='both', alpha=0.2, linestyle=gd_l, color=gd_c)
    
            if sb == 0:
                axs[sb].set_ylabel(ylabel)
                
            axs[sb].set_xlabel(x_label)
            axs[sb].axis(axeslims)
            
            if minor_ticks:
                axs[sb].minorticks_on()
        
        fig.suptitle(gen_title, fontsize = 20)
        fig.tight_layout()
        fig.savefig(fpath,dpi=dpiv)
        plt.close()
    
    return()


def ovp(title, fpath, dpi, xlims, ylims, size, yscale, 
        style, ltype, legend, colors, xlabel, ylabel,
        y, x, xerr, llim, ulim):

    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, #'axes.labelpad':12.,
                        'axes.labelsize':15, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':15., 'ytick.labelsize':15.}):         

        fig, ax = plt.subplots(figsize = size)
        
        ax.set_title(title, y=1.01)
    
        for i in range(len(legend)):
            mask_lim = (y[i] > llim[i]/yscale[i]) & (y[i] < ulim[i]/yscale[i])
    
            if style[i] == 'line' and len(x[i][mask_lim]) > 0:
                ax.plot(x[i][mask_lim], y[i][mask_lim]*yscale[i], ltype[i],
                         color = colors[i])
                
            if style[i] == 'step' and len(x[i][mask_lim]) > 0:
                ax.step(x[i][mask_lim], y[i][mask_lim]*yscale[i], ltype[i], 
                         where = 'post', color = colors[i])
        
        ax.legend(legend)
        
        for i in range(len(legend)):
            mask_lim = (y[i] > llim[i]/yscale[i]) & (y[i] < ulim[i]/yscale[i])
    
            if style[i] == 'line' and len(xerr[i][mask_lim]) > 0:
                ax.fill_betweenx(y[i][mask_lim]*yscale[i],
                                  x[i][mask_lim] - xerr[i][mask_lim], 
                                  x[i][mask_lim] + xerr[i][mask_lim],
                                  alpha = 0.2, color = colors[i])
                
            if style[i] == 'step' and len(xerr[i][mask_lim]) > 0:    
                ax.step(x[i][mask_lim] - xerr[i][mask_lim], 
                         y[i][mask_lim]*yscale[i], '--', 
                         where = 'post', color = colors[i])
                ax.step(x[i][mask_lim] + xerr[i][mask_lim], 
                         y[i][mask_lim]*yscale[i], '--', 
                         where = 'post', color = colors[i])
        
        profllim = xlims[0]
        if profllim != 0:
            ax.axvline(ls='-', lw='1.5', c='k', alpha = 0.6)

        ax.set(xlim = xlims, ylim=ylims, xlabel = xlabel, ylabel = ylabel)
    
        if minor_ticks:
            ax.minorticks_on()
        if gd:
            ax.grid(which='major',axis='both', alpha=0.2, linestyle=gd_l, color=gd_c)
    
        fig.tight_layout()
        fig.savefig(fpath, dpi = dpi)
        plt.close()
    
    return()

def lineplots_comp(eve_pack, eve_alt, eve_msk, 
                    lid_name, lid_pack, lid_alt, lid_msk, lid_nf_msk, 
                    nsigma, axes_lims, lid_loc,
                    evet_start, evet_end, lid_start, lid_end,
                    fig_type,x_label, ver, save_dir, dpival=150):

    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         


        # eve timeframe - Time info 
        t = 0
        eve_start, eve_end, out_t, str_d = handle_time(evet_start, evet_end, t)
    
        # eve    
        eve_prof = eve_pack[0][t,eve_msk]
        eve_err = eve_pack[1][t,eve_msk]
        eve_h = eve_alt[eve_msk]
       
        # lidar for comparison 
        lid_prof = lid_pack[0][lid_msk]
        lid_err = lid_pack[1][lid_msk]
        
        # near field channels
        inf = False
        if len(lid_pack)>2 and len(lid_nf_msk)>0:
            inf=True
            lid_nf_prof = lid_pack[2][lid_nf_msk]
            lid_nf_err = lid_pack[3][lid_nf_msk]
            lid_nf_alt = lid_pack[4]
            lid_nf_h = lid_nf_alt[lid_nf_msk]
        
        # figure
        fig, ax = plt.subplots(figsize=(x_inch,y_inch))    
    
        eve_line = ax.plot(eve_prof, eve_h, color=eve_c, linestyle=eve_ls, label='$eVe$')
        lid_line = ax.plot(lid_prof, lid_alt[lid_msk], color=lid_c, linestyle=lid_ls, label=lid_name)
        
        ax.fill_betweenx(eve_h,
                          eve_prof - nsigma*eve_err,
                          eve_prof + nsigma*eve_err,
                          alpha=0.2, color=eve_c)
        ax.fill_betweenx(lid_alt[lid_msk],
                          (lid_prof - nsigma*lid_err),
                          (lid_prof + nsigma*lid_err),
                          alpha=0.2, color=lid_c)
        nf=''
        if inf:
            lid_line_nf = ax.plot(lid_nf_prof, lid_nf_h, color=lid_nf_c, linestyle=lid_nf_ls, label=lid_name+' $nf$')    
            ax.fill_betweenx(lid_nf_h,
                              (lid_nf_prof - nsigma*lid_nf_err),
                              (lid_nf_prof + nsigma*lid_nf_err),
                              alpha=0.2, color=lid_nf_c)
            nf = '_withnf'
        
        ax.legend(loc=1)       
        
        fig_title = f'Location: {lid_loc}\n'+\
            f'Date: {str_d}\n'+\
            f'eVe: {eve_start}-{eve_end} UTC\n'+\
                f'{lid_name}: {lid_start}-{lid_end} UTC'
        ax.set_title(fig_title, fontsize=20)
        
        ax.axis(axes_lims)
        if gd:
            ax.grid(which='major',axis='both', alpha=0.2, linestyle=gd_l, color=gd_c)    
        
        if minor_ticks:
            ax.minorticks_on()
    
        ax.set(xlabel = x_label, ylabel = 'Altitude a.s.l. [Km]')
        fig.tight_layout()                
        fig_name = f'{fig_type}_{out_t}{nf}{ver}.png'
        fig.savefig(os.path.join(save_dir, fig_name), dpi=dpival)
        plt.close()    
    
    return()

def lineplots_comp_all_in_one(eve_mpack, eve_alt, eve_msk, 
                              lid_name, lid_mpack, lid_alt, lid_msk, lid_nf_msk, 
                              nsigma, axes_lims_pck, lid_loc,
                              evet_start, evet_end, lid_start, lid_end,
                              x_label_pck, ver, save_dir, dpival=150, fig_size=[23., 10.]):

    with mpl.rc_context({'lines.linewidth':ln_wdth,
                        'axes.titlesize':ax_tlsz, 'axes.labelpad':ax_lblpd,
                        'axes.labelsize':ax_lblsz, 'legend.fontsize':lgd_fsz,
                        'xtick.labelsize':xtck_lblsz, 'ytick.labelsize':ytck_lblsz}):         

        sbplt_num = len(eve_mpack)
        fig, axs = plt.subplots(nrows=1, ncols=sbplt_num, sharey=False, figsize=fig_size)
    
        # eve timeframe - Time info 
        t = 0
        eve_start, eve_end, out_t, str_d = handle_time(evet_start, evet_end, t)
    
        for sb in range(sbplt_num):
            eve_pack = eve_mpack[sb]
            lid_pack = lid_mpack[sb]
            axes_lims = axes_lims_pck[sb]
            x_label = x_label_pck[sb]
            
            # eve    
            eve_prof = eve_pack[0][t,eve_msk]
            eve_err = eve_pack[1][t,eve_msk]
            eve_h = eve_alt[eve_msk]
               
            # lidar for comparison 
            lid_prof = lid_pack[0][lid_msk]
            lid_err = lid_pack[1][lid_msk]
            
            # near field channels
            inf = False
            if len(lid_pack)>2 and len(lid_nf_msk)>0:
                inf=True
                lid_nf_prof = lid_pack[2][lid_nf_msk]
                lid_nf_err = lid_pack[3][lid_nf_msk]
                lid_nf_alt = lid_pack[4]
                lid_nf_h = lid_nf_alt[lid_nf_msk]
            
            eve_line = axs[sb].plot(eve_prof, eve_h, color=eve_c, linestyle=eve_ls, label='$eVe$')
            lid_line = axs[sb].plot(lid_prof, lid_alt[lid_msk], color=lid_c, linestyle=lid_ls, label=lid_name)
            
            axs[sb].fill_betweenx(eve_h,
                              eve_prof - nsigma*eve_err,
                              eve_prof + nsigma*eve_err,
                              alpha=0.2, color=eve_c)
            axs[sb].fill_betweenx(lid_alt[lid_msk],
                              (lid_prof - nsigma*lid_err),
                              (lid_prof + nsigma*lid_err),
                              alpha=0.2, color=lid_c)
            if inf:
                lid_line_nf = axs[sb].plot(lid_nf_prof, lid_nf_h, color=lid_nf_c, linestyle=lid_nf_ls, label=lid_name+' $nf$')    
                axs[sb].fill_betweenx(lid_nf_h,
                                  (lid_nf_prof - nsigma*lid_nf_err),
                                  (lid_nf_prof + nsigma*lid_nf_err),
                                  alpha=0.2, color=lid_nf_c)
            
            axs[sb].legend(loc=1)       
            
            axs[sb].axis(axes_lims)
            if gd:
                axs[sb].grid(which='major',axis='both', alpha=0.2, linestyle=gd_l, color=gd_c)    
            
            if minor_ticks:
                axs[sb].minorticks_on()
                
            axs[sb].set(xlabel = x_label, ylabel = 'Altitude a.s.l. [Km]')
            axs[sb].set_title('profiles at 355 nm', fontsize=18)
            
        nf=''
        if len(lid_nf_msk)>0:
            nf = '_withnf'
            
        fig_title = f'Location: {lid_loc}\n'+\
            f'Date: {str_d}\n'+\
            f'eVe: {eve_start}-{eve_end} UTC\n'+\
                f'{lid_name}: {lid_start}-{lid_end} UTC'
        fig.suptitle(fig_title, fontsize=20)
            
        fig.tight_layout()                
        fig_name = f'all_profiles_lid_comp_{out_t}{nf}{ver}.png'
        fig.savefig(os.path.join(save_dir, fig_name), dpi=dpival)
        plt.close()    
    
    return()

def make_title(lid_s_dt, lid_e_dt, lid_m_dt, lat, lon, st_long):
    
    title = []
    
    if isinstance(lid_m_dt,list) == False:
        lid_s_t = dt.datetime.strftime(lid_s_dt,'%H:%M:%S')
        lid_e_t = dt.datetime.strftime(lid_e_dt,'%H:%M:%S')
        date = dt.datetime.strftime(lid_m_dt,'%d %b %Y')        
                # \n\n
        title = f'{st_long} ({lat}$^o$ N, {lon}$^o$ E, 2m)  \n' +\
                f'{date}\neVe: {lid_s_t} - {lid_e_t} UTC'
        
    return(title)

def make_fname(dir_out, sat_dt, label):
    
    full_date = dt.datetime.strftime(sat_dt,'%Y%m%dT%H%M%S')
    
    fpath = os.path.join(dir_out,f'{label}_{full_date}')
    
    return(fpath)

def handle_time(t_start, t_end, t_ind):
    
    t1 = t_start[t_ind]
    start_dt = dt.datetime(int(t1[:4]), int(t1[4:6]), int(t1[6:8]),
                           int(t1[9:11]), int(t1[11:13]), int(t1[13:15]))
    start_t = start_dt.strftime('%H:%M')
    
    t2 = t_end[t_ind]
    end_dt = dt.datetime(int(t2[:4]), int(t2[4:6]), int(t2[6:8]),
                         int(t2[9:11]), int(t2[11:13]), int(t2[13:15]))
    end_t = end_dt.strftime('%H:%M')
    
    out_t = start_dt.strftime('%d%m%YT%H%M')
    
    str_d = start_dt.strftime('%d/%m/%Y')
    
    return(start_t, end_t, out_t, str_d)    