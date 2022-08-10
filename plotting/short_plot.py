"""
@author: nick
"""
import numpy as np
from plotting import official_plots
# from plotting import official_plots_old


x_inch = 5. #5. #4
y_inch = 10.5 #9. #8.5 #fig_size = [5., 10.5]Aeolus algo

dpiv = 300

yaxis_lims = [0, 10]
datallim = 0.4 #0.45
dataulim = 10.

# colors
lin_cl = 'blue' #'tab:blue'
cir_cl = 'tab:blue'#'dodgerblue'
aeol_cl= 'tab:cyan' #'m'# 'tab:red' #
earl_cl= 'slateblue'
conv_cir_cl= 'tab:olive'#'tab:cyan'#


def plot_bsc(dir_out, lid_dt, title, 
             alt_c, bsc_c, bsc_c_err, 
             alt_l, bsc_l, bsc_l_err, 
             alt_p_c, bsc_p_c, bsc_p_c_err,
             alt_p_l, bsc_p_l, bsc_p_l_err):


    xaxis_lims = [-1, 15.]   

# BSC LIN-CIR Plot     
    if len(bsc_c) > 0 and len(bsc_p_c) > 0 \
        and len(bsc_l) > 0 and len(bsc_p_l) > 0:

        #  only total bsc
        path = official_plots.make_fname(dir_out, lid_dt, 
                                         label = 'bsc_eve') #bsc_comp
    
        legend = ['Total - Lin. Pol. Emission', 'Total - Cir. Pol. Emission']
        colors = [lin_cl, cir_cl]
        style =  ['line', 'line']
        ltype =  ['-', '--']
        yscale = [1e-3, 1e-3]
        llim = list(np.ones(shape = len(legend))*datallim)
        ulim = list(np.ones(shape = len(legend))*dataulim)
        
        xlabel = 'Particle Backscatter 355nm [$Mm^{-1} sr^{-1}$]'
        ylabel = 'Height [km]'        
        bscs_y = [alt_l, alt_c]
        bscs_x = [bsc_l, bsc_c]
        bscs_xerr = [bsc_l_err, bsc_c_err]
        
        official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
                           xlims = xaxis_lims, ylims = yaxis_lims, 
                           size = [x_inch, y_inch], style = style, ltype = ltype,
                           legend = legend, colors = colors,
                           xlabel = xlabel, ylabel = ylabel,
                           y = bscs_y, x = bscs_x, xerr = bscs_xerr,
                           llim = llim, ulim = ulim)

#  only total cir bsc and aeol like

        path = official_plots.make_fname(dir_out, lid_dt, 
                                         label = 'bsc_eve_Aeol') #bsc_comp
        
        legend = ['Total - Cir. Pol. Emission', 'AEOLUS Like - Cir. Pol. Emission']
        colors = [cir_cl, aeol_cl]
        style =  ['line', 'line']
        ltype =  ['-', '--']
        yscale = [1e-3, 1e-3]
        llim = list(np.ones(shape = len(legend))*datallim)
        ulim = list(np.ones(shape = len(legend))*dataulim)
        
        xlabel = 'Particle Backscatter 355nm [$Mm^{-1} sr^{-1}$]'
        ylabel = 'Height [km]'        
        bscs_y = [alt_c, alt_p_c]
        bscs_x = [bsc_c, bsc_p_c]
        bscs_xerr = [bsc_c_err, bsc_p_c_err]

        
        official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
                           xlims = xaxis_lims, ylims = yaxis_lims, 
                           size = [x_inch, y_inch], style = style, ltype = ltype,
                           legend = legend, colors = colors,
                           xlabel = xlabel, ylabel = ylabel,
                           y = bscs_y, x = bscs_x, xerr = bscs_xerr,
                           llim = llim, ulim = ulim)

    
    return()

def plot_ext(dir_out, lid_dt, title, 
             alt_c, ext_c, ext_c_err, 
             alt_l, ext_l, ext_l_err):

    xaxis_lims = [-50,400] #900
# EXT LIN-CIR Plot     
    if len(ext_c) > 0 and len(ext_l) > 0:

        path = official_plots.make_fname(dir_out, lid_dt, 
                                         label = 'ext_eVe')
    
        legend = ['Lin. Pol. Emission', 'Cir. Pol. Emission']
        colors = [lin_cl, cir_cl]
        style =  ['line', 'line']
        ltype =  ['-', '--']
        yscale = [1e-3, 1e-3]
        llim = list(np.ones(shape = len(legend))*datallim)
        ulim = list(np.ones(shape = len(legend))*dataulim)
        
        xlabel = 'Particle Extinction 355nm [$Mm^{-1}$]'
        ylabel = 'Height [km]'        
        exts_y = [alt_l, alt_c]
        exts_x = [ext_l, ext_c]
        exts_xerr = [ext_l_err, ext_c_err]
        
        official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
                           xlims = xaxis_lims, ylims = yaxis_lims, 
                           size = [x_inch, y_inch], style = style, ltype = ltype,
                           legend = legend, colors = colors,
                           xlabel = xlabel, ylabel = ylabel,
                           y = exts_y, x = exts_x, xerr = exts_xerr,
                           llim = llim, ulim = ulim)        

        # path = official_plots.make_fname(dir_out, lid_dt, 
        #                                  label = 'ext_eVe_cir')
        # legend = ['Cir. Pol. Emission']
        # colors = [ cir_cl]
        # style =  [ 'line']
        # ltype =  [ '-']
        # yscale = [ 1e-3]
        # llim = [datallim]
        # ulim = [dataulim]
        
        # exts_y = [alt_c]
        # exts_x = [ext_c]
        # exts_xerr = [ext_c_err]
        
        # official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
        #                    xlims = xaxis_lims, ylims = yaxis_lims, 
        #                    size = [x_inch, y_inch], style = style, ltype = ltype,
        #                    legend = legend, colors = colors,
        #                    xlabel = xlabel, ylabel = ylabel,
        #                    y = exts_y, x = exts_x, xerr = exts_xerr,
        #                    llim = llim, ulim = ulim)        
        
    return()


def plot_lr(dir_out, lid_dt, title, 
            alt_c, lr_c, lr_c_err, 
            alt_l, lr_l, lr_l_err,
            alt_p_c, lr_p_c, lr_p_err,
            alt_p_l, lr_p_l, lr_p_l_err):

    xaxis_lims = [0, 150]
    lr_ulim = 5.4
# LR Plot     
    if len(lr_c) > 0 and len(lr_p_c) > 0 \
        and len(lr_l) > 0 and len(lr_p_l) > 0:


        path = official_plots.make_fname(dir_out, lid_dt, 
                                         label = 'lr_eVe')
    
        legend = ['Lin. Pol. Emission', 'Cir. Pol. Emission']
        colors = [lin_cl, cir_cl]
        style =  ['line', 'line']
        ltype =  ['-', '--']
        yscale = [1e-3, 1e-3]
        llim = list(np.ones(shape = len(legend))*datallim)
        ulim = list(np.ones(shape = len(legend))*lr_ulim)

        
        xlabel = 'Lidar Ratio 355nm [$sr$]'
        ylabel = 'Height [km]'        
        lrs_y = [alt_l, alt_c]
        lrs_x = [lr_l, lr_c]
        lrs_xerr = [lr_l_err, lr_c_err]
        
        official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
                           xlims = xaxis_lims, ylims = yaxis_lims, 
                           size = [x_inch, y_inch], style = style, ltype = ltype,
                           legend = legend, colors = colors,
                           xlabel = xlabel, ylabel = ylabel,
                           y = lrs_y, x = lrs_x, xerr = lrs_xerr,
                           llim = llim, ulim = ulim)  

        # eVe total vs Aeolus like
        path = official_plots.make_fname(dir_out, lid_dt, 
                                         label = 'lr_eVe_Aeol')

        legend = ['Cir. Pol. Emission', 'AEOLUS Like - Cir. Pol. Emission']
        colors = [cir_cl, aeol_cl]
        style =  ['line', 'line']
        ltype =  ['-','--']
        yscale = [1e-3, 1e-3]
        llim = list(np.ones(shape = len(legend))*datallim)
        ulim = list(np.ones(shape = len(legend))*lr_ulim)#dataulim)

        lrs_y = [alt_c, alt_p_c]
        lrs_x = [lr_c, lr_p_c]
        lrs_xerr = [lr_c_err, lr_p_err]
        
        official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
                           xlims = xaxis_lims, ylims = yaxis_lims, 
                           size = [x_inch, y_inch], style = style, ltype = ltype,
                           legend = legend, colors = colors,
                           xlabel = xlabel, ylabel = ylabel,
                           y = lrs_y, x = lrs_x, xerr = lrs_xerr,
                           llim = llim, ulim = ulim)  
        
    return()

def plot_vdr(dir_out, lid_dt, title, 
             alt_c, vdr_c, vdr_c_err, 
             alt_l, vdr_l, vdr_l_err,
             alt_lc, vdr_lc, vdr_lc_err):
    
# VDR Plot     
    if len(vdr_l) > 0 and len(vdr_c) > 0 and len(vdr_lc) > 0:

        path = official_plots.make_fname(dir_out, lid_dt, 
                                         label = 'vdr_comp')
    
        legend = ['VLDR', 'VCDR'] 
                  # 'Converted VLDR to VCDR']
        # legend = ['Lin. Pol. Emission', 'Cir. Pol. Emission', 
        #           'Converted Lin. to Cir.']
        colors = [lin_cl, cir_cl, conv_cir_cl] 
        # colors = ['tab:blue', 'tab:blue', 'tab:cyan']
        style =  ['line', 'line', 'line']
        ltype =  ['-', '--', ':']
        yscale = [1e-3, 1e-3, 1e-3]
        llim = list(np.ones(shape = len(legend))*datallim)
        ulim = list(np.ones(shape = len(legend))*dataulim)
        
        xlabel = 'Volume Depolarization Ratio 355nm'
        ylabel = 'Height [km]'        
        bscs_y = [alt_l, alt_c, alt_lc]
        bscs_x = [vdr_l, vdr_c, vdr_lc]
        bscs_xerr = [vdr_l_err, vdr_c_err, vdr_lc_err]
        
        official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
                           xlims = [-0.02, 0.4], ylims = yaxis_lims, 
                           size = [x_inch, y_inch], style = style, ltype = ltype,
                           legend = legend, colors = colors,
                           xlabel = xlabel, ylabel = ylabel,
                           y = bscs_y, x = bscs_x, xerr = bscs_xerr,
                           llim = llim, ulim = ulim)
    return()
        
def plot_pdr(dir_out, lid_dt, title, 
             alt_c, pdr_c, pdr_c_err, 
             alt_l, pdr_l, pdr_l_err,
             alt_lc, pdr_lc, pdr_lc_err):
# PDR Plot      
    pdr_ulim = 7.
    
    if len(pdr_l) > 0 and len(pdr_c) > 0 and len(pdr_lc) > 0:

        path = official_plots.make_fname(dir_out, lid_dt, 
                                         label = 'pdr_comp')
    
        legend = ['PLDR', 'PCDR']
                  # 'Converted PLDR to PCDR']
        # legend = ['Lin. Pol. Emission', 'Cir. Pol. Emission', 
        #           'Converted Lin. to Cir.']
        colors = [lin_cl, cir_cl, conv_cir_cl] 
        style =  ['line', 'line', 'line']
        ltype =  ['-', '--', ':']
        yscale = [1e-3, 1e-3, 1e-3]
        llim = list(np.ones(shape = len(legend))*datallim)
        # ulim = list(np.ones(shape = len(legend))*dataulim)
        ulim = list(np.ones(shape = len(legend))*pdr_ulim)
        
        xlabel = 'Particle Depolarization Ratio 355nm'
        ylabel = 'Height [km]'        
        bscs_y = [alt_l, alt_c, alt_lc]
        bscs_x = [pdr_l, pdr_c, pdr_lc]
        bscs_xerr = [pdr_l_err, pdr_c_err, pdr_lc_err]
        
        official_plots.ovp(title, path, dpi = dpiv, yscale = yscale,
                           xlims = [-0.05, 1.2], ylims = yaxis_lims, 
                           size = [x_inch, y_inch], style = style, ltype = ltype,
                           legend = legend, colors = colors,
                           xlabel = xlabel, ylabel = ylabel,
                           y = bscs_y, x = bscs_x, xerr = bscs_xerr,
                           llim = llim, ulim = ulim)
    
    return()