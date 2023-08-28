#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:26:22 2022

@author: nick
"""

from ..plotting import color_lib, make_colormap
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import os
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from subprocess import run

def quicklook(dir_out, fname, title, dpi_val, color_reduction, use_log,
              x_vals, t_vals, y_vals, z_vals, 
              x_lbin, x_ubin, y_lbin, y_ubin, 
              y_llim, y_ulim, z_llim, z_ulim,
              x_label, t_label, y_label, 
              x_tick, t_tick, y_tick, nodes):

    z_vals = z_vals.T
    
    # Convert to NaN wherever files are missing
    if len(nodes) > 0:
        z_vals[:,nodes] = np.nan
        z_vals[:,nodes+1] = np.nan

    # Define the colorscale
    rgb = color_lib.volkers_rgb()
    my_cmap = make_colormap.custom_rgb(rgb, name = 'volkers')

    # Create the variables to be plotted X, Y, Z
    [X, Y] = np.meshgrid(t_vals[slice(x_lbin, x_ubin+1)], 
                         y_vals[slice(y_lbin, y_ubin+1)])
    Z = z_vals[slice(y_lbin, y_ubin+1), slice(x_lbin, x_ubin+1)]
    
    # Create the figure
    fig = plt.figure(figsize=(10. , 5.))
    ax = fig.add_axes([0.07,0.13,0.99,0.68])
    
    ax.set_title(title, pad = 15) 
    
    ax.set_xlabel(t_label)
    ax.set_ylabel(y_label)

    y_ticks = np.arange(y_tick * np.ceil(y_llim / y_tick), 
                        y_tick * (np.floor(y_ulim / y_tick) + 1.), 
                        y_tick)
    
    if np.abs(y_llim - y_ticks[0]) < y_tick * 0.25:
        y_ticks[0] = y_llim
    else:
        y_ticks = np.hstack((y_llim, y_ticks))

    if np.abs(y_ulim - y_ticks[-1]) < y_tick * 0.25:
        y_ticks[-1] = y_ulim
    else:
        y_ticks = np.hstack((y_ticks, y_ulim))
        
    y_ticks = np.round(y_ticks, decimals = 2)
    
    ax.set_yticks(y_ticks, labels = y_ticks)
    ax.set_ylim([y_llim, y_ulim])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    if t_tick >= 4.:
        major_locator = mdates.MinuteLocator(byminute = np.arange(0,60, t_tick).astype(int))
        major_locator.MAXTICKS = 10000
        
        ax.xaxis.set_major_locator(major_locator)
        
        minor_locator = mdates.MinuteLocator(byminute = np.arange(0,60, t_tick/4.).astype(int))
        minor_locator.MAXTICKS = 10000
        
        ax.xaxis.set_minor_locator(minor_locator)
    else:
        major_locator = mdates.MinuteLocator(byminute = np.arange(0, 60, t_tick).astype(int))
        major_locator.MAXTICKS = 10000
        
        ax.xaxis.set_major_locator(major_locator)
        
        minor_locator = mdates.SecondLocator(bysecond = np.arange(0,60, 60. * t_tick/4.).astype(int))
        minor_locator.MAXTICKS = 10000
        
        ax.xaxis.set_minor_locator(minor_locator)


    plt.xticks(rotation = 35)

    if len(nodes) == 0:
        ax1 = ax.twiny()
    
        x_ticks = np.arange(x_lbin, x_ubin + x_tick, x_tick).astype(int)
        ax1.set_xticks(x_ticks, labels = x_ticks)
    
    if use_log == False:
        qck = ax.pcolormesh(X, Y, Z, vmin = z_llim, vmax = z_ulim, cmap = my_cmap)
    else:
        qck = ax.pcolormesh(X, Y, Z, vmin = z_llim, vmax = z_ulim, cmap = my_cmap, norm = 'log')
    
    fig.colorbar(qck, label = 'Range-corrected Signal [A.U.]', 
                 extend = 'both', pad = 0.02)
    
    fpath = os.path.join(dir_out, 'plots', fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()

    if color_reduction == True:
        run(["convert", fpath, "-posterize", "8", fpath])
    
    return(fpath)

def rayleigh(dir_out, fname, title, dpi_val, color_reduction, use_lin, norm_region,
             x_vals, y1_vals, y2_vals, y1_errs, coef, rsem, rslope, pval,
             x_lbin, x_ubin, x_llim, x_ulim, y_llim, y_ulim, 
             x_label, y_label, x_tick):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]
    Y1 = coef * y1_vals[slice(x_lbin, x_ubin)]
    Y2 = y2_vals[slice(x_lbin, x_ubin)]
    Y1E = coef * y1_errs[slice(x_lbin, x_ubin)]
    
    # Create the figure
    fig = plt.figure(figsize=(12. , 4.))
    fig.suptitle(title)

    ax = fig.add_axes([0.07,0.13,0.50,0.7])
        
    ax.plot(X, Y1, color = 'tab:blue', label = 'measured')
    ax.plot(X, Y2, color = 'tab:red', label = 'molecular')

    if np.isnan(Y1E).all() == False:
        ax.fill_between(X, Y1 - Y1E, Y1 + Y1E, color = 'tab:blue', alpha = 0.3)
    
    x_ticks = np.arange(x_tick * np.ceil(x_llim / x_tick), 
                        x_tick * (np.floor(x_ulim / x_tick) + 1.), 
                        x_tick)
    
    if x_tick >= x_ulim - x_llim:
        raise Exception(f"The x_tick ({x_tick}) must be smaller than the width of the normalization_region ({norm_region}) for the Rayleigh fit test. Please revise the settings_file.ini ")
    
    if x_llim > norm_region[1] or x_ulim < norm_region[0]:
        raise Exception(f"The normalization_region ({norm_region}) for the Rayleigh fit is out of the provided x_lims ([{x_llim}, {x_ulim}]). Please revise the settings_file.ini ")
    
    
    if np.abs(x_llim - x_ticks[0]) < x_tick * 0.25:
        x_ticks[0] = x_llim
    else:
        x_ticks = np.hstack((x_llim, x_ticks))

    if np.abs(x_ulim - x_ticks[-1]) < x_tick * 0.25:
        x_ticks[-1] = x_ulim
    else:
        x_ticks = np.hstack((x_ticks, x_ulim))

    x_ticks = np.round(x_ticks, decimals = 2)

    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label)
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.))

    ax.set_ylim([y_llim, y_ulim])
    ax.set_ylabel(y_label)
    if use_lin == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    ax.legend(loc = 'lower left')

    ax.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')

    n_llim = np.round(norm_region[0], decimals = 2)
    n_ulim = np.round(norm_region[1], decimals = 2)
    
    if rsem > 0.01:
        c_rsem = 'tab:red'
    else:
        c_rsem = 'tab:green'

    if pval < 0.05:
        c_pval = 'tab:red'
    else:
        c_pval = 'tab:green'
        
    if use_lin == False:

        ax.text(0.55 * x_ulim, 0.60 * y_ulim, 
                f'norm. region: {n_llim} - {n_ulim} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
            
        ax.text(0.55 * x_ulim, 0.30 * y_ulim, 
                f'rsem: {np.round(rsem, decimals = 4)}',
                bbox = dict(facecolor = c_rsem, alpha = 0.1, zorder = 3))
        
        ax.text(0.55 * x_ulim, 0.15 * y_ulim, 
                f'rslope: {np.round(rslope, decimals = 4)}',
                bbox = dict(facecolor = c_pval, alpha = 0.1, zorder = 3))
    else:
        ax.text(0.55 * x_ulim, 0.9 * y_ulim, 
                f'norm. region: {n_llim} - {n_ulim} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))

        ax.text(0.55 * x_ulim, 0.82 * y_ulim, 
                f'rsem: {np.round(rsem, decimals = 4)}',
                bbox = dict(facecolor = c_rsem, alpha = 0.1, zorder = 3))

        ax.text(0.55 * x_ulim, 0.74 * y_ulim, 
                f'rslope: {np.round(rslope, decimals = 4)}',
                bbox = dict(facecolor = c_pval, alpha = 0.1, zorder = 3))

    ax2 = fig.add_axes([0.65,0.13,0.30,0.7])
    
    if np.isnan(Y1E).all() == False:
        ax2.fill_between(X, (Y1 - Y1E - Y2) / Y2, 
                         (Y1 + Y1E - Y2) / Y2, color = 'tab:blue', 
                         alpha = 0.3, label = 'sem')
        ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue',label = 'mean')
        ax2.legend(loc = 'lower left')

    else:
        ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue')
        
    
    ax2.axhline(c = 'k')

    x_tick_2 = 2. * x_tick 
    x_ticks_2 = np.arange(x_tick_2 * np.floor(x_llim / x_tick_2), 
                          x_tick_2 * (np.ceil(x_ulim / x_tick_2) + 1.), 
                          x_tick_2)
    
    ax2.set_xticks(x_ticks_2, labels = x_ticks_2)
    ax2.set_xlim([x_llim, x_ulim])
    ax2.set_xlabel(x_label)
    ax2.xaxis.set_minor_locator(MultipleLocator(x_tick_2 / 2.))

    y_ticks = np.round(np.arange(-0.40, 0.40 + 0.10, 0.10), decimals = 2)
    ax2.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax2.set_ylim([y_ticks[0], y_ticks[-1]])
    ax2.set_ylabel('Relative Attenuated Bsc. Diff. ')
    
    ax2.grid(which = 'both')
    
    ax2.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')
    
    fpath = os.path.join(dir_out, 'plots', fname)
            
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    if color_reduction == True:
        run(["convert", fpath, "-posterize", "8", fpath])
    
    return(fpath)

def rayleigh_mask(dir_out, fname, title, dpi_val, color_reduction,
                  mmol, mder, cder, mshp, mneg, rsem):

    rgb = color_lib.volkers_rgb()
    vlk_cmap = make_colormap.custom_rgb(rgb, name = 'volkers')

    [X, Y] = np.meshgrid(mmol.lower_limit.values, mmol.window)
    fig = plt.figure(figsize=(12. , 8.))
    fig.suptitle(title)
    
    x_llim = 0.
    x_ulim = 16.
    
    y_llim = 1.
    y_ulim = 4.
    
    fig_x = 0.44
    fig_y = 0.23
    
    fig_edg1_x = 0.06
    fig_edg2_x = 0.54
    
    fig_edg1_y = 0.07
    fig_edg2_y = 0.36
    fig_edg3_y = 0.65
    
    ax = fig.add_axes([fig_edg1_x, fig_edg3_y, fig_x, fig_y])
    ax.pcolormesh(X, Y, mder.values)
    ax.set_title('Derivative mask', pad = 5)
    ax.set_ylabel('Window [km]')
    ax.set_ylim([y_llim, y_ulim])
    ax.set_xlim([x_llim, x_ulim])

    ax2 = fig.add_axes([fig_edg2_x, fig_edg3_y, fig_x, fig_y])
    ax2.pcolormesh(X, Y, rsem.values <= 0.02)
    ax2.set_title('Relative SEM mask', pad = 5)
    ax2.set_ylim([y_llim, y_ulim])
    ax2.set_xlim([x_llim, x_ulim])
    
    ax3 = fig.add_axes([fig_edg1_x, fig_edg2_y, fig_x, fig_y])
    ax3.pcolormesh(X, Y, cder.values)
    ax3.set_title('Curvature mask', pad = 5)
    ax3.set_ylabel('Window [km]')
    ax3.set_ylim([y_llim, y_ulim])
    ax3.set_xlim([x_llim, x_ulim])       

    ax4 = fig.add_axes([fig_edg2_x, fig_edg2_y, fig_x, fig_y])
    ax4.pcolormesh(X, Y, mshp.values)
    ax4.set_title('Shapiro-Wilk mask', pad = 5)
    ax4.set_ylim([y_llim, y_ulim])
    ax4.set_xlim([x_llim, x_ulim])

    ax5 = fig.add_axes([fig_edg1_x, fig_edg1_y, fig_x, fig_y])
    ax5.pcolormesh(X, Y, mneg.values)
    ax5.set_title('Cross-check mask', pad = 5)
    ax5.set_xlabel('Lower Limit [km]')
    ax5.set_ylabel('Window [km]')
    ax5.set_ylim([y_llim, y_ulim])
    ax5.set_xlim([x_llim, x_ulim])

    ax6 = fig.add_axes([fig_edg2_x, fig_edg1_y, fig_x, fig_y])
    plot6 = ax6.pcolormesh(X, Y, 100. * rsem.where(mmol).values, vmin = 0., vmax = 2., cmap = vlk_cmap)
    ax6.set_title('Masked Relative SEM (%)', pad = 5)
    ax6.grid(which = 'both')
    ax6.axes.xaxis.set_minor_locator(MultipleLocator(1))
    ax6.set_xlabel('Lower Limit [km]')
    ax6.set_ylim([y_llim, y_ulim])
    ax6.set_xlim([x_llim, x_ulim])
    
    cax = fig.add_axes([fig_edg2_x + 0.01, fig_edg1_y +0.01, 0.01, fig_y -0.02])

    fig.colorbar(plot6, cax=cax, orientation='vertical')
    
    fpath = os.path.join(dir_out, 'plots', fname)

    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    if color_reduction == True:
        run(["convert", fpath, "-posterize", "8", fpath])
    
    return(fpath)


def telecover_sec(dir_out, fname, title, dpi_val, color_reduction, 
                  norm_region, use_nonrc, x_vals, 
                  y1_raw, y2_raw, y3_raw, y4_raw,
                  y1_vals, y2_vals, y3_vals, y4_vals,
                  y1_extr, y2_extr, y3_extr, y4_extr,
                  y1_extr_raw, y2_extr_raw, y3_extr_raw, y4_extr_raw,
                  y1_lvar, y2_lvar, y3_lvar, y4_lvar, 
                  y1_uvar, y2_uvar, y3_uvar, y4_uvar, 
                  coef_1, coef_2, coef_3, coef_4, 
                  coef_extra_1, coef_extra_2, coef_extra_3, coef_extra_4, 
                  extra_sec, ranges,
                  x_lbin, x_ubin, x_llim, x_ulim, 
                  y_llim, y_ulim, y_llim_nr, y_ulim_nr, 
                  x_label, x_tick, use_last, iters):

    # Create the variables to be plotted X, Y
    X = x_vals
    
    if use_nonrc == True:
        R = ranges
    
        Y1 = y1_raw / np.power(R, 2)
        Y2 = y2_raw / np.power(R, 2)
        Y3 = y3_raw / np.power(R, 2)
        Y4 = y4_raw / np.power(R, 2)
        
        y_label_1 = 'Non RC Signals - Raw Units'
    
    else:    
        Y1 = y1_raw
        Y2 = y2_raw
        Y3 = y3_raw
        Y4 = y4_raw
        
        y_label_1 = 'RC Signals [A.U.]'

    
    Y1_N = coef_1 * y1_vals
    Y2_N = coef_2 * y2_vals
    Y3_N = coef_3 * y3_vals
    Y4_N = coef_4 * y4_vals

    Y1_NL = coef_1 * y1_lvar
    Y2_NL = coef_2 * y2_lvar
    Y3_NL = coef_3 * y3_lvar
    Y4_NL = coef_4 * y4_lvar
    
    Y1_NU = coef_1 * y1_uvar
    Y2_NU = coef_2 * y2_uvar
    Y3_NU = coef_3 * y3_uvar
    Y4_NU = coef_4 * y4_uvar
    
    Y_NM = np.mean([Y1_N, Y2_N, Y3_N, Y4_N], axis = 0)
    Y_RMSE = np.sqrt((np.power((Y1_N - Y_NM) / Y_NM, 2) +\
                      np.power((Y2_N - Y_NM) / Y_NM, 2) +\
                      np.power((Y3_N - Y_NM) / Y_NM, 2) +\
                      np.power((Y4_N - Y_NM) / Y_NM, 2)) / 4)

    Y_E = np.nan * Y_NM
    Y_O = np.nan * Y_NM
    
    Y_E_N = np.nan * Y_NM
    
    if extra_sec['N']:
        Y_E = y1_extr_raw
        Y_E_N = coef_extra_1 * y1_extr
        Y_O = Y1_N
        extra_label = f'north{iters+1}'
        extra_label_short = f'N{iters+1}'
        extra_exists = True
    elif extra_sec['E']:
        Y_E = y2_extr_raw
        Y_E_N = coef_extra_2 * y2_extr
        Y_O = Y2_N
        extra_label = f'east{iters+1}'
        extra_label_short = f'E{iters+1}'
        extra_exists = True
    elif extra_sec['S']:
        Y_E = y3_extr_raw
        Y_E_N = coef_extra_3 * y3_extr
        Y_O = Y3_N
        extra_label = f'south{iters+1}'
        extra_label_short = f'S{iters+1}'
        extra_exists = True
    elif extra_sec['W']:
        Y_E = y4_extr_raw
        Y_E_N = coef_extra_4 * y4_extr
        Y_O = Y4_N
        extra_label = f'west{iters+1}'
        extra_label_short = f'W{iters+1}'
        extra_exists = True
    else:
        extra_exists = False

    Y_DIFF = (Y_E_N - Y_O) / Y_O
    
    x_ticks = np.arange(x_tick * np.ceil(x_llim / x_tick), 
                        x_tick * (np.floor(x_ulim / x_tick) + 1.), 
                        x_tick)
        
    if x_tick >= x_ulim - x_llim:
        raise Exception(f"The x_tick ({x_tick}) must be smaller than the width of the normalization_region ({norm_region}) for the telecover test. Please revise the settings_file.ini ")
    
    if x_llim > norm_region[1] or x_ulim < norm_region[0]:
        raise Exception(f"The normalization_region ({norm_region}) for the telecover fit is out of the provided x_lims ([{x_llim}, {x_ulim}]). Please revise the settings_file.ini ")
    
    if np.abs(x_llim - x_ticks[0]) < x_tick * 0.25:
        x_ticks[0] = x_llim
    else:
        x_ticks = np.hstack((x_llim, x_ticks))

    if np.abs(x_ulim - x_ticks[-1]) < x_tick * 0.25:
        x_ticks[-1] = x_ulim
    else:
        x_ticks = np.hstack((x_ticks, x_ulim))
    
    x_ticks = np.round(x_ticks, decimals = 2)

    plt.rc('font', size = 14) 

    # Create the figure
    fig = plt.figure(figsize=(20. , 4.))
    
    fig.suptitle(title)

    # Subplot: Raw Signals - Near range
    ax = fig.add_axes([0.04,0.13,0.19,0.7])
    ax2 = fig.add_axes([0.24,0.13,0.07,0.7])
    ax3 = fig.add_axes([0.36,0.13,0.19,0.7])
    ax4 = fig.add_axes([0.56,0.13,0.07,0.7])
    ax5 = fig.add_axes([0.69,0.13,0.29,0.7])
        
    ax.plot(X, Y1, color = 'tab:blue', label = 'north', alpha = 0.7)
    ax.plot(X, Y2, color = 'tab:orange', label = 'east', alpha = 0.7)
    ax.plot(X, Y3, color = 'tab:green', label = 'south', alpha = 0.7)
    ax.plot(X, Y4, color = 'tab:red', label = 'west', alpha = 0.7)

    if use_last == True and extra_exists:
        ax.plot(X, Y_E, color = 'tab:purple', label = extra_label, alpha = 0.7)
        
    ax.legend()

    ax.plot([X[0], X[-1]], [0., 0.], color = 'black')
    
    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label, loc = "right")
    
    # raw_llim = np.nanmin([Y1[x_lbin:x_ubin], Y2[x_lbin:x_ubin],
    #                       Y3[x_lbin:x_ubin], Y4[x_lbin:x_ubin]])
    raw_ulim = np.nanmax([Y1[x_lbin:x_ubin], Y2[x_lbin:x_ubin],
                          Y3[x_lbin:x_ubin], Y4[x_lbin:x_ubin]])
    
    # if raw_llim > 0.:
    #     ax.set_ylim([0, 1.1 * raw_ulim])
    # else:
    #     ax.set_ylim([1.1 * raw_llim, 1.1 * raw_ulim])
    ax.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])

    ax.set_ylabel(y_label_1)
    ax.ticklabel_format(axis = 'y', useMathText = True, style='sci', scilimits=(-1,1))
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.))
    
    ax.grid(which = 'both')

    # Subplot: Raw Signals - Far range
    
    ax2.plot(X, Y1, color = 'tab:blue', alpha = 0.7)
    ax2.plot(X, Y2, color = 'tab:orange', alpha = 0.7)
    ax2.plot(X, Y3, color = 'tab:green', alpha = 0.7)
    ax2.plot(X, Y4, color = 'tab:red', alpha = 0.7)

    if use_last == True and extra_exists:
        ax2.plot(X, Y_E, color = 'tab:purple', alpha = 0.7)
                
    ax2.plot([X[0], X[-1]], [0., 0.], color = 'black')

    ax2.set_xlim([x_ulim, 20.])
    # ax2.set_xlabel(x_label)
    ax2.xaxis.set_minor_locator(MultipleLocator(2.))

    # if raw_llim > 0.:
    #     ax2.set_ylim([0, 1.1 * raw_ulim])
    # else:
    #     ax2.set_ylim([1.1 * raw_llim, 1.1 * raw_ulim])
    ax2.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])

    ax2.set(yticklabels=[])

    ax2.grid(which = 'both')
    
    # Subplot: Normalized Signals - Near range
        
    ax3.plot(X, Y1_N, color = 'tab:blue')
    ax3.plot(X, Y2_N, color = 'tab:orange')
    ax3.plot(X, Y3_N, color = 'tab:green')
    ax3.plot(X, Y4_N, color = 'tab:red')
    
    if use_last == True and extra_exists:
        ax3.plot(X, Y_E_N, color = 'purple')
        
    ax3.fill_between(X, Y1_NL, Y1_NU, color = 'tab:blue', alpha = 0.3)
    ax3.fill_between(X, Y2_NL, Y2_NU, color = 'tab:orange', alpha = 0.3)
    ax3.fill_between(X, Y3_NL, Y3_NU, color = 'tab:green', alpha = 0.3)
    ax3.fill_between(X, Y4_NL, Y4_NU, color = 'tab:red', alpha = 0.3)

    ax3.plot([X[0], X[-1]], [0., 0.], color = 'black')

    ax3.set_xticks(x_ticks, labels = x_ticks)
    ax3.set_xlim([x_llim, x_ulim])
    ax3.set_xlabel(x_label, loc = "right")

    ax3.set_ylim([y_llim_nr, y_ulim_nr])
    ax3.set_ylabel('Normalized RC Signals [A.U.]')
    ax3.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.))

    ax3.grid(which = 'both')

    ax3.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')

    n_llim = np.round(norm_region[0], decimals = 2)
    n_ulim = np.round(norm_region[1], decimals = 2)
    
    ax3.text(0.30 * x_ulim, 0.90 * y_ulim_nr, 
             f'norm. region: {n_llim} - {n_ulim} km',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 9))

    # Subplot: Normalized Signals - Far Range
        
    ax4.plot(X, Y1_N, color = 'tab:blue')
    ax4.plot(X, Y2_N, color = 'tab:orange')
    ax4.plot(X, Y3_N, color = 'tab:green')
    ax4.plot(X, Y4_N, color = 'tab:red')

    if use_last == True and extra_exists:
        ax4.plot(X, Y_E_N, color = 'purple')
        
    ax4.fill_between(X, Y1_NL, Y1_NU, color = 'tab:blue', alpha = 0.3)
    ax4.fill_between(X, Y2_NL, Y2_NU, color = 'tab:orange', alpha = 0.3)
    ax4.fill_between(X, Y3_NL, Y3_NU, color = 'tab:green', alpha = 0.3)
    ax4.fill_between(X, Y4_NL, Y4_NU, color = 'tab:red', alpha = 0.3)

    ax4.plot([X[0], X[-1]], [0., 0.], color = 'black')

    ax4.set_xlim([x_ulim, 20.])
    # ax4.set_xlabel(x_label)
    ax4.xaxis.set_minor_locator(MultipleLocator(2.))

    ax4.set_ylim([y_llim_nr, y_ulim_nr])
    ax4.set(yticklabels=[])

    ax4.grid(which = 'both')

    # Subplot: Normalized Deviations
    ax5.plot(X, (Y1_N - Y_NM) / Y_NM, color = 'tab:blue', label='_nolegend_')
    ax5.plot(X, (Y2_N - Y_NM) / Y_NM, color = 'tab:orange', label='_nolegend_')
    ax5.plot(X, (Y3_N - Y_NM) / Y_NM, color = 'tab:green', label='_nolegend_')
    ax5.plot(X, (Y4_N - Y_NM) / Y_NM, color = 'tab:red', label='_nolegend_')
    ax5.plot(X, Y_RMSE, color = 'yellow', label = 'RMS Sector Diff.')
    
    if use_last == True and extra_exists:
        ax5.plot(X, Y_DIFF, color = 'tab:purple', label = f'{extra_label_short} - {extra_label_short[0]}')
    
    ax5.legend()

    ax5.fill_between(X, (Y1_NL - Y_NM) / Y_NM, (Y1_NU - Y_NM) / Y_NM, 
                     color = 'tab:blue', alpha = 0.3)
    ax5.fill_between(X, (Y2_NL - Y_NM) / Y_NM, (Y2_NU - Y_NM) / Y_NM, 
                     color = 'tab:orange', alpha = 0.3)
    ax5.fill_between(X, (Y3_NL - Y_NM) / Y_NM, (Y3_NU - Y_NM) / Y_NM, 
                     color = 'tab:green', alpha = 0.3)
    ax5.fill_between(X, (Y4_NL - Y_NM) / Y_NM, (Y4_NU - Y_NM) / Y_NM, 
                     color = 'tab:red', alpha = 0.3)

    ax5.plot(X, np.zeros(X.shape), '--', color = 'black', zorder = 10)

    ax5.plot(X, -0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 10, alpha = 0.7)
    ax5.plot(X, 0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 10, alpha = 0.7)

    ax5.set_xticks(x_ticks, labels = x_ticks)
    ax5.set_xlim([x_llim, x_ulim])
    ax5.set_xlabel(x_label)
    
    y_ticks = np.round(np.arange(-0.20, 0.20 + 0.05, 0.05), decimals = 2)
    ax5.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax5.set_ylim([y_ticks[0], y_ticks[-1]])
    ax5.set_ylabel('Relative Sector Deviation')
    ax5.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.))

    ax5.grid(which = 'both')

    fpath = os.path.join(dir_out, 'plots', fname)
   
    fig.savefig(fpath, dpi = dpi_val)

    fig.clf()
    
    plt.close()
    
    plt.rcParams.update(plt.rcParamsDefault)

    if color_reduction == True:
        run(["convert", fpath, "-posterize", "8", fpath])
        
    return(fpath)


def telecover_rin(dir_out, fname, title, dpi_val, color_reduction, 
                  norm_region, use_nonrc, x_vals, 
                  y1_raw, y2_raw, 
                  y1_vals, y2_vals,
                  y1_extr, y2_extr,
                  y1_extr_raw, y2_extr_raw,
                  y1_lvar, y2_lvar, 
                  y1_uvar, y2_uvar,
                  coef_1, coef_2, 
                  coef_extra_1, coef_extra_2, extra_rin, ranges,
                  x_lbin, x_ubin, x_llim, x_ulim, 
                  y_llim, y_ulim, y_llim_nr, y_ulim_nr, 
                  x_label, x_tick, use_last, iters):

    # Create the variables to be plotted X, Y
    X = x_vals
    
    if use_nonrc == True:
        R = ranges
    
        Y1 = y1_raw / np.power(R, 2)
        Y2 = y2_raw / np.power(R, 2)
        
        y_label_1 = 'Non RC Signals - Raw Units'
    
    else:    
        Y1 = y1_raw
        Y2 = y2_raw
        
        y_label_1 = 'RC Signals [A.U.]'

    
    Y1_N = coef_1 * y1_vals
    Y2_N = coef_2 * y2_vals

    Y1_NL = coef_1 * y1_lvar
    Y2_NL = coef_2 * y2_lvar
    
    Y1_NU = coef_1 * y1_uvar
    Y2_NU = coef_2 * y2_uvar

    Y_NM = np.mean([Y1_N, Y2_N], axis = 0) 
    
    Y_E = np.nan * Y_NM
    Y_O = np.nan * Y_NM

    if extra_rin['O']:
        Y_E = y1_extr_raw
        Y_E_N = coef_extra_1 * y1_extr
        Y_O = Y1_N
        extra_label = f'outer{iters+1}'
        extra_label_short = f'O{iters+1}'
        extra_exists = True
    elif extra_rin['I']:
        Y_E = y2_extr_raw
        Y_E_N = coef_extra_2 * y2_extr
        Y_O = Y2_N
        extra_label = f'inner{iters+1}'
        extra_label_short = f'I{iters+1}'
        extra_exists = True
    else:
        extra_exists = False

    Y_DIFF = (Y_E_N - Y_O) / Y_O
    
    x_ticks = np.arange(x_tick * np.ceil(x_llim / x_tick), 
                        x_tick * (np.floor(x_ulim / x_tick) + 1.), 
                        x_tick)
        
    if x_tick >= x_ulim - x_llim:
        raise Exception(f"The x_tick ({x_tick}) must be smaller than the width of the normalization_region ({norm_region}) for the telecover test. Please revise the settings_file.ini ")
    
    if x_llim > norm_region[1] or x_ulim < norm_region[0]:
        raise Exception(f"The normalization_region ({norm_region}) for the telecover fit is out of the provided x_lims ([{x_llim}, {x_ulim}]). Please revise the settings_file.ini ")
    
    if np.abs(x_llim - x_ticks[0]) < x_tick * 0.25:
        x_ticks[0] = x_llim
    else:
        x_ticks = np.hstack((x_llim, x_ticks))

    if np.abs(x_ulim - x_ticks[-1]) < x_tick * 0.25:
        x_ticks[-1] = x_ulim
    else:
        x_ticks = np.hstack((x_ticks, x_ulim))
    
    x_ticks = np.round(x_ticks, decimals = 2)

    plt.rc('font', size = 14) 

    # Create the figure
    fig = plt.figure(figsize=(20. , 4.))
    
    fig.suptitle(title)

    # Subplot: Raw Signals - Near range
    ax = fig.add_axes([0.04,0.13,0.19,0.7])
    ax2 = fig.add_axes([0.24,0.13,0.07,0.7])
    ax3 = fig.add_axes([0.36,0.13,0.19,0.7])
    ax4 = fig.add_axes([0.56,0.13,0.07,0.7])
    ax5 = fig.add_axes([0.69,0.13,0.29,0.7])
        
    ax.plot(X, Y1, color = 'tab:green', label = 'outer', alpha = 0.7)
    ax.plot(X, Y2, color = 'tab:orange', label = 'inner', alpha = 0.7)

    if use_last == True and extra_exists:
        ax.plot(X, Y_E, color = 'tab:purple', label = extra_label, alpha = 0.7)
        
    ax.legend()
    
    ax.plot([X[0], X[-1]], [0., 0.], color = 'black')

    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label, loc = "right")
    
    # raw_llim = np.nanmin([Y1[x_lbin:x_ubin], Y2[x_lbin:x_ubin]])
    raw_ulim = np.nanmax([Y1[x_lbin:x_ubin], Y2[x_lbin:x_ubin]])
    
    # if raw_llim > 0.:
    #     ax.set_ylim([0, 1.1 * raw_ulim])
    # else:
    #     ax.set_ylim([1.1 * raw_llim, 1.1 * raw_ulim])
    ax.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])

    ax.set_ylabel(y_label_1)
    ax.ticklabel_format(axis = 'y', useMathText = True, style='sci', scilimits=(-1,1))
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.))
    
    ax.grid(which = 'both')

    # Subplot: Raw Signals - Far range
    
    ax2.plot(X, Y1, color = 'tab:green', alpha = 0.7)
    ax2.plot(X, Y2, color = 'tab:orange', alpha = 0.7)

    if use_last == True and extra_exists:
        ax2.plot(X, Y_E, color = 'tab:purple', alpha = 0.7)
                
    ax2.plot([X[0], X[-1]], [0., 0.], color = 'black')
    
    ax2.set_xlim([x_ulim, 20.])
    # ax2.set_xlabel(x_label)
    ax2.xaxis.set_minor_locator(MultipleLocator(2.))

    # if raw_llim > 0.:
    #     ax2.set_ylim([0, 1.1 * raw_ulim])
    # else:
    #     ax2.set_ylim([1.1 * raw_llim, 1.1 * raw_ulim])

    ax2.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])

    ax2.set(yticklabels=[])

    ax2.grid(which = 'both')
    
    # Subplot: Normalized Signals - Near range
        
    ax3.plot(X, Y1_N, color = 'tab:green')
    ax3.plot(X, Y2_N, color = 'tab:orange')
    
    if use_last == True and extra_exists:
        ax3.plot(X, Y_E_N, color = 'purple')
        
    ax3.fill_between(X, Y1_NL, Y1_NU, color = 'tab:green', alpha = 0.3)
    ax3.fill_between(X, Y2_NL, Y2_NU, color = 'tab:orange', alpha = 0.3)

    ax3.plot([X[0], X[-1]], [0., 0.], color = 'black')

    ax3.set_xticks(x_ticks, labels = x_ticks)
    ax3.set_xlim([x_llim, x_ulim])
    ax3.set_xlabel(x_label, loc = "right")

    ax3.set_ylim([y_llim_nr, y_ulim_nr])
    ax3.set_ylabel('Normalized RC Signals [A.U.]')
    ax3.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.))

    ax3.grid(which = 'both')

    ax3.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')

    n_llim = np.round(norm_region[0], decimals = 2)
    n_ulim = np.round(norm_region[1], decimals = 2)
    
    ax3.text(0.30 * x_ulim, 0.90 * y_ulim_nr, 
             f'norm. region: {n_llim} - {n_ulim} km',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 9))

    # Subplot: Normalized Signals - Far Range
        
    ax4.plot(X, Y1_N, color = 'tab:green')
    ax4.plot(X, Y2_N, color = 'tab:orange')

    if use_last == True:
        ax4.plot(X, Y_E_N, color = 'purple')
        
    ax4.fill_between(X, Y1_NL, Y1_NU, color = 'tab:green', alpha = 0.3)
    ax4.fill_between(X, Y2_NL, Y2_NU, color = 'tab:orange', alpha = 0.3)
    
    ax4.plot([X[0], X[-1]], [0., 0.], color = 'black')

    ax4.set_xlim([x_ulim, 20.])
    # ax4.set_xlabel(x_label)
    ax4.xaxis.set_minor_locator(MultipleLocator(2.))

    ax4.set_ylim([y_llim_nr, y_ulim_nr])
    ax4.set(yticklabels=[])

    ax4.grid(which = 'both')

    # Subplot: Normalized Deviations
    ax5.plot(X, (Y1_N - Y_NM) / Y_NM, color = 'tab:green', label='_nolegend_')
    ax5.plot(X, (Y2_N - Y_NM) / Y_NM, color = 'tab:orange', label='_nolegend_')

    if use_last == True and extra_exists:
        ax5.plot(X, Y_DIFF, color = 'tab:purple', label = f'{extra_label_short}-{extra_label_short[0]}')
    
    ax5.legend()

    ax5.fill_between(X, (Y1_NL - Y_NM) / Y_NM, (Y1_NU - Y_NM) / Y_NM, 
                     color = 'tab:green', alpha = 0.3)
    ax5.fill_between(X, (Y2_NL - Y_NM) / Y_NM, (Y2_NU - Y_NM) / Y_NM, 
                     color = 'tab:orange', alpha = 0.3)

    ax5.plot(X, np.zeros(X.shape), '--', color = 'black', zorder = 10)

    ax5.plot(X, -0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 10, alpha = 0.7)
    ax5.plot(X, 0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 10, alpha = 0.7)

    ax5.set_xticks(x_ticks, labels = x_ticks)
    ax5.set_xlim([x_llim, x_ulim])
    ax5.set_xlabel(x_label)
    
    y_ticks = np.round(np.arange(-0.20, 0.20 + 0.05, 0.05), decimals = 2)
    ax5.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax5.set_ylim([y_ticks[0], y_ticks[-1]])
    ax5.set_ylabel('Relative Sector Deviation')
    ax5.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.))

    ax5.grid(which = 'both')

    fpath = os.path.join(dir_out, 'plots', fname)
   
    fig.savefig(fpath, dpi = dpi_val)

    fig.clf()
    
    plt.close()
    
    plt.rcParams.update(plt.rcParamsDefault)

    if color_reduction == True:
        run(["convert", fpath, "-posterize", "8", fpath])
        
    return(fpath)

def intercomparison(dir_out, fname, title, dpi_val, color_reduction, 
                    use_lin, x_refr, refr_hwin,
                    x_vals, y1_vals, y2_vals, y1_errs, y2_errs, y3_vals, 
                    coef1, coef2, use_molecular, 
                    x_lbin, x_ubin, x_llim, x_ulim, y_llim, y_ulim, 
                    x_label, y_label, x_tick, lidars):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]
    Y1 = coef1 * y1_vals[slice(x_lbin, x_ubin)]
    Y2 = coef2 * y2_vals[slice(x_lbin, x_ubin)]
    Y1E = coef1 * y1_errs[slice(x_lbin, x_ubin)]
    Y2E = coef2 * y2_errs[slice(x_lbin, x_ubin)]
    Y3 = y3_vals[slice(x_lbin, x_ubin)]
    
    # Create the figure
    fig = plt.figure(figsize=(12. , 4.))
    fig.suptitle(title)

    ax = fig.add_axes([0.08,0.13,0.50,0.7])
    
    # ax.set_title(title, pad = 15)
    
    ax.plot(X, Y1, color = 'tab:blue', label = lidars[0])
    ax.plot(X, Y2, color = 'tab:red', label = lidars[1])
    if use_molecular:
        ax.plot(X, Y3, color = 'tab:green', label = 'molecular', zorder = 0)

    x_ticks = np.arange(x_tick * np.ceil(x_llim / x_tick), 
                        x_tick * (np.floor(x_ulim / x_tick) + 1.), 
                        x_tick)
        
    if np.abs(x_llim - x_ticks[0]) < x_tick * 0.25:
        x_ticks[0] = x_llim
    else:
        x_ticks = np.hstack((x_llim, x_ticks))

    if np.abs(x_ulim - x_ticks[-1]) < x_tick * 0.25:
        x_ticks[-1] = x_ulim
    else:
        x_ticks = np.hstack((x_ticks, x_ulim))
        
    x_ticks = np.round(x_ticks, decimals = 2)

    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label)
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.))

    ax.set_ylim([y_llim, y_ulim])
    if use_molecular:
        ax.set_ylabel(y_label)
        
    else:
        ax.set_ylabel('Attenuated Backscatter [A.U.]')

    if use_lin == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    ax.legend()

    ax.fill_between(X, Y1 - Y1E, Y1 + Y1E, color = 'tab:blue', alpha = 0.3)
    ax.fill_between(X, Y2 - Y2E, Y2 + Y2E, color = 'tab:red', alpha = 0.3)
    
    x_refr_bin = np.where(x_vals >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))

    ax.axvspan(x_vals[x_refr_bin - refr_hbin], x_vals[x_refr_bin + refr_hbin + 1],
               alpha = 0.2, facecolor = 'tab:grey')

    n_llim = np.round(x_refr - refr_hwin * 1e-3, decimals = 2)
    n_ulim = np.round(x_refr + refr_hwin * 1e-3, decimals = 2)
    
    if use_lin == False:
        ax.text(0.55 * x_ulim, 0.60 * y_ulim, 
                f'norm. region: {n_llim} - {n_ulim} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
    else:
        ax.text(0.55 * x_ulim, 0.9 * y_ulim, 
                f'norm. region: {n_llim} - {n_ulim} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
        
    ax2 = fig.add_axes([0.65,0.13,0.30,0.7])

    # ax2.set_title('', pad = 15)

    frac_bias = 2. * (Y1 - Y2) / (Y1 + Y2)
    ax2.plot(X, frac_bias, color = 'tab:blue')
    
    x_tick_2 = 2. * x_tick
    x_ticks_2 = np.arange(x_tick_2 * np.ceil(x_llim / x_tick_2), 
                          x_tick_2 * (np.floor(x_ulim / x_tick_2) + 1.), 
                          x_tick_2)
        
    if np.abs(x_llim - x_ticks_2[0]) < x_tick_2 * 0.25:
        x_ticks_2[0] = x_llim
    else:
        x_ticks_2 = np.hstack((x_llim, x_ticks_2))

    if np.abs(x_ulim - x_ticks[-1]) < x_tick_2 * 0.25:
        x_ticks_2[-1] = x_ulim
    else:
        x_ticks_2 = np.hstack((x_ticks_2, x_ulim))
        
    x_ticks_2 = np.round(x_ticks_2, decimals = 2)
    
    ax2.set_xticks(x_ticks_2, labels = x_ticks_2)
    ax2.set_xlim([x_llim, x_ulim])
    ax2.set_xlabel(x_label)
    ax2.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.))

    y_ticks = np.round(np.arange(-0.40, 0.40 + 0.10, 0.10), decimals = 2)
    ax2.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax2.set_ylim([y_ticks[0], y_ticks[-1]])
    ax2.set_ylabel(f'Fractional Bias ({lidars[0]} - {lidars[1]})')
    
    ax2.grid(which = 'both')
    
    ax2.axhline(c = 'k')

    ax2.axvspan(x_vals[x_refr_bin - refr_hbin], x_vals[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')
    
        
    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()

    plt.close()
    
    if color_reduction == True:
        run(["convert", fpath, "-posterize", "8", fpath])
    
    return(fpath)

def polarization_calibration(dir_out, fname, title, dpi_val, color_reduction, 
                             cal_region, vdr_region,
                             x_vals_cal, x_vals_vdr,
                             y1_vals, y2_vals, y3_vals, y4_vals, y5_vals, y6_vals,
                             eta, eta_f_s, eta_s,
                             delta_m, delta_c_def, delta_c, epsilon,
                             eta_err, eta_f_s_err, eta_s_err,
                             delta_c_def_err, delta_c_err, epsilon_err,
                             x_lbin_cal, x_ubin_cal, 
                             x_llim_cal, x_ulim_cal, 
                             y_llim_cal, y_ulim_cal, 
                             x_lbin_vdr, x_ubin_vdr, 
                             x_llim_vdr, x_ulim_vdr, 
                             y_llim_vdr, y_ulim_vdr, 
                             K, G_R, G_T, H_R, H_T,
                             y_label_cal, x_label_cal, x_tick_cal,
                             y_label_vdr, x_label_vdr, x_tick_vdr):
        
    delta_l = (delta_c - delta_m) / (1. - delta_c * delta_m)
    delta_l_err = delta_c_err * (1. - delta_m) * (1. + delta_c) / \
        (1. - delta_m * delta_c)**2
    
    # Create the variables to be plotted X, Y
    XA = x_vals_cal[slice(x_lbin_cal, x_ubin_cal)]
    XB = x_vals_vdr[slice(x_lbin_vdr, x_ubin_vdr)]
    
    Y1 = y1_vals[slice(x_lbin_cal, x_ubin_cal)]
    Y2 = y2_vals[slice(x_lbin_cal, x_ubin_cal)]
    Y3 = y3_vals[slice(x_lbin_cal, x_ubin_cal)]
    Y4 = y4_vals[slice(x_lbin_vdr, x_ubin_vdr)]
    Y5 = y5_vals[slice(x_lbin_vdr, x_ubin_vdr)]
    Y6 = y6_vals[slice(x_lbin_vdr, x_ubin_vdr)]

    Y1E = np.nan * Y1
    Y2E = np.nan * Y2
    Y3E = np.nan * Y3
    Y4E = np.nan * Y4
    Y5E = np.nan * Y5
    
    # Create the figure
    fig = plt.figure(figsize=(12. , 4.))
    fig.suptitle(title)

    ax = fig.add_axes([0.07,0.13,0.40,0.65])
        
    ax.plot(XA, Y1, color = 'tab:purple', label = 'η')
    ax.plot(XA, Y2, color = 'tab:red', label = '$η_{+45}$')
    ax.plot(XA, Y3, color = 'tab:cyan', label = '$η_{-45}$')

    ax.fill_between(XA, Y1 - Y1E, Y1 + Y1E, color = 'tab:purple', alpha = 0.3)
    ax.fill_between(XA, Y2 - Y2E, Y2 + Y2E, color = 'tab:red', alpha = 0.3)
    ax.fill_between(XA, Y3 - Y3E, Y3 + Y3E, color = 'tab:cyan', alpha = 0.3)
    
    if x_tick_cal >= x_ulim_cal - x_llim_cal:
        raise Exception(f"The x_tick_calibration ({x_tick_cal}) must be smaller than the width of the calibration_region ({cal_region}) for the polarization calibration test. Please revise the settings_file.ini ")
    
    if x_llim_cal > cal_region[1] or x_ulim_cal < cal_region[0]:
        raise Exception(f"The calibration_region ({cal_region}) is out of the provided x_lims_calibration ([{x_llim_cal}, {x_ulim_cal}]) for the polarization calibration test. Please revise the settings_file.ini ")

    x_ticks_cal = np.arange(x_tick_cal * np.ceil(x_llim_cal / x_tick_cal), 
                            x_tick_cal * (np.floor(x_ulim_cal / x_tick_cal) + 1.), 
                            x_tick_cal)
        
    if np.abs(x_llim_cal - x_ticks_cal[0]) < x_tick_cal * 0.25:
        x_ticks_cal[0] = x_llim_cal
    else:
        x_ticks_cal = np.hstack((x_llim_cal, x_ticks_cal))

    if np.abs(x_ulim_cal - x_ticks_cal[-1]) < x_tick_cal * 0.25:
        x_ticks_cal[-1] = x_ulim_cal
    else:
        x_ticks_cal = np.hstack((x_ticks_cal, x_ulim_cal))

    x_ticks_cal = np.round(x_ticks_cal, decimals = 2)

    ax.set_xticks(x_ticks_cal, labels = x_ticks_cal)
    ax.set_xlim([x_llim_cal, x_ulim_cal])
    ax.set_xlabel(x_label_cal)

    ax.set_ylim([y_llim_cal, y_ulim_cal])
    ax.set_ylabel(y_label_cal)

    ax.grid(which = 'both')
    ax.legend(loc = 'upper right')

    ax.axvspan(cal_region[0], cal_region[1],
               alpha = 0.2, facecolor = 'tab:grey')

    c_llim = np.round(cal_region[0], decimals = 2)
    c_ulim = np.round(cal_region[1], decimals = 2)
    
    ax.text(0.05 * x_ulim_cal, 0.95 * y_ulim_cal, 
            f'cal. region: {c_llim} - {c_ulim} km',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))  
    # ax.text(0.55 * x_ulim_cal, 0.87 * y_ulim_cal, 
    #         f'ε: {round_it(epsilon,2)}'+'${}^o \pm$ '+f'{round_it(epsilon_err,2)}'+'${}^o$,',
    #         bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.05 * x_ulim_cal, 0.87 * y_ulim_cal, 
            f'ε: {round_it(epsilon,2)}' +'${}^o$'+ ' $\pm$ ' + f'{round_it(epsilon_err, 2)}' +'${}^o$' + f', K: {round_it(K, 4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.05 * x_ulim_cal, 0.79 * y_ulim_cal, 
            r'$η^{\star}_{f}$'+f': {round_it(eta_f_s, 3)}' + ' $\pm$ ' + f'{round_it(eta_f_s_err, 2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.05 * x_ulim_cal, 0.71 * y_ulim_cal, 
            r'$η^{\star}$'+f': {round_it(eta_s, 3)}' + ' $\pm$ ' + f'{round_it(eta_s_err, 2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    if eta_err / eta <= 0.02:
        color_eta = 'tab:green'
    else:
        color_eta = 'tab:red'

    ax.text(0.05 * x_ulim_cal, 0.63 * y_ulim_cal, 
            f'η: {round_it(eta, 3)}' + ' $\pm$ ' + f'{round_it(eta_err, 2)}',
            bbox = dict(facecolor = color_eta, alpha = 0.1, zorder = 3)) 

    ax2 = fig.add_axes([0.55,0.13,0.40,0.65])

    ax2.plot(XB, Y4, color = 'tab:blue', label = 'measured')
    ax2.plot(XB, Y5, color = 'tab:orange', label = 'corrected')
    ax2.plot(XB, Y6, color = 'tab:green', label = 'molecular')

    ax2.fill_between(XB, Y4 - Y4E, Y4 + Y4E, color = 'tab:blue', alpha = 0.3)
    ax2.fill_between(XB, Y5 - Y5E, Y5 + Y5E, color = 'tab:orange', alpha = 0.3)
    
    x_ticks_vdr = np.arange(x_tick_vdr * np.ceil(x_llim_vdr / x_tick_vdr), 
                            x_tick_vdr * (np.floor(x_ulim_vdr / x_tick_vdr) + 1.), 
                            x_tick_vdr)

    if x_tick_vdr >= x_ulim_vdr - x_llim_vdr:
        raise Exception(f"The x_tick_rayleigh ({x_tick_vdr}) must be smaller than the width of the rayleigh_region ({vdr_region}). Please revise the settings_file.ini ")
    
    if x_llim_vdr > vdr_region[1] or x_ulim_vdr < vdr_region[0]:
        raise Exception(f"The rayleigh_region ({vdr_region}) is out of the provided x_lims_rayleigh ([{x_llim_vdr}, {x_ulim_vdr}]). Please revise the settings_file.ini ")
        
    if np.abs(x_llim_vdr - x_ticks_vdr[0]) < x_tick_vdr * 0.25:
        x_ticks_vdr[0] = x_llim_vdr
    else:
        x_ticks_vdr = np.hstack((x_llim_vdr, x_ticks_vdr))

    if np.abs(x_ulim_vdr - x_ticks_vdr[-1]) < x_tick_vdr * 0.25:
        x_ticks_vdr[-1] = x_ulim_vdr
    else:
        x_ticks_vdr = np.hstack((x_ticks_vdr, x_ulim_vdr))

    x_ticks_vdr = np.round(x_ticks_vdr, decimals = 2)

    ax2.set_xticks(x_ticks_vdr, labels = x_ticks_vdr)
    ax2.set_xlim([x_llim_vdr, x_ulim_vdr])
    ax2.set_xlabel(x_label_vdr)

    ax2.set_ylim([y_llim_vdr, y_ulim_vdr])
    ax2.set_ylabel(y_label_vdr)

    ax2.grid(which = 'both')
    ax2.legend(loc = 'upper right')

    ax2.axvspan(vdr_region[0], vdr_region[1],
                alpha = 0.2, facecolor = 'tab:grey')

    m_llim = np.round(vdr_region[0], decimals = 2)
    m_ulim = np.round(vdr_region[1], decimals = 2)

    ax2.text(0.05 * x_ulim_vdr, 0.65 * y_ulim_vdr, 
             f'mol. cal. region: {m_llim} - {m_ulim} km',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))  
    ax2.text(0.05 * x_ulim_vdr, 0.55 * y_ulim_vdr, 
            r'$δ_m$: '+f'{np.round(delta_m,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.45 * y_ulim_vdr, 
            r'$δ^{\star}$'+f': {np.round(delta_c_def,4)}' + ' $\pm$ ' + f'{round_it(delta_c_def_err,2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.35 * y_ulim_vdr, 
            r'$δ_{c}$'+f': {np.round(delta_c,4)}' + ' $\pm$ ' + f'{round_it(delta_c_err,2)}', 
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.25 * y_ulim_vdr, 
            r'$δ_{res}$: '+f'{np.round(delta_l,4)}' + ' $\pm$ ' + f'{round_it(delta_l_err,2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.15 * y_ulim_vdr, 
            r'$G_R$: '+f'{np.round(G_R,4)}, $G_T$: '+f'{round_it(G_T,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.05 * y_ulim_vdr, 
            r'$H_R$: '+f'{np.round(H_R,4)}, $H_T$: '+f'{round_it(H_T,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 

    fpath = os.path.join(dir_out, 'plots', fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    if color_reduction == True:
        run(["convert", fpath, "-posterize", "8", fpath])
    
    return(fpath)

def round_it(x, sig):
    
    if not np.isfinite(x) or np.isnan(x):
        x = -999.
        sig = 3
        
    if x != 0.:
        x_out = np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)
    else:
        x_out = 0.
        
    return x_out