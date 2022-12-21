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
from matplotlib.ticker import MultipleLocator

def quicklook(dir_out, fname, title, dpi_val, use_log,
              x_vals, t_vals, y_vals, z_vals, 
              x_lbin, x_ubin, y_lbin, y_ubin, 
              y_llim, y_ulim, z_llim, z_ulim,
              x_label, t_label, y_label, z_label, 
              x_tick, t_tick, y_tick, nodes):
        
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
    ax = fig.add_axes([0.09,0.13,0.99,0.7])
    
    ax.set_title(title, pad = 15)
    
    qck = ax.pcolormesh(X, Y, Z, vmin = z_llim, vmax = z_ulim, cmap = my_cmap)
    fig.colorbar(qck, label = z_label, extend = 'both')
    
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
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval = int(t_tick)))

    plt.xticks(rotation = 35)

    if len(nodes) == 0:
        ax1 = ax.twiny()
    
        x_ticks = np.arange(x_lbin, x_ubin + x_tick, x_tick).astype(int)
        ax1.set_xticks(x_ticks, labels = x_ticks)
        
    fpath = os.path.join(dir_out, 'plots', fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)

def rayleigh(dir_out, fname, title, dpi_val, use_lin, x_refr, refr_hwin,
             x_vals, y1_vals, y2_vals, y1_errs, coef, rsem,
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

    ax.fill_between(X, Y1 - Y1E, Y1 + Y1E, color = 'tab:blue', alpha = 0.3)
    
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
    ax.set_ylabel(y_label)
    if use_lin == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    ax.legend(loc = 'lower left')

    x_refr_bin = np.where(x_vals >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))

    ax.axvspan(x_vals[x_refr_bin - refr_hbin], x_vals[x_refr_bin + refr_hbin + 1],
               alpha = 0.2, facecolor = 'tab:grey')

    n_llim = np.round(x_refr - refr_hwin * 1e-3, decimals = 2)
    n_ulim = np.round(x_refr + refr_hwin * 1e-3, decimals = 2)
    
    if use_lin == False:

        ax.text(0.55 * x_ulim, 0.60 * y_ulim, 
                f'norm: {n_llim} - {n_ulim} Km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))

        ax.text(0.55 * x_ulim, 0.30 * y_ulim, 
                 f'rsem: {np.round(rsem, decimals = 4)}',
                 bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
    else:
        ax.text(0.55 * x_ulim, 0.9 * y_ulim, 
                f'norm: {n_llim} - {n_ulim} Km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))

        ax.text(0.55 * x_ulim, 0.82 * y_ulim, 
                 f'rsem: {np.round(rsem, decimals = 4)}',
                 bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))

    ax2 = fig.add_axes([0.65,0.13,0.30,0.7])

    ax2.fill_between(X, (Y1 - Y1E - Y2) / Y2, 
                     (Y1 + Y1E - Y2) / Y2, color = 'tab:blue', 
                     alpha = 0.3, label = 'sem')
    ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue',label = 'mean')
    
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
    ax2.legend(loc = 'lower left')
    
    ax2.axvspan(x_vals[x_refr_bin - refr_hbin], x_vals[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')
    
    fpath = os.path.join(dir_out, 'plots', fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)

def telecover_sec(dir_out, fname, title, dpi_val, 
                  x_refr, refr_hwin, use_nonrc,x_vals, 
                  y1_raw, y2_raw, y3_raw, y4_raw,
                  y1_lraw, y2_lraw, y3_lraw, y4_lraw,
                  y1_uraw, y2_uraw, y3_uraw, y4_uraw,
                  y1_vals, y2_vals, y3_vals, y4_vals,
                  y1_extr, y2_extr, y3_extr, y4_extr,
                  y1_lvar, y2_lvar, y3_lvar, y4_lvar, 
                  y1_uvar, y2_uvar, y3_uvar, y4_uvar, 
                  coef_1, coef_2, coef_3, coef_4, ranges,
                  x_lbin, x_ubin, x_llim, x_ulim, 
                  y_llim, y_ulim, y_llim_nr, y_ulim_nr, 
                  x_label, x_tick, use_last):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]
    
    if use_nonrc == True:
        R = ranges[slice(x_lbin, x_ubin)]
    
        Y1 = y1_raw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y2 = y2_raw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y3 = y3_raw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y4 = y4_raw[slice(x_lbin, x_ubin)] / np.power(R, 2)
    
        Y1_L = y1_lraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y2_L = y2_lraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y3_L = y3_lraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y4_L = y4_lraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        
        Y1_U = y1_uraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y2_U = y2_uraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y3_U = y3_uraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y4_U = y4_uraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        
        y_label_1 = 'Non RC Signals - Raw Units'
    
    else:    
        Y1 = y1_raw[slice(x_lbin, x_ubin)]
        Y2 = y2_raw[slice(x_lbin, x_ubin)]
        Y3 = y3_raw[slice(x_lbin, x_ubin)]
        Y4 = y4_raw[slice(x_lbin, x_ubin)]
    
        Y1_L = y1_lraw[slice(x_lbin, x_ubin)]
        Y2_L = y2_lraw[slice(x_lbin, x_ubin)]
        Y3_L = y3_lraw[slice(x_lbin, x_ubin)]
        Y4_L = y4_lraw[slice(x_lbin, x_ubin)]
        
        Y1_U = y1_uraw[slice(x_lbin, x_ubin)]
        Y2_U = y2_uraw[slice(x_lbin, x_ubin)]
        Y3_U = y3_uraw[slice(x_lbin, x_ubin)]
        Y4_U = y4_uraw[slice(x_lbin, x_ubin)]
        
        y_label_1 = 'RC Signals [A.U.]'

    
    Y1_N = coef_1 * y1_vals[slice(x_lbin, x_ubin)]
    Y2_N = coef_2 * y2_vals[slice(x_lbin, x_ubin)]
    Y3_N = coef_3 * y3_vals[slice(x_lbin, x_ubin)]
    Y4_N = coef_4 * y4_vals[slice(x_lbin, x_ubin)]

    Y1_NL = coef_1 * y1_lvar[slice(x_lbin, x_ubin)]
    Y2_NL = coef_2 * y2_lvar[slice(x_lbin, x_ubin)]
    Y3_NL = coef_3 * y3_lvar[slice(x_lbin, x_ubin)]
    Y4_NL = coef_4 * y4_lvar[slice(x_lbin, x_ubin)]
    
    Y1_NU = coef_1 * y1_uvar[slice(x_lbin, x_ubin)]
    Y2_NU = coef_2 * y2_uvar[slice(x_lbin, x_ubin)]
    Y3_NU = coef_3 * y3_uvar[slice(x_lbin, x_ubin)]
    Y4_NU = coef_4 * y4_uvar[slice(x_lbin, x_ubin)]
    
    Y_NM = np.mean([Y1_N, Y2_N, Y3_N, Y4_N], axis = 0)
    Y_RMSE = np.sqrt(np.power((Y1_N - Y_NM) / Y_NM, 2) +\
                     np.power((Y2_N - Y_NM) / Y_NM, 2) +\
                     np.power((Y3_N - Y_NM) / Y_NM, 2) +\
                     np.power((Y4_N - Y_NM) / Y_NM, 2)) / 4

    Y_E = np.nan * Y_NM
    Y_O = np.nan * Y_NM

    if len(y1_extr) > 0:
        Y_E = coef_1 * y1_extr[slice(x_lbin, x_ubin)]
        Y_O = Y1_N
    if len(y2_extr) > 0:
        Y_E = coef_2 * y2_extr[slice(x_lbin, x_ubin)]
        Y_O = Y2_N
    if len(y3_extr) > 0:
        Y_E = coef_3 * y3_extr[slice(x_lbin, x_ubin)]
        Y_O = Y3_N
    if len(y4_extr) > 0:
        Y_E = coef_4 * y4_extr[slice(x_lbin, x_ubin)]
        Y_O = Y4_N

    Y_DIFF = 2. * (Y_E - Y_O) / (Y_E + Y_O)
    
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

    x_refr_bin = np.where(x_vals >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))
    

    plt.rc('font', size = 14) 

    # Create the figure
    fig = plt.figure(figsize=(20. , 4.))
    
    fig.suptitle(title)

    # Subplot: Raw Signals
    ax = fig.add_axes([0.07,0.13,0.26,0.7])
    
    # ax2.set_title(title, pad = 15)
    
    ax.plot(X, Y1, color = 'tab:blue')
    ax.plot(X, Y2, color = 'tab:orange')
    ax.plot(X, Y3, color = 'tab:green')
    ax.plot(X, Y4, color = 'tab:red')
    
    ax.legend(['north','east','south','west'])
    
    ax.fill_between(X, Y1_L, Y1_U, color = 'tab:blue', alpha = 0.3)
    ax.fill_between(X, Y2_L, Y2_U, color = 'tab:orange', alpha = 0.3)
    ax.fill_between(X, Y3_L, Y3_U, color = 'tab:green', alpha = 0.3)
    ax.fill_between(X, Y4_L, Y4_U, color = 'tab:red', alpha = 0.3)
    
    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label)
    
    ax.set_ylabel(y_label_1)
    ax.ticklabel_format(axis = 'y', useMathText = True)
    
    ax.grid(which = 'both')

    # Subplot: Normalized Signals
    ax2 = fig.add_axes([0.39,0.13,0.26,0.7])
        
    ax2.plot(X, Y1_N, color = 'tab:blue')
    ax2.plot(X, Y2_N, color = 'tab:orange')
    ax2.plot(X, Y3_N, color = 'tab:green')
    ax2.plot(X, Y4_N, color = 'tab:red')

    ax2.fill_between(X, Y1_NL, Y1_NU, color = 'tab:blue', alpha = 0.3)
    ax2.fill_between(X, Y2_NL, Y2_NU, color = 'tab:orange', alpha = 0.3)
    ax2.fill_between(X, Y3_NL, Y3_NU, color = 'tab:green', alpha = 0.3)
    ax2.fill_between(X, Y4_NL, Y4_NU, color = 'tab:red', alpha = 0.3)

    ax2.set_xticks(x_ticks, labels = x_ticks)
    ax2.set_xlim([x_llim, x_ulim])
    ax2.set_xlabel(x_label)

    ax2.set_ylim([y_llim_nr, y_ulim_nr])
    ax2.set_ylabel('Normalized RC Signals [A.U.]')

    ax2.grid(which = 'both')

    ax2.axvspan(x_vals[x_refr_bin - refr_hbin], x_vals[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')

    n_llim = np.round(x_refr - refr_hwin * 1e-3, decimals = 2)
    n_ulim = np.round(x_refr + refr_hwin * 1e-3, decimals = 2)
    
    ax2.text(0.55 * x_ulim, 0.9 * y_ulim_nr, 
             f'norm: {n_llim} - {n_ulim} Km',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 9))

    # Subplot: Normalized Deviations
    ax3 = fig.add_axes([0.72,0.13,0.26,0.7])
    
    ax3.plot(X, (Y1_N - Y_NM) / Y_NM, color = 'tab:blue', label='_nolegend_')
    ax3.plot(X, (Y2_N - Y_NM) / Y_NM, color = 'tab:orange', label='_nolegend_')
    ax3.plot(X, (Y3_N - Y_NM) / Y_NM, color = 'tab:green', label='_nolegend_')
    ax3.plot(X, (Y4_N - Y_NM) / Y_NM, color = 'tab:red', label='_nolegend_')
    ax3.plot(X, Y_RMSE, color = 'yellow', label = 'Sector RMS')
    if use_last == True:
        ax3.plot(X, Y_DIFF, color = 'black', label = 'Atm. Dif.')
    
    ax3.legend()

    ax3.fill_between(X, (Y1_NL - Y_NM) / Y_NM, (Y1_NU - Y_NM) / Y_NM, 
                     color = 'tab:blue', alpha = 0.3)
    ax3.fill_between(X, (Y2_NL - Y_NM) / Y_NM, (Y2_NU - Y_NM) / Y_NM, 
                     color = 'tab:orange', alpha = 0.3)
    ax3.fill_between(X, (Y3_NL - Y_NM) / Y_NM, (Y3_NU - Y_NM) / Y_NM, 
                     color = 'tab:green', alpha = 0.3)
    ax3.fill_between(X, (Y4_NL - Y_NM) / Y_NM, (Y4_NU - Y_NM) / Y_NM, 
                     color = 'tab:red', alpha = 0.3)

    ax3.plot(X, np.zeros(X.shape), '--', color = 'black', zorder = 10)

    ax3.plot(X, -0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 10, alpha = 0.7)
    ax3.plot(X, 0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 10, alpha = 0.7)

    ax3.set_xticks(x_ticks, labels = x_ticks)
    ax3.set_xlim([x_llim, x_ulim])
    ax3.set_xlabel(x_label)
    
    y_ticks = np.round(np.arange(-0.20, 0.20 + 0.05, 0.05), decimals = 2)
    ax3.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax3.set_ylim([y_ticks[0], y_ticks[-1]])
    ax3.set_ylabel('Relative Sector Deviation')

    ax3.grid(which = 'both')

    fpath = os.path.join(dir_out, 'plots', fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    plt.rcParams.update(plt.rcParamsDefault)
    
    return(fpath)

def telecover_rin(dir_out, fname, title, dpi_val, 
                  x_refr, refr_hwin, x_vals, use_nonrc, 
                  y1_raw, y2_raw, 
                  y1_lraw, y2_lraw, 
                  y1_uraw, y2_uraw, 
                  y1_vals, y2_vals, 
                  y1_lvar, y2_lvar, y1_uvar, y2_uvar,
                  coef_1, coef_2, ranges,
                  x_lbin, x_ubin, x_llim, x_ulim, 
                  y_llim, y_ulim, y_llim_nr, y_ulim_nr, 
                  x_label, x_tick):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]
    
    if use_nonrc == True:
        R = ranges[slice(x_lbin, x_ubin)]
    
        Y1 = y1_raw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y2 = y2_raw[slice(x_lbin, x_ubin)] / np.power(R, 2)
    
        Y1_L = Y1 * y1_lraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y2_L = Y2 * y2_lraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        
        Y1_U = Y1 * y1_uraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        Y2_U = Y2 * y2_uraw[slice(x_lbin, x_ubin)] / np.power(R, 2)
        
        y_label_1 = 'Non RC Signals - Raw Units'
    
    else:
        Y1 = y1_raw[slice(x_lbin, x_ubin)] 
        Y2 = y2_raw[slice(x_lbin, x_ubin)] 
    
        Y1_L = Y1 * y1_lraw[slice(x_lbin, x_ubin)]
        Y2_L = Y2 * y2_lraw[slice(x_lbin, x_ubin)]
        
        Y1_U = Y1 * y1_uraw[slice(x_lbin, x_ubin)]
        Y2_U = Y2 * y2_uraw[slice(x_lbin, x_ubin)]
        
        y_label_1 = 'RC Signals [A.U.]'

    
    Y1_N = coef_1 * y1_vals[slice(x_lbin, x_ubin)]
    Y2_N = coef_2 * y2_vals[slice(x_lbin, x_ubin)]

    Y1_NL = coef_1 * y1_lvar[slice(x_lbin, x_ubin)]
    Y2_NL = coef_2 * y2_lvar[slice(x_lbin, x_ubin)]
    
    Y1_NU = coef_1 * y1_uvar[slice(x_lbin, x_ubin)]
    Y2_NU = coef_2 * y2_uvar[slice(x_lbin, x_ubin)]
    
    Y_NM = np.mean([Y1_N, Y2_N], axis = 0) 

    x_ticks = np.arange(x_tick * np.ceil(x_llim / x_tick), 
                        x_tick * (np.floor(x_ulim / x_tick) + 1.), 
                        x_tick)
    
    x_ticks = np.round(x_ticks, decimals = 2)

        
    if np.abs(x_llim - x_ticks[0]) < x_tick * 0.25:
        x_ticks[0] = x_llim
    else:
        x_ticks = np.hstack((x_llim, x_ticks))

    if np.abs(x_ulim - x_ticks[-1]) < x_tick * 0.25:
        x_ticks[-1] = x_ulim
    else:
        x_ticks = np.hstack((x_ticks, x_ulim))
    
    x_refr_bin = np.where(x_vals >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))
    
    plt.rc('font', size = 14) 

    # Create the figure
    fig = plt.figure(figsize=(20. , 4.))
    fig.suptitle(title)

    # Subplot: Raw Signals
    ax = fig.add_axes([0.07,0.13,0.26,0.7])
    
    # ax2.set_title(title, pad = 15)
    
    ax.plot(X, Y1, color = 'tab:purple')
    ax.plot(X, Y2, color = 'tab:cyan')
    
    ax.legend(['inner','outer'])
    
    ax.fill_between(X, Y1_L, Y1_U, color = 'tab:purple', alpha = 0.3)
    ax.fill_between(X, Y2_L, Y2_U, color = 'tab:cyan', alpha = 0.3)
    
    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label)
    
    ax.set_ylabel(y_label_1)
    
    ax.grid(which = 'both')

    ax2 = fig.add_axes([0.39,0.13,0.26,0.7])
        
    ax2.plot(X, Y1_N, color = 'tab:purple')
    ax2.plot(X, Y2_N, color = 'tab:cyan')

    ax2.legend(['outer','inner'])

    ax2.fill_between(X, Y1_NL, Y1_NU, color = 'tab:purple', alpha = 0.3)
    ax2.fill_between(X, Y2_NL, Y2_NU, color = 'tab:cyan', alpha = 0.3)

    ax2.set_xticks(x_ticks, labels = x_ticks)
    ax2.set_xlim([x_llim, x_ulim])
    ax2.set_xlabel(x_label)

    ax2.set_ylim([y_llim_nr, y_ulim_nr])
    ax2.set_ylabel('Normalized RC Signals')

    ax2.grid(which = 'both')
    
    ax2.axvspan(x_vals[x_refr_bin - refr_hbin], x_vals[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')

    n_llim = np.round(x_refr - refr_hwin * 1e-3, decimals = 2)
    n_ulim = np.round(x_refr + refr_hwin * 1e-3, decimals = 2)
    
    ax2.text(0.55 * x_ulim, 0.9 * y_ulim_nr, 
             f'norm: {n_llim} - {n_ulim} Km',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 9))
    
    # Subplot: Normalized Deviations
    ax3 = fig.add_axes([0.72,0.13,0.26,0.7])
    
    ax3.plot(X, (Y1_N - Y_NM) / Y_NM, color = 'tab:purple')
    ax3.plot(X, (Y2_N - Y_NM) / Y_NM, color = 'tab:cyan')
        
    ax3.fill_between(X, (Y1_NL - Y_NM) / Y_NM, (Y1_NU - Y_NM) / Y_NM, 
                     color = 'tab:purple', alpha = 0.3)
    ax3.fill_between(X, (Y2_NL - Y_NM) / Y_NM, (Y2_NU - Y_NM) / Y_NM, 
                     color = 'tab:cyan', alpha = 0.3)
   
    ax3.plot(X, np.zeros(X.shape), '--', color = 'black', zorder = 8)
    
    ax3.plot(X, -0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 8, alpha = 0.7)
    ax3.plot(X, 0.05 * np.ones(X.shape), '--', 
             color = 'black', zorder = 8, alpha = 0.7)
    
    ax3.set_xticks(x_ticks, labels = x_ticks)
    ax3.set_xlim([x_llim, x_ulim])
    ax3.set_xlabel(x_label)
    
    y_ticks = np.round(np.arange(-0.20, 0.20 + 0.05, 0.05), decimals = 2)
    ax3.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax3.set_ylim([y_ticks[0], y_ticks[-1]])
    ax3.set_ylabel('Relative Sector Deviation')
    
    ax3.grid(which = 'both')
    
    fpath = os.path.join(dir_out, 'plots', fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    plt.rcParams.update(plt.rcParamsDefault)
    
    return(fpath)


def intercomparison(dir_out, fname, title, dpi_val, use_lin, x_refr, refr_hwin,
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
                f'norm: {n_llim} - {n_ulim} Km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
    else:
        ax.text(0.55 * x_ulim, 0.9 * y_ulim, 
                f'norm: {n_llim} - {n_ulim} Km',
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
    
    return(fpath)

def polarization_calibration(dir_out, fname, title, dpi_val, 
                             x_cal, cal_hwin, x_vdr, vdr_hwin,
                             x_vals_cal, x_vals_vdr,
                             y1_vals, y2_vals, y3_vals, y4_vals, y5_vals, y6_vals,
                             eta, eta_f_s, eta_s,
                             delta_m, delta_c, delta, epsilon,
                             eta_err, eta_f_s_err, eta_s_err,
                             delta_c_err, delta_err, epsilon_err,
                             x_lbin_cal, x_ubin_cal, 
                             x_llim_cal, x_ulim_cal, 
                             y_llim_cal, y_ulim_cal, 
                             x_lbin_vdr, x_ubin_vdr, 
                             x_llim_vdr, x_ulim_vdr, 
                             y_llim_vdr, y_ulim_vdr, 
                             K, G_R, G_T, H_R, H_T,
                             y_label_cal, x_label_cal, x_tick_cal,
                             y_label_vdr, x_label_vdr, x_tick_vdr):
        
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

    cal_bin = np.where(x_vals_cal >= x_cal)[0][0]
    cal_hwin_bins = int(1E-3 * cal_hwin / (x_vals_cal[1] - x_vals_cal[0]))

    ax.axvspan(x_vals_cal[cal_bin - cal_hwin_bins], x_vals_cal[cal_bin + cal_hwin_bins + 1],
               alpha = 0.2, facecolor = 'tab:grey')

    c_llim = np.round(x_vals_cal[cal_bin] - cal_hwin * 1e-3, decimals = 2)
    c_ulim = np.round(x_vals_cal[cal_bin] + cal_hwin * 1e-3, decimals = 2)
        
    ax.text(0.45 * x_ulim_cal, 0.95 * y_ulim_cal, 
            f'cal.: {c_llim} - {c_ulim} Km',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))  
    # ax.text(0.55 * x_ulim_cal, 0.87 * y_ulim_cal, 
    #         f'ε: {round_it(epsilon,2)}'+'${}^o \pm$ '+f'{round_it(epsilon_err,2)}'+'${}^o$,',
    #         bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.45 * x_ulim_cal, 0.87 * y_ulim_cal, 
            f'ε: {round_it(epsilon,2)}'+'${}^o$'+ f', K: {round_it(K, 3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.45 * x_ulim_cal, 0.79 * y_ulim_cal, 
            r'$η^{\star}_{f}$'+f': {round_it(eta_f_s, 2)}' + ' $\pm$ ' + f'{round_it(eta_f_s_err, 3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.45 * x_ulim_cal, 0.71 * y_ulim_cal, 
            r'$η^{\star}$'+f': {round_it(eta_s, 3)}' + ' $\pm$ ' + f'{round_it(eta_s_err, 3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.45 * x_ulim_cal, 0.63 * y_ulim_cal, 
            f'η: {round_it(eta, 3)}' + ' $\pm$ ' + f'{round_it(eta_err, 3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 

    ax2 = fig.add_axes([0.55,0.13,0.40,0.65])

    ax2.plot(XB, Y4, color = 'tab:blue', label = 'measured')
    ax2.plot(XB, Y5, color = 'tab:orange', label = 'corrected')
    ax2.plot(XB, Y6, color = 'tab:green', label = 'molecular')

    ax2.fill_between(XB, Y4 - Y4E, Y4 + Y4E, color = 'tab:blue', alpha = 0.3)
    ax2.fill_between(XB, Y5 - Y5E, Y5 + Y5E, color = 'tab:orange', alpha = 0.3)
    
    x_ticks_vdr = np.arange(x_tick_vdr * np.ceil(x_llim_vdr / x_tick_vdr), 
                            x_tick_vdr * (np.floor(x_ulim_vdr / x_tick_vdr) + 1.), 
                            x_tick_vdr)
        
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

    vdr_bin = np.where(x_vals_vdr >= x_vdr)[0][0]
    vdr_hwin_bins = int(1E-3 * vdr_hwin / (x_vals_vdr[1] - x_vals_vdr[0]))

    ax2.axvspan(x_vals_vdr[vdr_bin - vdr_hwin_bins], x_vals_vdr[vdr_bin + vdr_hwin_bins + 1],
                alpha = 0.2, facecolor = 'tab:grey')
    
    m_llim = np.round(x_vals_vdr[vdr_bin] - vdr_hwin * 1e-3, decimals = 2)
    m_ulim = np.round(x_vals_vdr[vdr_bin] + vdr_hwin * 1e-3, decimals = 2)
        
    ax2.text(0.40 * x_ulim_vdr, 0.93 * y_ulim_vdr, 
             f'mol.: {m_llim} - {m_ulim}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))  
    ax2.text(0.40 * x_ulim_vdr, 0.83 * y_ulim_vdr, 
            r'$δ^{\star}$'+f': {round_it(delta_c,3)}' + ' $\pm$ ' + f'{round_it(delta_c_err,2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.40 * x_ulim_vdr, 0.73 * y_ulim_vdr, 
            r'$δ_{c}$'+f': {round_it(delta,3)}' + ' $\pm$ ' + f'{round_it(delta_err,2)}', 
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.40 * x_ulim_vdr, 0.63 * y_ulim_vdr, 
            r'$δ_m$: '+f'{round_it(delta_m,3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.40 * x_ulim_vdr, 0.53 * y_ulim_vdr, 
            r'$G_R$: '+f'{round_it(G_R,4)}, $G_T$: '+f'{round_it(G_T,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.40 * x_ulim_vdr, 0.43 * y_ulim_vdr, 
            r'$H_R$: '+f'{round_it(H_R,4)}, $H_T$: '+f'{round_it(H_T,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 

    fpath = os.path.join(dir_out, 'plots', fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
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