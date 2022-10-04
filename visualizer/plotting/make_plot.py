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
import os, glob
from matplotlib import rcParams

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
    ax = fig.add_axes([0.09,0.13,1.,0.7])
    
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
        
    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)

def rayleigh(dir_out, fname, title, dpi_val, use_lin, x_refr, refr_hwin,
             x_vals, y1_vals, y2_vals, y1_errs, coef,
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

    ax.set_ylim([y_llim, y_ulim])
    ax.set_ylabel(y_label)
    if use_lin == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    ax.legend(loc = 'lower left')

    x_refr_bin = np.where(X >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))

    ax.axvspan(X[x_refr_bin - refr_hbin], X[x_refr_bin + refr_hbin + 1],
               alpha = 0.2, facecolor = 'tab:grey')

    if use_lin == False:
        ax.text(0.55 * x_ulim, 0.60 * y_ulim, f'normalize: {x_refr} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
        ax.text(0.55 * x_ulim, 0.30 * y_ulim, f'window: {2. * refr_hwin} m',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))  
        
    else:
        ax.text(0.55 * x_ulim, 0.9 * y_ulim, f'normalize: {x_refr} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
        ax.text(0.55 * x_ulim, 0.75 * y_ulim, f'window: {2. * refr_hwin} m',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))   
    
    ax2 = fig.add_axes([0.65,0.13,0.30,0.7])

    ax2.fill_between(X, (Y1 - Y1E - Y2) / Y2, 
                     (Y1 + Y1E - Y2) / Y2, color = 'tab:blue', 
                     alpha = 0.3, label = 'sem')
    ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue',label = 'mean')
    
    ax2.axhline(c = 'k')

    x_ticks = np.arange(2. * x_tick * np.floor(x_llim / x_tick), 
                        2. * x_tick * (np.ceil(x_ulim / x_tick) + 1.), 
                        2. * x_tick)
    
    ax2.set_xticks(x_ticks, labels = x_ticks)
    ax2.set_xlim([x_llim, x_ulim])
    ax2.set_xlabel(x_label)

    y_ticks = np.round(np.arange(-0.40, 0.40 + 0.10, 0.10), decimals = 2)
    ax2.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax2.set_ylim([y_ticks[0], y_ticks[-1]])
    ax2.set_ylabel('Relative Attenuated Bck. Diff. ')
    
    ax2.grid(which = 'both')
    ax2.legend(loc = 'lower left')
    
    ax2.axvspan(X[x_refr_bin - refr_hbin], X[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')

    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)

def telecover_sec(dir_out, fname, title, dpi_val, x_refr, refr_hwin, x_vals, 
                  y1_vals, y2_vals, y3_vals, y4_vals,
                  y1_extr, y2_extr, y3_extr, y4_extr,
                  y1_lvar, y2_lvar, y3_lvar, y4_lvar, 
                  y1_uvar, y2_uvar, y3_uvar, y4_uvar, 
                  coef_1, coef_2, coef_3, coef_4,
                  x_lbin, x_ubin, x_llim, x_ulim, 
                  y_llim, y_ulim, y_llim_nr, y_ulim_nr, 
                  x_label, x_tick):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]

    Y1 = y1_vals[slice(x_lbin, x_ubin)]
    Y2 = y2_vals[slice(x_lbin, x_ubin)]
    Y3 = y3_vals[slice(x_lbin, x_ubin)]
    Y4 = y4_vals[slice(x_lbin, x_ubin)]

    Y1_L = y1_lvar[slice(x_lbin, x_ubin)]
    Y2_L = y2_lvar[slice(x_lbin, x_ubin)]
    Y3_L = y3_lvar[slice(x_lbin, x_ubin)]
    Y4_L = y4_lvar[slice(x_lbin, x_ubin)]
    
    Y1_U = y1_uvar[slice(x_lbin, x_ubin)]
    Y2_U = y2_uvar[slice(x_lbin, x_ubin)]
    Y3_U = y3_uvar[slice(x_lbin, x_ubin)]
    Y4_U = y4_uvar[slice(x_lbin, x_ubin)]
    
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
        Y_E = y1_extr[slice(x_lbin, x_ubin)]
        Y_O = Y1
    if len(y2_extr) > 0:
        Y_E = y2_extr[slice(x_lbin, x_ubin)]
        Y_O = Y2
    if len(y3_extr) > 0:
        Y_E = y3_extr[slice(x_lbin, x_ubin)]
        Y_O = Y3
    if len(y4_extr) > 0:
        Y_E = y4_extr[slice(x_lbin, x_ubin)]
        Y_O = Y4

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

    x_refr_bin = np.where(X >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))
    

    plt.rc('font', size = 14) 

    # Create the figure
    fig = plt.figure(figsize=(20. , 4.))
    
    fig.suptitle(title)

    # Subplot: Raw Signals
    ax = fig.add_axes([0.07,0.13,0.23,0.7])
    
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
    
    ax.set_ylim([y_llim, y_ulim])
    ax.set_ylabel('RC Signals [A.U.]')
    ax.ticklabel_format(axis = 'y', useMathText = True)
    
    ax.grid(which = 'both')

    # Subplot: Normalized Signals
    ax2 = fig.add_axes([0.37,0.13,0.23,0.7])
        
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
    ax2.set_ylabel('Normalized RC Signals')

    ax2.grid(which = 'both')

    ax2.axvspan(X[x_refr_bin - refr_hbin], X[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')

    ax2.text(0.55 * x_ulim, 0.9 * y_ulim_nr, f'normalize: {x_refr} km',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 9))
    ax2.text(0.55 * x_ulim, 0.75 * y_ulim_nr, f'window: {2. * refr_hwin} m',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 9))   
    
    # Subplot: Normalized Deviations
    ax3 = fig.add_axes([0.67,0.13,0.23,0.7])
    
    ax3.plot(X, (Y1_N - Y_NM) / Y_NM, color = 'tab:blue', label='_nolegend_')
    ax3.plot(X, (Y2_N - Y_NM) / Y_NM, color = 'tab:orange', label='_nolegend_')
    ax3.plot(X, (Y3_N - Y_NM) / Y_NM, color = 'tab:green', label='_nolegend_')
    ax3.plot(X, (Y4_N - Y_NM) / Y_NM, color = 'tab:red', label='_nolegend_')
    ax3.plot(X, Y_RMSE, color = 'yellow', label = 'Sector RMS')
    ax3.plot(X, Y_DIFF, color = 'black', label = 'Atm. Dif.')
    
    ax3.legend(['all_dev', 'atm_dif'])

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

    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    plt.rcParams.update(plt.rcParamsDefault)
    
    return(fpath)

def telecover_rin(dir_out, fname, title, dpi_val, x_refr, refr_hwin, x_vals, 
                  y1_vals, y2_vals, 
                  y1_lvar, y2_lvar, y1_uvar, y2_uvar,
                  coef_1, coef_2, 
                  x_lbin, x_ubin, x_llim, x_ulim, 
                  y_llim, y_ulim, y_llim_nr, y_ulim_nr, 
                  x_label, x_tick):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]

    Y1 = y1_vals[slice(x_lbin, x_ubin)]
    Y2 = y2_vals[slice(x_lbin, x_ubin)]

    Y1_L = Y1 * y1_lvar[slice(x_lbin, x_ubin)]
    Y2_L = Y2 * y2_lvar[slice(x_lbin, x_ubin)]
    
    Y1_U = Y1 * y1_uvar[slice(x_lbin, x_ubin)]
    Y2_U = Y2 * y2_uvar[slice(x_lbin, x_ubin)]
    
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
    
    x_refr_bin = np.where(X >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))
    
    plt.rc('font', size = 14) 

    # Create the figure
    fig = plt.figure(figsize=(20. , 4.))
    fig.suptitle(title)

    # Subplot: Raw Signals
    ax = fig.add_axes([0.07,0.13,0.23,0.7])
    
    # ax2.set_title(title, pad = 15)
    
    ax.plot(X, Y1, color = 'tab:purple')
    ax.plot(X, Y2, color = 'tab:cyan')
    
    ax.legend(['inner','outer'])
    
    ax.fill_between(X, Y1_L, Y1_U, color = 'tab:purple', alpha = 0.3)
    ax.fill_between(X, Y2_L, Y2_U, color = 'tab:cyan', alpha = 0.3)
    
    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label)
    
    ax.set_ylim([y_llim, y_ulim])
    ax.set_ylabel('RC Signals [A.U.]')
    
    ax.grid(which = 'both')

    ax2 = fig.add_axes([0.37,0.13,0.23,0.7])
        
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
    
    ax2.axvspan(X[x_refr_bin - refr_hbin], X[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')

    ax2.text(0.55 * x_ulim, 0.9 * y_ulim_nr, f'normalize: {x_refr} km',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 5))
    ax2.text(0.55 * x_ulim, 0.75 * y_ulim_nr, f'window: {2. * refr_hwin} m',
             bbox=dict(facecolor='tab:cyan', alpha=0.1, zorder = 5))   
    
    # Subplot: Normalized Deviations
    ax3 = fig.add_axes([0.67,0.13,0.23,0.7])
    
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
    
    fpath = os.path.join(dir_out, fname)
    
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
    
    x_refr_bin = np.where(X >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))

    ax.axvspan(X[x_refr_bin - refr_hbin], X[x_refr_bin + refr_hbin + 1],
               alpha = 0.2, facecolor = 'tab:grey')

    if use_lin == False:
        ax.text(0.55 * x_ulim, 0.60 * y_ulim, f'normalize: {x_refr} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
        ax.text(0.55 * x_ulim, 0.30 * y_ulim, f'window: {2. * refr_hwin} m',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))  
        
    else:
        ax.text(0.55 * x_ulim, 0.9 * y_ulim, f'normalize: {x_refr} km',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
        ax.text(0.55 * x_ulim, 0.75 * y_ulim, f'window: {2. * refr_hwin} m',
                bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))   
    
        
    ax2 = fig.add_axes([0.65,0.13,0.30,0.7])

    # ax2.set_title('', pad = 15)

    frac_bias = 2. * (Y1 - Y2) / (Y1 + Y2)
    ax2.plot(X, frac_bias, color = 'tab:blue')

    x_ticks = np.arange(2. * x_tick * np.floor(x_llim / x_tick), 
                        2. * x_tick * (np.ceil(x_ulim / x_tick) + 1.), 
                        2. * x_tick)
    
    ax2.set_xticks(x_ticks, labels = x_ticks)
    ax2.set_xlim([x_llim, x_ulim])
    ax2.set_xlabel(x_label)

    y_ticks = np.round(np.arange(-0.40, 0.40 + 0.10, 0.10), decimals = 2)
    ax2.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax2.set_ylim([y_ticks[0], y_ticks[-1]])
    ax2.set_ylabel(f'Fractional Bias ({lidars[0]} - {lidars[1]})')
    
    ax2.grid(which = 'both')
    
    ax2.axhline(c = 'k')

    ax2.axvspan(X[x_refr_bin - refr_hbin], X[x_refr_bin + refr_hbin + 1],
                alpha = 0.2, facecolor = 'tab:grey')
    
        
    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()

    plt.close()
    
    return(fpath)

def polarization_calibration(dir_out, fname, title, dpi_val, 
                             y_cal, cal_hwin, y_vdr, vdr_hwin,
                             y_vals_cal, y_vals_vdr,
                             x1_vals, x2_vals, x3_vals, x4_vals, x5_vals, x6_vals,
                             eta, eta_f_s, eta_s,
                             delta_m, delta_c, delta, epsilon,
                             y_lbin_cal, y_ubin_cal, 
                             y_llim_cal, y_ulim_cal, 
                             x_llim_cal, x_ulim_cal, 
                             y_lbin_vdr, y_ubin_vdr, 
                             y_llim_vdr, y_ulim_vdr, 
                             x_llim_vdr, x_ulim_vdr, 
                             x_label_cal, y_label_cal, y_tick_cal,
                             x_label_vdr, y_label_vdr, y_tick_vdr):
        
    # Create the variables to be plotted X, Y
    YA = y_vals_cal[slice(y_lbin_cal, y_ubin_cal)]
    YB = y_vals_vdr[slice(y_lbin_vdr, y_ubin_vdr)]
    
    X1 = x1_vals[slice(y_lbin_cal, y_ubin_cal)]
    X2 = x2_vals[slice(y_lbin_cal, y_ubin_cal)]
    X3 = x3_vals[slice(y_lbin_cal, y_ubin_cal)]
    X4 = x4_vals[slice(y_lbin_vdr, y_ubin_vdr)]
    X5 = x5_vals[slice(y_lbin_vdr, y_ubin_vdr)]
    X6 = x6_vals[slice(y_lbin_vdr, y_ubin_vdr)]

    X1E = np.nan * X1
    X2E = np.nan * X2
    X3E = np.nan * X3
    X4E = np.nan * X4
    X5E = np.nan * X5
    
    # Create the figure
    fig = plt.figure(figsize=(8. , 6.))
    fig.suptitle(title)

    ax = fig.add_axes([0.07,0.13,0.40,0.70])
        
    ax.plot(X1, YA, color = 'tab:purple', label = 'η')
    ax.plot(X2, YA, color = 'tab:red', label = '$η_{+45}$')
    ax.plot(X3, YA, color = 'tab:cyan', label = '$η_{-45}$')

    ax.fill_betweenx(YA, X1 - X1E, X1 + X1E, color = 'tab:purple', alpha = 0.3)
    ax.fill_betweenx(YA, X2 - X2E, X2 + X2E, color = 'tab:red', alpha = 0.3)
    ax.fill_betweenx(YA, X3 - X3E, X3 + X3E, color = 'tab:cyan', alpha = 0.3)
    
    y_ticks_cal = np.arange(y_tick_cal * np.ceil(y_llim_cal / y_tick_cal), 
                            y_tick_cal * (np.floor(y_ulim_cal / y_tick_cal) + 1.), 
                            y_tick_cal)
        
    if np.abs(y_llim_cal - y_ticks_cal[0]) < y_tick_cal * 0.25:
        y_ticks_cal[0] = y_llim_cal
    else:
        y_ticks_cal = np.hstack((y_llim_cal, y_ticks_cal))

    if np.abs(y_ulim_cal - y_ticks_cal[-1]) < y_tick_cal * 0.25:
        y_ticks_cal[-1] = y_ulim_cal
    else:
        y_ticks_cal = np.hstack((y_ticks_cal, y_ulim_cal))

    y_ticks_cal = np.round(y_ticks_cal, decimals = 2)

    ax.set_yticks(y_ticks_cal, labels = y_ticks_cal)
    ax.set_ylim([y_llim_cal, y_ulim_cal])
    ax.set_ylabel(y_label_cal)

    ax.set_xlim([x_llim_cal, x_ulim_cal])
    ax.set_xlabel(x_label_cal)

    ax.grid(which = 'both')
    ax.legend(loc = 'upper right')

    cal_bin = np.where(YB >= y_cal)[0][0]
    cal_hwin_bins = int(1E-3 * cal_hwin / (y_vals_cal[1] - y_vals_cal[0]))

    ax.axhspan(YA[cal_bin - cal_hwin_bins], YA[cal_bin + cal_hwin_bins + 1],
               alpha = 0.2, facecolor = 'tab:grey')

    
    ax.text(0.65 * x_ulim_cal, 0.70 * y_ulim_cal, 
            f'cal.: {np.round(YA[cal_bin], decimals = 2)} km',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
    ax.text(0.65 * x_ulim_cal, 0.63 * y_ulim_cal, 
            f'win.: {np.round(2. * cal_hwin, decimals = 0)} m',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))   
    ax.text(0.65 * x_ulim_cal, 0.56 * y_ulim_cal, 
            f'ε: {round_it(epsilon,2)}'+'${}^o$',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.65 * x_ulim_cal, 0.49 * y_ulim_cal, 
            r'$η^{\star}_{f}$'+f': {round_it(eta_f_s, 3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.65 * x_ulim_cal, 0.42 * y_ulim_cal, 
            r'$η^{\star}$'+f': {round_it(eta_s, 3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax.text(0.65 * x_ulim_cal, 0.35 * y_ulim_cal, 
            f'η: {round_it(eta, 3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
        
    ax2 = fig.add_axes([0.56,0.13,0.40,0.70])

    ax2.plot(X4, YB, color = 'tab:blue', label = 'measured')
    ax2.plot(X5, YB, color = 'tab:orange', label = 'corrected')
    ax2.plot(X6, YB, color = 'tab:green', label = 'molecular')

    ax2.fill_betweenx(YB, X4 - X4E, X4 + X4E, color = 'tab:blue', alpha = 0.3)
    ax2.fill_betweenx(YB, X5 - X5E, X5 + X5E, color = 'tab:orange', alpha = 0.3)
    
    y_ticks_vdr = np.arange(y_tick_vdr * np.ceil(y_llim_vdr / y_tick_vdr), 
                            y_tick_vdr * (np.floor(y_ulim_vdr / y_tick_vdr) + 1.), 
                            y_tick_vdr)
        
    if np.abs(y_llim_vdr - y_ticks_vdr[0]) < y_tick_vdr * 0.25:
        y_ticks_vdr[0] = y_llim_vdr
    else:
        y_ticks_vdr = np.hstack((y_llim_vdr, y_ticks_vdr))

    if np.abs(y_ulim_vdr - y_ticks_vdr[-1]) < y_tick_vdr * 0.25:
        y_ticks_vdr[-1] = y_ulim_vdr
    else:
        y_ticks_vdr = np.hstack((y_ticks_vdr, y_ulim_vdr))

    y_ticks_vdr = np.round(y_ticks_vdr, decimals = 2)

    ax2.set_yticks(y_ticks_vdr, labels = y_ticks_vdr)
    ax2.set_ylim([y_llim_vdr, y_ulim_vdr])
    # ax2.set_ylabel(y_label_vdr)

    ax2.set_xlim([x_llim_vdr, x_ulim_vdr])
    ax2.set_xlabel(x_label_vdr)

    ax2.grid(which = 'both')
    ax2.legend(loc = 'upper right')

    vdr_bin = np.where(YB >= y_vdr)[0][0]
    vdr_hwin_bins = int(1E-3 * vdr_hwin / (y_vals_vdr[1] - y_vals_vdr[0]))

    ax2.axhspan(YB[vdr_bin - vdr_hwin_bins], YB[vdr_bin + vdr_hwin_bins + 1],
                alpha = 0.2, facecolor = 'tab:grey')
    
    ax2.text(0.65 * x_ulim_vdr, 0.70 * y_ulim_vdr, 
             f'mol.: {np.round(YB[vdr_bin], decimals = 2)} km',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))
    ax2.text(0.65 * x_ulim_vdr, 0.63 * y_ulim_vdr, 
             f'win.: {np.round(2. * vdr_hwin, decimals = 0)} m',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3))   
    ax2.text(0.65 * x_ulim_vdr, 0.56 * y_ulim_vdr, 
            r'$δ^{\star}$'+f': {round_it(delta_c,3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.65 * x_ulim_vdr, 0.49 * y_ulim_vdr, 
            r'$δ_{c}$'+f': {round_it(delta,3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 
    ax2.text(0.65 * x_ulim_vdr, 0.42 * y_ulim_vdr, 
            r'$δ_m$: '+f'{round_it(delta_m,3)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.1, zorder = 3)) 

    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)

def round_it(x, sig):
    
    if not np.isfinite(x) or np.isnan(x):
        x = -999.
        sig = 3
        
    return np.round(x, sig-int(np.floor(np.log10(abs(x))))-1)