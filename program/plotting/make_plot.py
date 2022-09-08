#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:26:22 2022

@author: nick
"""

from plotting import color_lib, make_colormap
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
    
    # Define the colorscale
    rgb = color_lib.volkers_rgb()
    my_cmap = make_colormap.custom_rgb(rgb, name = 'volkers')

    # Create the variables to be plotted X, Y, Z
    [X, Y] = np.meshgrid(t_vals[slice(x_lbin, x_ubin+1)], 
                         y_vals[slice(y_lbin, y_ubin+1)])
    Z = z_vals[slice(y_lbin, y_ubin+1), slice(x_lbin, x_ubin+1)]
    
    # Create the figure
    fig = plt.figure(figsize=(10. , 5.))
    ax = fig.add_axes([0.07,0.13,0.9,0.7])
    
    ax.set_title(title, pad = 15)
    
    qck = ax.pcolormesh(X, Y, Z, vmin = z_llim, vmax = z_ulim, cmap = my_cmap)
    fig.colorbar(qck, label = z_label, extend = 'both')
    
    ax.set_xlabel(t_label)
    ax.set_ylabel(y_label)

    y_ticks = np.arange(2. * y_tick * np.floor(y_llim / y_tick), 
                        2. * y_tick * (np.ceil(y_ulim / y_tick) + 1.), 
                        2. * y_tick)
    
    ax.set_yticks(y_ticks, labels = y_ticks)
    ax.set_ylim([y_llim, y_ulim])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval = int(t_tick)))

    plt.xticks(rotation = 35)

    if len(nodes) > 0:
        ax1 = ax.twiny()
    
        x_ticks = np.arange(x_lbin, x_ubin + x_tick, x_tick).astype(int)
        ax1.set_xticks(x_ticks, labels = x_ticks)
        
    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)

def rayleigh(dir_out, fname, title, dpi_val, use_lin, x_refr, refr_hwin,
             x_vals, y1_vals, y2_vals, y3_vals,
             x_lbin, x_ubin, x_llim, x_ulim, y_llim, y_ulim, 
             x_label, y_label, x_tick):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]
    Y1 = y1_vals[slice(x_lbin, x_ubin)]
    Y2 = y2_vals[slice(x_lbin, x_ubin)]
    Y3 = y3_vals[slice(x_lbin, x_ubin)]
    
    # Create the figure
    fig = plt.figure(figsize=(12. , 4.))
    fig.suptitle(title)

    ax = fig.add_axes([0.07,0.13,0.50,0.7])
    
    # ax.set_title(title, pad = 15)
    
    ax.plot(X, Y1, color = 'tab:blue')
    ax.plot(X, Y2, color = 'tab:red')

    x_ticks = np.arange(x_tick * np.floor(x_llim / x_tick), 
                        x_tick * (np.ceil(x_ulim / x_tick) + 1.), 
                        x_tick)

    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label)

    ax.set_ylim([y_llim, y_ulim])
    ax.set_ylabel(y_label)
    if use_lin == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    ax.legend(['measured','molecular'])

    x_refr_bin = np.where(X >= x_refr)[0][0]
    refr_hbin = int(1E-3 * refr_hwin / (x_vals[1] - x_vals[0]))

    ax.scatter(X[x_refr_bin], 
               np.mean(Y2[x_refr_bin-refr_hbin:x_refr_bin+refr_hbin+1]), 
               marker = '*', s = 300, color = 'black', zorder = 2)
        
    ax2 = fig.add_axes([0.65,0.13,0.30,0.7])

    # ax2.set_title('', pad = 15)

    ax2.plot(X, (Y3 - Y2) / Y2, color = 'tab:green', alpha = 0.5)
    ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue')

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
    ax2.legend(['unsmoothed','smoothed'])
    
    ax2.plot(X, 0. * Y2, color = 'black')
    
    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)

def telecover_sec(dir_out, fname, title, dpi_val, x_refr, refr_hwin, x_vals, 
                  y1_vals, y2_vals, y3_vals, y4_vals,
                  y1_norm, y2_norm, y3_norm, y4_norm,
                  y1_extr, y2_extr, y3_extr, y4_extr,
                  y1_lvar, y2_lvar, y3_lvar, y4_lvar, 
                  y1_uvar, y2_uvar, y3_uvar, y4_uvar, 
                  x_lbin, x_ubin, x_llim, x_ulim, 
                  y_llim, y_ulim, y_llim_nr, y_ulim_nr, 
                  x_label, x_tick):
        
    # Create the variables to be plotted X, Y
    X = x_vals[slice(x_lbin, x_ubin)]

    Y1 = y1_vals[slice(x_lbin, x_ubin)]
    Y2 = y2_vals[slice(x_lbin, x_ubin)]
    Y3 = y3_vals[slice(x_lbin, x_ubin)]
    Y4 = y4_vals[slice(x_lbin, x_ubin)]

    Y1_L = Y1 * y1_lvar[slice(x_lbin, x_ubin)]
    Y2_L = Y2 * y2_lvar[slice(x_lbin, x_ubin)]
    Y3_L = Y3 * y3_lvar[slice(x_lbin, x_ubin)]
    Y4_L = Y4 * y4_lvar[slice(x_lbin, x_ubin)]
    
    Y1_U = Y1 * y1_uvar[slice(x_lbin, x_ubin)]
    Y2_U = Y2 * y2_uvar[slice(x_lbin, x_ubin)]
    Y3_U = Y3 * y3_uvar[slice(x_lbin, x_ubin)]
    Y4_U = Y4 * y4_uvar[slice(x_lbin, x_ubin)]
    
    Y1_N = y1_norm[slice(x_lbin, x_ubin)]
    Y2_N = y2_norm[slice(x_lbin, x_ubin)]
    Y3_N = y3_norm[slice(x_lbin, x_ubin)]
    Y4_N = y4_norm[slice(x_lbin, x_ubin)]

    Y_M = np.mean([Y1_N, Y2_N, Y3_N, Y4_N], axis = 0)
    Y_RMSE = np.sqrt(np.power((Y1_N - Y_M) / Y_M, 2) +\
                     np.power((Y2_N - Y_M) / Y_M, 2) +\
                     np.power((Y3_N - Y_M) / Y_M, 2) +\
                     np.power((Y4_N - Y_M) / Y_M, 2)) / 4

    Y1_NL = Y1_N * y1_lvar[slice(x_lbin, x_ubin)]
    Y2_NL = Y2_N * y2_lvar[slice(x_lbin, x_ubin)]
    Y3_NL = Y3_N * y3_lvar[slice(x_lbin, x_ubin)]
    Y4_NL = Y4_N * y4_lvar[slice(x_lbin, x_ubin)]
    
    Y1_NU = Y1_N * y1_uvar[slice(x_lbin, x_ubin)]
    Y2_NU = Y2_N * y2_uvar[slice(x_lbin, x_ubin)]
    Y3_NU = Y3_N * y3_uvar[slice(x_lbin, x_ubin)]
    Y4_NU = Y4_N * y4_uvar[slice(x_lbin, x_ubin)]

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
    
        
    x_ticks = np.arange(x_tick * np.floor(x_llim / x_tick), 
                        x_tick * (np.ceil(x_ulim / x_tick) + 1.), 
                        x_tick)
    
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
    
    # ax2.set_title(title, pad = 15)
    
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

    ax2.scatter(X[x_refr_bin], 
               np.mean(Y_M[x_refr_bin-refr_hbin:x_refr_bin+refr_hbin+1]), 
               marker = '*', s = 300, color = 'black', zorder = 8)
    
    # Subplot: Normalized Deviations
    ax3 = fig.add_axes([0.67,0.13,0.23,0.7])
    
    ax3.plot(X, (Y1_N - Y_M) / Y_M, color = 'tab:blue', label='_nolegend_')
    ax3.plot(X, (Y2_N - Y_M) / Y_M, color = 'tab:orange', label='_nolegend_')
    ax3.plot(X, (Y3_N - Y_M) / Y_M, color = 'tab:green', label='_nolegend_')
    ax3.plot(X, (Y4_N - Y_M) / Y_M, color = 'tab:red', label='_nolegend_')
    ax3.plot(X, Y_RMSE, color = 'yellow', label = 'Sector RMS')
    ax3.plot(X, Y_DIFF, color = 'black', label = 'Atm. Dif.')
    
    ax3.legend(['all_dev', 'atm_dif'])

    ax3.fill_between(X, (Y1_NL - Y_M) / Y_M, (Y1_NU - Y_M) / Y_M, 
                     color = 'tab:blue', alpha = 0.3)
    ax3.fill_between(X, (Y2_NL - Y_M) / Y_M, (Y2_NU - Y_M) / Y_M, 
                     color = 'tab:orange', alpha = 0.3)
    ax3.fill_between(X, (Y3_NL - Y_M) / Y_M, (Y3_NU - Y_M) / Y_M, 
                     color = 'tab:green', alpha = 0.3)
    ax3.fill_between(X, (Y4_NL - Y_M) / Y_M, (Y4_NU - Y_M) / Y_M, 
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
                  y1_vals, y2_vals, y1_norm, y2_norm, 
                  y1_lvar, y2_lvar, y1_uvar, y2_uvar, 
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
    
    Y1_N = y1_norm[slice(x_lbin, x_ubin)]
    Y2_N = y2_norm[slice(x_lbin, x_ubin)]

    Y_M = np.mean([Y1_N, Y2_N], axis = 0)

    Y1_NL = Y1_N * y1_lvar[slice(x_lbin, x_ubin)]
    Y2_NL = Y2_N * y2_lvar[slice(x_lbin, x_ubin)]
    
    Y1_NU = Y1_N * y1_uvar[slice(x_lbin, x_ubin)]
    Y2_NU = Y2_N * y2_uvar[slice(x_lbin, x_ubin)]        

    x_ticks = np.arange(x_tick * np.floor(x_llim / x_tick), 
                        x_tick * (np.ceil(x_ulim / x_tick) + 1.), 
                        x_tick)
    
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
    
    # ax2.set_title(title, pad = 15)
    
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

    ax2.scatter(X[x_refr_bin], 
                np.mean(Y_M[x_refr_bin-refr_hbin:x_refr_bin+refr_hbin+1]), 
                marker = '*', s = 300, color = 'black', zorder = 4)
    
    # Subplot: Normalized Deviations
    ax3 = fig.add_axes([0.67,0.13,0.23,0.7])
    
    ax3.plot(X, (Y1_N - Y_M) / Y_M, color = 'tab:purple')
    ax3.plot(X, (Y2_N - Y_M) / Y_M, color = 'tab:cyan')
        
    ax3.fill_between(X, (Y1_NL - Y_M) / Y_M, (Y1_NU - Y_M) / Y_M, 
                     color = 'tab:purple', alpha = 0.3)
    ax3.fill_between(X, (Y2_NL - Y_M) / Y_M, (Y2_NU - Y_M) / Y_M, 
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