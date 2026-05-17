#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:26:22 2022

@author: nick
"""

import os
import numpy as np
import pandas as pd
from PIL import Image
from xarray import DataArray
from PIL import PngImagePlugin
from utils.toolbox import round_it
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from visualizer.smoothing import sliding_average_2D
from visualizer.smoothing import sliding_average_2D_fast
from visualizer.smoothing import sliding_average_1D_fast as smooth_1D

def telecover_sec(dir_out, fname, title, dpi_val, color_reduction, 
                  norm_region,
                  use_nonrc, x_vals, 
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
    
    if norm_region[0] > 20. or norm_region[1] < x_llim:
        raise Exception(f"The normalization_region ({norm_region}) for the telecover fit is out of the x axis limits limits ([{x_llim}, 20]). Please revise the settings_file.ini ")
    
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
        
    if ax.get_legend_handles_labels() != ([], []):
        ax.legend()

    ax.plot([X[0], X[-1]], [0., 0.], color = 'black')
    
    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label, loc = "right")
    
    raw_ulim = np.nanmax([Y1[x_lbin:x_ubin], Y2[x_lbin:x_ubin],
                          Y3[x_lbin:x_ubin], Y4[x_lbin:x_ubin]])
    
    if np.isfinite(raw_ulim):    
        ax.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])
    else:
        ax.set_ylim([0, 1])

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

    if np.isfinite(raw_ulim):    
        ax2.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])
    else:
        ax2.set_ylim([0, 1])
        
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
    ax3.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')

    ax3.set_xticks(x_ticks, labels = x_ticks)
    ax3.set_xlim([x_llim, x_ulim])
    ax3.set_xlabel(x_label, loc = "right")

    ax3.set_ylim([y_llim_nr, y_ulim_nr])
    ax3.set_ylabel('Normalized RC Signals [A.U.]')
    ax3.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.))

    ax3.grid(which = 'both')

    n_llim = np.round(norm_region[0], decimals = 2)
    n_ulim = np.round(norm_region[1], decimals = 2)
    
    c_norm = 'tab:green'
        
    ax3.text(0.30 * x_ulim, 0.90 * y_ulim_nr, 
             f'norm. region: {n_llim} - {n_ulim} km',
             bbox=dict(facecolor=c_norm, alpha=0.1, zorder = 9))
    
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

    ax4.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')

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
    ax5.plot(X, Y_RMSE, color = 'tab:cyan', label = 'RMS Sector Diff.')
    
    if use_last == True and extra_exists:
        ax5.plot(X, Y_DIFF, color = 'tab:purple', label = f'{extra_label_short} - {extra_label_short[0]}')
    
    if ax5.get_legend_handles_labels() != ([], []):
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

    fpath = os.path.join(dir_out, fname)
   
    fig.savefig(fpath, dpi = dpi_val)

    fig.clf()
    
    plt.close()
    
    plt.rcParams.update(plt.rcParamsDefault)

    perform_color_reduction(color_reduction, fpath)

    return(fpath)


def telecover_rin(dir_out, fname, title, dpi_val, color_reduction,
                  norm_region,
                  use_nonrc, x_vals, 
                  y1_raw, y2_raw, 
                  y1_vals, y2_vals,
                  y1_extr, y2_extr,
                  y1_extr_raw, y2_extr_raw,
                  y1_lvar, y2_lvar, 
                  y1_uvar, y2_uvar,
                  coef_1, coef_2, 
                  coef_extra_1, coef_extra_2, extra_sec, ranges,
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

    Y_E_N = np.nan * Y_NM
    
    Y_E = np.nan * Y_NM
    Y_O = np.nan * Y_NM

    if extra_sec['O']:
        Y_E = y1_extr_raw
        Y_E_N = coef_extra_1 * y1_extr
        Y_O = Y1_N
        extra_label = f'outer{iters+1}'
        extra_label_short = f'O{iters+1}'
        extra_exists = True
    elif extra_sec['I']:
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
    
    if norm_region[0] > 20. or norm_region[1] < x_llim:
        raise Exception(f"The normalization_region ({norm_region}) for the telecover fit is out of the x axis limits limits ([{x_llim}, 20]). Please revise the settings_file.ini ")
    
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
        
    
    if ax.get_legend_handles_labels() != ([], []):
        ax.legend()
    
    ax.plot([X[0], X[-1]], [0., 0.], color = 'black')

    ax.set_xticks(x_ticks, labels = x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label, loc = "right")
    
    raw_ulim = np.nanmax([Y1[x_lbin:x_ubin], Y2[x_lbin:x_ubin]])
    
    if np.isfinite(raw_ulim):    
        ax.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])
    else:
        ax.set_ylim([0, 1])
        
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

    ax2.xaxis.set_minor_locator(MultipleLocator(2.))

    if np.isfinite(raw_ulim):    
        ax2.set_ylim([-0.1*raw_ulim, 1.1 * raw_ulim])
    else:
        ax2.set_ylim([0, 1])

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
    
    c_norm = 'tab:green'
        
    ax3.text(0.30 * x_ulim, 0.90 * y_ulim_nr, 
             f'norm. region: {n_llim} - {n_ulim} km',
             bbox=dict(facecolor=c_norm, alpha=0.1, zorder = 9))

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
    
    if ax5.get_legend_handles_labels() != ([], []):
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

    fpath = os.path.join(dir_out, fname)
   
    fig.savefig(fpath, dpi = dpi_val)

    fig.clf()
    
    plt.close()
    
    plt.rcParams.update(plt.rcParamsDefault)

    perform_color_reduction(color_reduction, fpath)
            
    return(fpath)

def polarization_calibration(dir_out, fname, title, dpi_val, color_reduction,
                             cal_region, vdr_region,
                             x_vals_cal, x_vals_vdr,
                             y1_vals, y2_vals, y3_vals, y4_vals, y5_vals, y6_vals,
                             eta, eta_f_s, eta_s,
                             delta_m, delta_c_def, delta_c, epsilon,
                             eta_err, eta_f_s_err, eta_s_err,
                             delta_c_def_err, delta_c_err, 
                             delta_l, delta_l_err, epsilon_err, sr_lim, err_p,
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
    fig = plt.figure(figsize=(16. , 3.2))
    fig.suptitle(title)

    ax = fig.add_axes([0.045,0.14,0.44,0.65])
        
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
    
    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc = 'upper right')

    ax.axvspan(cal_region[0], cal_region[1],
               alpha = 0.2, facecolor = 'tab:grey')

    c_llim = np.round(cal_region[0], decimals = 2)
    c_ulim = np.round(cal_region[1], decimals = 2)
    
    c_cal = 'tab:green'
        
    ax.text(0.05 * x_ulim_cal, 0.94 * y_ulim_cal, 
            f'cal. region: {c_llim} - {c_ulim} km',
            bbox = dict(facecolor = c_cal, alpha = 0.22, zorder = 3))  
    ax.text(0.05 * x_ulim_cal, 0.85 * y_ulim_cal, 
            f'ε: {round_it(epsilon,2)}' +'${}^o$'+ ' $\pm$ ' + f'{round_it(epsilon_err, 2)}' +'${}^o$' + f', K: {round_it(K, 4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax.text(0.05 * x_ulim_cal, 0.76 * y_ulim_cal, 
            r'$η^{\star}_{f}$'+f': {round_it(eta_f_s, 3)}' + ' $\pm$ ' + f'{round_it(eta_f_s_err, 2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax.text(0.05 * x_ulim_cal, 0.67 * y_ulim_cal, 
            r'$η^{\star}$'+f': {round_it(eta_s, 3)}' + ' $\pm$ ' + f'{round_it(eta_s_err, 2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    if eta_err / eta <= 0.02:
        color_eta = 'tab:green'
    else:
        color_eta = 'tab:red'

    ax.text(0.05 * x_ulim_cal, 0.59 * y_ulim_cal, 
            f'η: {round_it(eta, 3)}' + ' $\pm$ ' + f'{round_it(eta_err, 2)}',
            bbox = dict(facecolor = color_eta, alpha = 0.22, zorder = 3)) 

    ax2 = fig.add_axes([0.545,0.14,0.44,0.65])

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
    
    if ax2.get_legend_handles_labels() != ([], []):
        ax2.legend(loc = 'upper right')

    ax2.axvspan(vdr_region[0], vdr_region[1],
                alpha = 0.2, facecolor = 'tab:grey')

    m_llim = np.round(vdr_region[0], decimals = 2)
    m_ulim = np.round(vdr_region[1], decimals = 2)

    c_ray = 'tab:green'
        
    ax2.text(0.05 * x_ulim_vdr, 0.90 * y_ulim_vdr, 
             f'mol. cal. region: {m_llim} - {m_ulim} km',
            bbox = dict(facecolor = c_ray, alpha = 0.22, zorder = 3))  
    ax2.text(0.05 * x_ulim_vdr, 0.78 * y_ulim_vdr, 
            r'$δ_m$: '+f'{np.round(delta_m,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.66 * y_ulim_vdr, 
            r'$δ^{\star}$'+f': {np.round(delta_c_def,4)}' + ' $\pm$ ' + f'{round_it(delta_c_def_err,2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.54 * y_ulim_vdr, 
            r'$δ_{c}$'+f': {np.round(delta_c,4)}' + ' $\pm$ ' + f'{round_it(delta_c_err,2)}', 
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.42 * y_ulim_vdr, 
            r'$δ_{res}$: '+f'{np.round(delta_l,4)}' + ' $\pm$ ' + f'{round_it(delta_l_err,2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.30 * y_ulim_vdr, 
            r'$G_R$: '+f'{np.round(G_R,4)}, $G_T$: '+f'{round_it(G_T,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.18 * y_ulim_vdr, 
            r'$H_R$: '+f'{np.round(H_R,4)}, $H_T$: '+f'{round_it(H_T,4)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 
    ax2.text(0.05 * x_ulim_vdr, 0.06 * y_ulim_vdr, 
            r'$SR$ > '+f'{np.round(sr_lim,3)}, ' + r'$Δδ_p$ < ' + f'{np.round(err_p,decimals = 2)}',
            bbox = dict(facecolor = 'tab:cyan', alpha = 0.22, zorder = 3)) 

    fpath = os.path.join(dir_out, fname)
    
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    perform_color_reduction(color_reduction, fpath)
            
    return(fpath)

def make_filename(metadata, channel, meas_type, version, extra_type = '', extra_channel = None):
    
    if extra_channel == None:
        parts = [str(metadata['station_id']), str(metadata['lidar_id']), str(metadata['version_id']), str(metadata['config_id']), str(metadata['start_date']), str(metadata['start_time']), meas_type, str(channel), str(metadata['scc_channel_id'].loc[channel].values), str(extra_type), 'ATLAS', str(version)]
    else:
        parts = [str(metadata['station_id']), str(metadata['lidar_id']), str(metadata['version_id']), str(metadata['config_id']), str(metadata['start_date']), str(metadata['start_time']), meas_type, str(channel), str(metadata['scc_channel_id'].loc[extra_channel].values), str(extra_type), 'ATLAS', str(version)]
        
    fname = "_".join([part for part in parts if len(part) > 0])
    
    return(fname)

def make_filename_intercomparison(metadata_1, metadata_2, channel_1, channel_2, version, extra_type = ''):
        
    parts = [str(metadata_1['meas_id']), str(metadata_2['meas_id']), str(metadata_1['lidar_id']), str(metadata_2['lidar_id']), 'cmp', str(channel_1), str(channel_2), str(extra_type), 'ATLAS', str(version)]        

    fname = "_".join([part for part in parts if len(part) > 0])
    
    return(fname)

def get_plot_metadata(metadata, args, channel, meas_type, version, data_source_id = None):

    channel_d = dict(channel = channel)
    
    plot_metadata = dict()
    plot_metadata['processing_software'] = f"ATLAS_{version}"
    plot_metadata['measurement_type'] = meas_type
    
    if data_source_id == None:
        plot_metadata['atlas_channel_id'] = channel
    else:
        plot_metadata[f'atlas_channel_id_{data_source_id}'] = channel
    
    for key in metadata.keys():
        if data_source_id == None:
            key_new = key
        else:
            key_new = f'{key}_{data_source_id}'          
               
        if np.isscalar(metadata[key]):
            plot_metadata[key_new] = f"{metadata[key]}"
        elif type(metadata[key]) == type(DataArray()):
            if 'channel' in metadata[key].dims:
                plot_metadata[key_new] = f"{metadata[key].loc[channel_d].values}"
            
    for key in args.keys():
        plot_metadata[key] = f"{args[key]}"
    
    return(plot_metadata)

def add_plot_metadata(plot_path, plot_metadata, plot_metadata_extra=None):
    """
    Add metadata to a PNG file.

    Notes
    -----
    PNG metadata values must be text-like. Therefore:
    - None values are skipped
    - keys are converted to strings
    - values are converted to strings
    """

    im = Image.open(plot_path)
    meta = PngImagePlugin.PngInfo()

    def add_metadata_dict(metadata_dict):
        for key, value in metadata_dict.items():

            # Skip None values
            if value is None:
                continue

            # Skip empty lists if desired
            if value == []:
                continue

            meta.add_text(str(key), str(value))

    add_metadata_dict(plot_metadata)

    if plot_metadata_extra is not None:
        add_metadata_dict(plot_metadata_extra)

    im.save(plot_path, "png", pnginfo=meta)

    return None
    
def perform_color_reduction(color_reduction, plot_path):
    
    if color_reduction == True:
        im = Image.open(plot_path)
        im = im.convert('P', palette = Image.ADAPTIVE, colors = 255) 
        im.save(plot_path)

    return()

def export_plot(fig, args):
    
    dpi_val = args['dpi']

    dir_out = os.path.join(args['output_folder'],'plots') 

    fpath = os.path.join(dir_out, f"{args['filename']}.png")
            
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)


def clean_plots(plot_dir,pattern):
    
    if not os.path.isdir(plot_dir):
        return

    for filename in os.listdir(plot_dir):
        if pattern not in filename:
            continue

        file_path = os.path.join(plot_dir, filename)

        if os.path.isfile(file_path):
            os.remove(file_path)
            
def prepare_folder(caller_info, pattern):
    
    plot_dir = os.path.join(caller_info["output_folder"], "plots")
            
    os.makedirs(plot_dir, exist_ok = True)
            
    clean_plots(plot_dir, pattern = pattern)  

def collect_dict(data_list, data_keys, add_dicts = []):
    
    info = {}
    
    for i in range(len(data_keys)):
        info[data_keys[i]] = data_list[i]
    
    for d in add_dicts:
        info = d | info 
    
    return dict(sorted(info.items()))

def pass_to_args(args, data_list, data_keys):
    
    for i in range(len(data_keys)):
        args[data_keys[i]] = data_list[i]
        
    return(args)
    
def smoothing(args, x_vals, y_vals, err_type = "std"):
    
    if args['smooth'] and args['smoothing_window']:

        y_vals_sm, y_errs = \
            smooth_1D(y_vals = y_vals, 
                      x_vals = x_vals,
                      x_sm_lims = args['smoothing_range'],
                      x_sm_win = 1E3 * args['smoothing_window'],
                      expo = False,
                      err_type = err_type)
    
    else:
        y_vals_sm = y_vals.copy()
        y_errs = np.nan * y_vals.copy()
        
    return(y_vals_sm, y_errs)

def smoothing_2D(args, x_vals, y_vals, err_type = "std"):
    
    if args['smooth'] and args['smoothing_window']:
        
        if isinstance(args['smoothing_window'],list):
            smooth_2D = sliding_average_2D
        else:
            smooth_2D = sliding_average_2D_fast

        y_vals_sm, y_errs = smooth_2D(
            z_vals = y_vals, 
            y_vals = x_vals,
            y_sm_lims = args['smoothing_range'],
            y_sm_win = args['smoothing_window'],
            expo = args['smooth_exponential']
            )
    else:
        y_vals_sm = y_vals
        y_errs = np.nan * np.zeros_like(y_vals)
        
    return(y_vals_sm, y_errs)

def slice_by_vertical_scale(
    da,
    vertical_scale,
    x_lims,
    time_dim="time",
    bin_dim="bins",
):
    """
    Slice a DataArray along bins using vertical_scale values.

    Works for:
    - 1D arrays with dims ("bins",)
    - 2D arrays with dims ("time", "bins")

    Bins are kept if:
    - vertical_scale is finite and within x_lims
    - data are valid:
        - 1D: bin value is not NaN
        - 2D: bin has at least one non-NaN value over time
    """

    if bin_dim not in da.dims:
        raise ValueError(f"da must contain the '{bin_dim}' dimension.")

    if bin_dim not in vertical_scale.dims:
        raise ValueError(f"vertical_scale must contain the '{bin_dim}' dimension.")

    if time_dim in da.dims:
        valid_bins = da.notnull().any(dim=time_dim)
    else:
        valid_bins = da.notnull()

    scale_mask = (
        vertical_scale.notnull()
        & (vertical_scale >= x_lims[0])
        & (vertical_scale <= x_lims[1])
    )

    mask = valid_bins & scale_mask

    # Required when da/vertical_scale are Dask-backed and drop=True is used
    mask = mask.compute()

    da_sliced = da.where(mask, drop=True)
    vertical_scale_sliced = vertical_scale.where(mask, drop=True)

    return da_sliced, vertical_scale_sliced, mask

def add_extra_plot_metadata(plot_metadata, norm_region_flag, 
                            stats_norm_region, maximum_channel_height):
    
    plot_metadata['norm_region_flag'] = f"{norm_region_flag}"
    for key in stats_norm_region.keys():
        plot_metadata[f"stats_{key}"] = f"{stats_norm_region[key]}"
    for key in stats_norm_region.keys():
        plot_metadata[f"masks_{key}"] = f"{stats_norm_region[key]}"
        
    plot_metadata['maximum_channel_height'] = f"{maximum_channel_height}"
    
    return(plot_metadata)

def convert_m_to_km(x):
   
    return 1E-3 * x

def get_vertical_axis_label(vertical_scale):
   
    # Convert meters to kilometers and select ranges or heights for the x axis depending on the use_dis value 
    if vertical_scale == 'range':
        x_label = "Range from the lidar [km]"
        
    elif vertical_scale == 'height_agl':
        x_label = "Height [km asl]"
        
    elif vertical_scale == 'height_asl':
        x_label = "Height [km agl]"
        
    else:
        raise Exception("Vertical scale '{vertical_scale}' not supported")

    return x_label 

def multiply_y_values(sig, sig_err, coef):
    
    # Multiply the   
    y_vals  = coef * sig.copy()
    
    y_errs = coef * sig_err.copy()
    
    return(y_vals, y_errs)

def insert_nan_time_gaps(sig, dim="time", gap_factor=1.5):
    """
    Insert NaN profiles around large time gaps while keeping original time values.

    Parameters
    ----------
    sig : xr.DataArray
        Input DataArray with a time dimension.
    dim : str
        Name of the time dimension.
    gap_factor : float
        A gap is detected when the time difference between two adjacent profiles
        is larger than gap_factor * minimum_time_difference.

    Returns
    -------
    sig_gap : xr.DataArray
        DataArray with extra NaN profiles inserted if gaps were detected.
        If no gaps were detected, the original sorted array is returned.
    time_res : pandas.Timedelta or None
        Inferred minimum temporal resolution.
    changed : bool
        True if NaN profiles were inserted, False otherwise.
    """

    sig = sig.sortby(dim)

    time = pd.DatetimeIndex(sig[dim].values)

    if len(time) < 2:
        return sig, None, False

    dt = time.to_series().diff().dropna()
    time_res = dt.min()

    gap_limit = gap_factor * time_res

    extra_times = []

    for t0, t1 in zip(time[:-1], time[1:]):
        if (t1 - t0) > gap_limit:
            extra_times.append(t0 + time_res)
            extra_times.append(t1 - time_res)

    if len(extra_times) == 0:
        return sig, time_res, False

    extra_times = pd.DatetimeIndex(extra_times)

    # Keep only valid extra times inside the measurement period
    extra_times = extra_times[
        (extra_times > time[0]) &
        (extra_times < time[-1])
    ]

    # Remove possible duplicates
    extra_times = extra_times.difference(time)

    if len(extra_times) == 0:
        return sig, time_res, False

    new_time = time.union(extra_times).sort_values()

    sig_gap = sig.reindex({dim: new_time})

    return sig_gap, time_res, True

def slice_time(da, t_lims, dim="time"):
    """
    Slice an xarray DataArray/Dataset along time and report if slicing happened.

    Parameters
    ----------
    da : xr.DataArray or xr.Dataset
        Input data with a time coordinate.
    t_lims : list
        [] or None means no slicing.
        Otherwise: ["yyyymmdd_hhmm", "yyyymmdd_hhmm"].
    dim : str
        Name of the time dimension.

    Returns
    -------
    da_out : xr.DataArray or xr.Dataset
        Sliced or original data.
    sliced : bool
        True if the returned array was actually sliced.
        False if the original array was returned.
    """

    if t_lims is None or len(t_lims) == 0:
        return da, False

    if len(t_lims) != 2:
        raise ValueError("t_lims must be [], None, or [start, end].")

    da_sorted = da.sortby(dim)

    time_min = pd.Timestamp(da_sorted[dim].values[0])
    time_max = pd.Timestamp(da_sorted[dim].values[-1])

    t_start = pd.to_datetime(t_lims[0], format="%Y%m%d_%H%M")
    t_end = pd.to_datetime(t_lims[1], format="%Y%m%d_%H%M")

    if t_start > t_end:
        raise ValueError("The start time in t_lims must be before the end time.")

    # Requested limits fully include the data, so nothing needs to be sliced
    if t_start <= time_min and t_end >= time_max:
        return da_sorted, False

    # Requested limits do not overlap the data at all, so keep original data
    if t_end < time_min or t_start > time_max:
        return da_sorted, False

    da_sliced = da_sorted.sel({dim: slice(t_start, t_end)})

    # If the number of time entries did not change, consider it untouched
    sliced = da_sliced.sizes[dim] != da_sorted.sizes[dim]

    return da_sliced, sliced

def add_fitting_suptitle(
    fig,
    title,
    y=1.,
    max_fontsize=12,
    min_fontsize=7,
    margin=0.04,
):
    """
    Add a suptitle and shrink its font size if it is too wide.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    title : str
        Title text.
    y : float
        Vertical position of the title in figure coordinates.
    max_fontsize : int
        Starting font size.
    min_fontsize : int
        Smallest allowed font size.
    margin : float
        Fractional horizontal margin on each side.

    Returns
    -------
    title_obj : matplotlib.text.Text
        The title text object.
    """

    title_obj = fig.suptitle(title, y=y, fontsize=max_fontsize)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    fig_width = fig.bbox.width
    allowed_width = fig_width * (1.0 - 2.0 * margin)

    fontsize = max_fontsize

    while fontsize > min_fontsize:
        title_width = title_obj.get_window_extent(renderer=renderer).width

        if title_width <= allowed_width:
            break

        fontsize -= 1
        title_obj.set_fontsize(fontsize)

        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

    return title_obj