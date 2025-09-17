#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 18:39:28 2025

@author: nikos
"""

import os
import numpy as np
from matplotlib import pyplot as plt
from plot_utils import round_it
from matplotlib.ticker import MultipleLocator

def generate_plot(dir_out, fname, title, norm_region, norm_region_flag, 
                  X, Y1, Y2, Y1E, x_lims, y_lims, 
                  stats_norm_region, masks_norm_region, args):
    
    # Unpack parameters from the settings file
    dpi_val = args['dpi']
    use_lin = args['use_lin_scale']
    x_tick = args['x_tick']

    # Create the figure
    fig = plt.figure(figsize=(15 , 3))

    fig.suptitle(title)

    ax = fig.add_axes([0.05, 0.145, 0.52, 0.69])
        
    ax.plot(X, Y1, color = 'tab:blue', label = 'measured')
    ax.plot(X, Y2, color = 'tab:red', label = 'molecular')

    if np.isnan(Y1E).all() == False:
        ax.fill_between(X, Y1 - Y1E, Y1 + Y1E, color = 'tab:blue', alpha = 0.3)

    x_label = get_x_label(use_dis = args['use_dis'])
    y_label = get_y_label()
    
    x_ticks_1 = get_x_ticks_1(x_lims = args['x_lims'], 
                              x_tick = args['x_tick'])
    x_ticks_2 = get_x_ticks_2(x_lims = args['x_lims'], 
                              x_tick = args['x_tick'])


    ax.set_xticks(x_ticks_1, labels = x_ticks_1)
    ax.set_xlim([x_lims[0], x_lims[1]])
    ax.set_xlabel(x_label)
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.))

    ax.set_ylim([y_lims[0], y_lims[1]])
    ax.set_ylabel(y_label)
    if use_lin == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    
    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc = 'lower left')

    ax.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')
    
    box_colors = \
        mask_to_color(msem_c = msem_c, 
                      mder_c = mder_c, 
                      msec_c = msec_c, 
                      mshp_c = mshp_c, 
                      mcrc_c = mcrc_c,
                      norm_region_flag = norm_region_flag)
        
    box_edge_x, box_edge_y = box_edges(x_ulim = x_lims[1], 
                                       y_ulim = y_lims[1], 
                                       use_lin = use_lin)
    
    box_text = get_box_text(norm_region = norm_region, rsem = rsem, rslope = rslope)
    
    ax = add_text_ax1(ax, 
                      box_colors = box_colors, 
                      box_edge_x = box_edge_x, 
                      box_edge_y = box_edge_y, 
                      box_text = box_text)

    ax2 = fig.add_axes([0.625,0.145,0.36,0.69])
    
    if np.isnan(Y1E).all() == False and np.isnan(Y2E).all() == True:
        ax2.fill_between(X, (Y1 - Y1E - Y2) / Y2, 
                         (Y1 + Y1E - Y2) / Y2, color = 'tab:blue', 
                         alpha = 0.3, label = 'sem')
        ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue',label = 'mean')

    if np.isnan(Y1E).all() == False and np.isnan(Y2E).all() == False:
        ax2.fill_between(X, (Y1 - Y1E - Y2) / Y2, 
                         (Y1 + Y1E - Y2) / Y2, color = 'tab:blue', 
                         alpha = 0.3, label = 'sem1')
        ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue',label = 'mean')
        ax2.fill_between(X, (Y1 - Y2E - Y2) / Y2, 
                         (Y1 + Y2E - Y2) / Y2, color = 'tab:red', 
                         alpha = 0.3, label = 'sem2')
        
        if ax2.get_legend_handles_labels() != ([], []):
            ax2.legend(loc = 'lower left')

    else:
        ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue')
        
    
    ax2.axhline(c = 'k')

    ax2.set_xticks(x_ticks_2, labels = x_ticks_2)
    ax2.set_xlim([x_lims[0], x_lims[1]])
    ax2.set_xlabel(x_label)
    ax2.xaxis.set_minor_locator(MultipleLocator(x_tick))

    y_ticks = np.round(np.arange(-0.40, 0.40 + 0.10, 0.10), decimals = 2)
    ax2.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax2.set_ylim([y_ticks[0], y_ticks[-1]])
    ax2.set_ylabel('Relative Diff.')
    
    ax2.grid(which = 'both')
    
    ax2.axvspan(norm_region[0], norm_region[1], alpha = 0.2, facecolor = 'tab:grey')
    
    fpath = os.path.join(dir_out, fname)
            
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
        
    return(fpath)

def get_x_label(use_dis):

    # Get the x axis label
    if use_dis:
        x_label = 'Range above the lidar [km]'
        
    else:
        x_label = 'Height above the lidar [km]'
    
    return(x_label)

def get_y_label():
    
    # Get the y axis labels
    y_label = 'Attn. Bsc. rel. to fit range [$m^{-1} sr^{-1}$]'
    
    return(y_label)

def get_x_ticks_1(x_lims, x_tick):

    x_llim = x_lims[0]
    x_ulim = x_lims[1]
    
    # Check if the given x axis limits are compatible with the x axis tick
    if x_tick >= x_ulim - x_llim:
        raise Exception(f"The region between the provided x axis limits ({x_llim} to {x_ulim} km) is shorter than the x_tick parameter ({x_tick} km)")
    
    # Get the ticks for the first plot
    x_ticks_1 = np.arange(x_tick * np.ceil(x_llim / x_tick), 
                          x_tick * (np.floor(x_ulim / x_tick) + 1.), 
                          x_tick)
    
    if np.abs(x_llim - x_ticks_1[0]) < x_tick * 0.25:
        x_ticks_1[0] = x_llim
    else:
        x_ticks_1 = np.hstack((x_llim, x_ticks_1))

    if np.abs(x_ulim - x_ticks_1[-1]) < x_tick * 0.25:
        x_ticks_1[-1] = x_ulim
    else:
        x_ticks_1 = np.hstack((x_ticks_1, x_ulim))

    x_ticks_1 = np.round(x_ticks_1, decimals = 2)
    
    return(x_ticks_1)

def get_x_ticks_2(x_lims, x_tick):

    x_llim = x_lims[0]
    x_ulim = x_lims[1]
    
    x_tick_2 = 2. * x_tick 

    # Check if the given x axis limits are compatible with the x axis tick
    if x_tick_2 >= x_ulim - x_llim:
        raise Exception(f"The region between the provided x axis limits ({x_llim} to {x_ulim} km) is shorter than two times the x_tick parameter ({x_tick_2} km)")
    
    # Get the ticks for the second plot
    x_tick_2 = 2. * x_tick 
    x_ticks_2 = np.arange(x_tick_2 * np.floor(x_llim / x_tick_2), 
                          x_tick_2 * (np.ceil(x_ulim / x_tick_2) + 1.), 
                          x_tick_2)
    
    return(x_ticks_2)


def add_text_ax1(ax, box_colors, box_edge_x, box_edge_y, box_text):

    for key in ["norm_region", "rsem", "rslope"]:
        ax.text(box_edge_x['column_1'], box_edge_y[key], 
                box_text[key],
                bbox = dict(facecolor = box_colors[key], alpha = 0.22, zorder = 3))
        
    for key in ["curvature", "gaussian_noise", "cross_crit"]:
        ax.text(box_edge_x['column_1'], box_edge_y[key], 
                box_text[key],
                bbox = dict(facecolor = box_colors[key], alpha = 0.22, zorder = 3))
        
    return(ax)


def mask_to_color(msem_c, mder_c, msec_c, mshp_c, mcrc_c, norm_region_flag):
    
    if msem_c: c_msem = 'tab:yellow'
    else: c_msem = 'tab:green'
        
    if mder_c: c_mder = 'tab:green'
    else: c_mder = 'tab:red'
    
    if norm_region_flag == 'auto': c_norm = 'tab:green'
    elif norm_region_flag == 'external': c_norm = 'tab:yellow'
    else: c_norm = 'tab:red'
        
    if msec_c: c_msec = 'tab:green'
    else: c_msec = 'tab:red'
        
    if mshp_c: c_mshp = 'tab:green'
    else: c_mshp = 'tab:red'
    
    if mcrc_c: c_mcrc = 'tab:green'
    else: c_mcrc = 'tab:red'
    
    box_colors = {'sem':c_msem,
                  'der':c_mder,
                  'sec':c_msec,
                  'shp':c_mshp,
                  'crc':c_mcrc,
                  'norm':c_norm}
        
    return(box_colors)

def get_box_text(norm_region, rsem, rslope):
        
    n_llim = np.round(norm_region[0], decimals = 2)
    n_ulim = np.round(norm_region[1], decimals = 2)
    
    box_text = {'normalization_region':f'norm. region: {n_llim} - {n_ulim} km',
                'rsem':f'rsem: {round_it(rsem, 3)}',
                'rslope':f'rslope: {round_it(rslope, 3)}',
                'curvature':'Curvature',
                'gaussian_noise':'Gaussian noise',
                'cross_crit':'Cross crit'}
    
    return(box_text)

def box_edges(x_ulim, y_ulim, use_lin):
    
    box_edge_x = {'column_1':0.55 * x_ulim,
                  'column_2':0.75 * x_ulim}
    
    if use_lin == False:
        box_edge_y = {'normalization_region':0.60 * y_ulim,
                      'rsem':0.30 * y_ulim,
                      'rslope':0.15 * y_ulim}       
    else:
        box_edge_y = {'curvature':0.90 * y_ulim,
                      'gaussian_noise':0.82 * y_ulim,
                      'cross_crit':0.74 * y_ulim} 
        
    return(box_edge_x, box_edge_y)