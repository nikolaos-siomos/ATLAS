#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 18:39:28 2025

@author: nikos
"""


import numpy as np
from utils.toolbox import round_it
from matplotlib import pyplot as plt
from visualizer.plot_utils import export_plot
from matplotlib.ticker import MultipleLocator
from visualizer.generate_rayleigh_fit_axis import get_rayleigh_fit_y_limits
from visualizer.plot_utils import get_vertical_axis_label

def generate_plot(X, Y1, Y2, Y1E, args):

    # Left panel coordinates
    ax1_coords = [0.05, 0.145, 0.52, 0.69]

    # Right panel coordinates
    ax2_coords = [0.625,0.145,0.36,0.69]
    
    # Create the figure
    fig = plt.figure(figsize=(15 , 3))

    fig.suptitle(args["title"])
    
    left_panel(fig = fig, ax1_coords = ax1_coords, 
               X = X, Y1 = Y1, Y1E = Y1E, Y2 = Y2, 
               args = args)
    
    right_panel(fig = fig, ax2_coords = ax2_coords, 
                X = X, Y1 = Y1, Y1E = Y1E, Y2 = Y2, 
                args = args)
    
    fpath = export_plot(fig, args)
        
    return(fpath)

def left_panel(fig, ax1_coords, X, Y1, Y1E, Y2, args):
    
    ax = fig.add_axes(ax1_coords)
        
    ax.plot(X, Y1, color = 'tab:blue', label = 'measured')
    ax.plot(X, Y2, color = 'tab:red', label = 'molecular')

    if np.isnan(Y1E).all() == False:
        ax.fill_between(X, Y1 - Y1E, Y1 + Y1E, color = 'tab:blue', alpha = 0.3)

    y1_label = get_y1_label()
    
    x_ticks_1, x_tick_labels_1 = get_x_ticks_1(x_lims = args['x_lims'], 
                                               x_tick = args['x_tick'])


    ax.set_xticks(x_ticks_1, labels = x_tick_labels_1)
    
    ax.set_xlim([args["x_lims"][0], args["x_lims"][1]])
    
    # Get the x axis label depending on use_dis
    x_label = get_vertical_axis_label(args['vertical_scale'])
    ax.set_xlabel(x_label)
    
    x_tick = args['x_tick']
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.))

    # Get the limits of the y axis in case they are not provided by the user (default)
    y_lims = get_rayleigh_fit_y_limits(
        y1_vals = Y1, 
        y2_vals = Y2, 
        y_lims = args['y_lims'], 
        wavelength = args['detected_wavelength'], 
        use_lin_scale = args['use_lin_y_scale']
        )
    
    ax.set_ylim([y_lims[0], y_lims[1]])
    ax.set_ylabel(y1_label)
    
    use_lin_scale = args['use_lin_y_scale']
    
    if use_lin_scale == False:
        ax.set_yscale('log')

    ax.grid(which = 'both')
    
    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc = 'lower left')

    ax.axvspan(args['norm_region'][0], 
               args['norm_region'][1], 
               alpha = 0.2, facecolor = 'tab:grey')
    
    box_colors = \
        mask_to_color(masks_norm_region = args["masks_norm_region"],
                      norm_region_flag = args['norm_flag'])
        
    box_edge_x, box_edge_y = box_edges(x_ulim = args['x_lims'][1], 
                                       y_ulim = y_lims[1], 
                                       use_lin_scale = use_lin_scale)
    
    box_text = get_box_text(norm_region = args['norm_region'], 
                            stats_norm_region = args["stats_norm_region"])
    
    ax = add_text_ax1(ax, 
                      box_colors = box_colors, 
                      box_edge_x = box_edge_x, 
                      box_edge_y = box_edge_y, 
                      box_text = box_text)
    
    return(ax)

def right_panel(fig, ax2_coords, X, Y1, Y1E, Y2, args):
    
    x_ticks_2, x_tick_labels_2 = get_x_ticks_2(x_lims = args['x_lims'], 
                                               x_tick = args['x_tick'])
    
    ax2 = fig.add_axes(ax2_coords)
    
    if np.isnan(Y1E).all() == False:
        ax2.fill_between(X, (Y1 - Y1E - Y2) / Y2, 
                         (Y1 + Y1E - Y2) / Y2, color = 'tab:blue', 
                         alpha = 0.3, label = 'sem')
        
    ax2.plot(X, (Y1 - Y2) / Y2, color = 'tab:blue',label = 'mean')
    
    ax2.axhline(c = 'k')
    
    y2_label = get_y2_label()
    
    ax2.set_xticks(x_ticks_2, labels = x_tick_labels_2)
    ax2.set_xlim([args['x_lims'][0], args['x_lims'][1]])
    
    # Get the x axis label depending on use_dis
    x_label = get_vertical_axis_label(args['vertical_scale'])
    ax2.set_xlabel(x_label)
    ax2.xaxis.set_minor_locator(MultipleLocator(args['x_tick']))
    
    y_ticks = np.round(np.arange(-0.40, 0.40 + 0.10, 0.10), decimals = 2)
    ax2.set_yticks(y_ticks, labels = ["%.2f" % tick for tick in y_ticks])
    ax2.set_ylim([y_ticks[0], y_ticks[-1]])
    ax2.set_ylabel(y2_label)
    
    ax2.grid(which = 'both')
    
    ax2.axvspan(args['norm_region'][0], 
                args['norm_region'][1], 
                alpha = 0.2, facecolor = 'tab:grey')
    
    return(ax2)

def get_y1_label():
    
    # Get the y axis labels
    y_label = 'Attn. Bsc. rel. to fit range [$m^{-1} sr^{-1}$]'
    
    return(y_label)


def get_y2_label():
    
    # Get the y axis labels
    y_label = 'Relative Diff.'
    
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
    
    x_tick_labels_1 = [str(float(tick)).rstrip('0').rstrip('.') for tick in x_ticks_1]
    
    return(x_ticks_1, x_tick_labels_1)

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
    
    x_tick_labels_2 = [str(float(tick)).rstrip('0').rstrip('.') for tick in x_ticks_2]

    return(x_ticks_2, x_tick_labels_2)


def add_text_ax1(ax, box_colors, box_edge_x, box_edge_y, box_text):
 


    for key in ["norm_region", "max_channel_height"]:
        ax.text(box_edge_x['low_mid'], box_edge_y[key], 
                box_text[key],
                transform = ax.transAxes,
                bbox = dict(facecolor = box_colors[key], alpha = 0.22, zorder = 3))
        
    for key in ["relative_sem", "first_derivative"]:
        ax.text(box_edge_x['column_1'], box_edge_y[key], 
                box_text[key],
                transform = ax.transAxes,
                bbox = dict(facecolor = box_colors[key], alpha = 0.22, zorder = 3))
        
    for key in ["second_derivative", "shapiro_wilk","cross_criterion"]:
        ax.text(box_edge_x['column_2'], box_edge_y[key], 
                box_text[key],
                transform = ax.transAxes,
                bbox = dict(facecolor = box_colors[key], alpha = 0.22, zorder = 3))

    for key in ["is_positive", "durbin_watson", "residual_extinction"]:
        ax.text(box_edge_x['column_3'], box_edge_y[key], 
                box_text[key],
                transform = ax.transAxes,
                bbox = dict(facecolor = box_colors[key], alpha = 0.22, zorder = 3))
                
    return(ax)


def mask_to_color(masks_norm_region, norm_region_flag):
    
    if norm_region_flag == 'auto': 
        c_nrmg = 'tab:blue'
        c_maxh = 'tab:blue'
        c_norm = 'tab:green'
        
    elif norm_region_flag == 'external': 
        c_nrmg = 'tab:orange'
        c_maxh = 'tab:orange'
        c_norm = 'tab:red'
        
    else: 
        c_nrmg = 'tab:orange'
        c_maxh = 'tab:orange'
        c_norm = 'tab:red'
        
        
    if masks_norm_region['relative_sem']: 
        c_msem = 'tab:green'
    else: 
        c_msem = 'tab:orange'
        
    if masks_norm_region['first_derivative']: 
        c_mder = 'tab:green'
    else: 
        c_mder = 'tab:red'
    
    if masks_norm_region['second_derivative']: 
        c_msec = 'tab:green'
    else: 
        c_msec = 'tab:red'
        
    if masks_norm_region['shapiro_wilk']: 
        c_mshp = 'tab:green'
    else: 
        c_mshp = 'tab:red'
    
    if masks_norm_region['cross_criterion']: 
        c_mcrc = 'tab:green'
    else: 
        c_mcrc = 'tab:red'

    if masks_norm_region['is_positive']: 
        c_mpos = 'tab:green'
    else: 
        c_mpos = 'tab:red'

    if masks_norm_region['durbin_watson']: 
        c_mdbw = 'tab:green'
    else: 
        c_mdbw = 'tab:red'

    if masks_norm_region['residual_extinction']: 
        c_mext = 'tab:green'
    else: 
        c_mext = 'tab:red'
        
    box_colors = {'norm_region': c_nrmg,
                  'max_channel_height': c_maxh,
                  'relative_sem':c_msem,
                  'first_derivative':c_mder,
                  'second_derivative':c_msec,
                  'shapiro_wilk':c_mshp,
                  'cross_criterion':c_mcrc,
                  'is_positive':c_mpos,
                  'durbin_watson':c_mdbw,
                  'normalization_factor':c_norm,
                  'residual_extinction':c_mext}
        
    return(box_colors)

def get_box_text(norm_region, stats_norm_region):
        
    n_llim = np.round(norm_region[0], decimals = 2)
    n_ulim = np.round(norm_region[1], decimals = 2)

    max_ch_h = np.round((norm_region[1] + norm_region[0]) / 2., decimals = 2)
    
    rsem = stats_norm_region['relative_sem']
    rslope = stats_norm_region['first_derivative']
    
    box_text = {
        'norm_region':f'norm. window: {n_llim} - {n_ulim} km',
        'max_channel_height':f'max ch height: {max_ch_h} km',
        'relative_sem':f'rsem: {round_it(rsem, 3)}',
        'first_derivative':f'rslope: {round_it(rslope, 3)}',
        'second_derivative':'Curvature',
        'shapiro_wilk':'Shapiro-Wilk',
        'cross_criterion':'Cross crit',
        'is_positive':'Positive sig',
        'durbin_watson':'Durbin Watson',
        'residual_extinction':'Residual Extinction'
        }
    
    return(box_text)

def box_edges(x_ulim, y_ulim, use_lin_scale):
    
    box_edge_x = {
        'column_1':0.45,
        'column_2':0.63,
        'column_3':0.80,
        'low_mid': 0.20
        }
    
    # box_edge_y = {
    #     'norm_region':0.90,
    #     'max_channel_height':0.77,
    #     'relative_sem':0.64,
    #     'first_derivative':0.90,
    #     'second_derivative':0.77,
    #     'shapiro_wilk':0.64,
    #     'cross_criterion':0.90,
    #     'is_positive':0.77,
    #     'durbin_watson':0.64
    #     'durbin_watson':0.64
    #     }   

    box_edge_y = {
        'norm_region':0.06,
        'max_channel_height':0.20,
        'relative_sem':0.90,
        'first_derivative':0.77,
        'second_derivative':0.90,
        'shapiro_wilk':0.77,
        'cross_criterion':0.64,
        'is_positive':0.90,
        'durbin_watson':0.77,
        'residual_extinction':0.64
        }       
        
    return(box_edge_x, box_edge_y)