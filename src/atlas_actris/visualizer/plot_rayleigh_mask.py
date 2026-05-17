#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 21:31:32 2025

@author: nikos
"""

import numpy as np
from matplotlib import pyplot as plt
from visualizer.plot_utils import export_plot

def generate_plot(args, masks):
    
    # [X, Y] = np.meshgrid(masks['total'].middle_point.values, masks['total'].window)
    fig = plt.figure(figsize=(12. , 8.))

    fig.suptitle(args["title"])
    
    # fig_x = 0.44
    # fig_y = 0.16
    
    # fig_edg1_x = 0.06
    # fig_edg2_x = 0.54
    
    # fig_edg1_y = 0.72
    # fig_edg2_y = 0.50
    # fig_edg3_y = 0.28
    # fig_edg4_y = 0.06

    # ax1_coords = [fig_edg1_x, fig_edg1_y, fig_x, fig_y]
    # ax2_coords = [fig_edg2_x, fig_edg1_y, fig_x, fig_y]
    # ax3_coords = [fig_edg1_x, fig_edg2_y, fig_x, fig_y]
    # ax4_coords = [fig_edg2_x, fig_edg2_y, fig_x, fig_y]
    # ax5_coords = [fig_edg1_x, fig_edg3_y, fig_x, fig_y]
    # ax6_coords = [fig_edg2_x, fig_edg3_y, fig_x, fig_y]
    # ax7_coords = [fig_edg1_x, fig_edg4_y, fig_x, fig_y]
    # ax8_coords = [fig_edg2_x, fig_edg4_y, fig_x, fig_y]
    # ax8_coords = [fig_edg2_x, fig_edg4_y, fig_x, fig_y]
    
    fig_x = 0.28
    fig_y = 0.23
    
    fig_edg1_x = 0.06
    fig_edg2_x = 0.38
    fig_edg3_x = 0.70
    
    fig_edg1_y = 0.64
    fig_edg2_y = 0.35
    fig_edg3_y = 0.06
    
    ax1_coords = [fig_edg1_x, fig_edg1_y, fig_x, fig_y]
    ax2_coords = [fig_edg2_x, fig_edg1_y, fig_x, fig_y]
    ax3_coords = [fig_edg3_x, fig_edg1_y, fig_x, fig_y]
    ax4_coords = [fig_edg1_x, fig_edg2_y, fig_x, fig_y]
    ax5_coords = [fig_edg2_x, fig_edg2_y, fig_x, fig_y]
    ax6_coords = [fig_edg3_x, fig_edg2_y, fig_x, fig_y]
    ax7_coords = [fig_edg1_x, fig_edg3_y, fig_x, fig_y]
    ax8_coords = [fig_edg2_x, fig_edg3_y, fig_x, fig_y]
    ax9_coords = [fig_edg3_x, fig_edg3_y, fig_x, fig_y]

    X = masks['total'].middle_point.values
    Y = masks['total'].window.values
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax1_coords,
                     X = X,
                     Y = Y,
                     Z = masks['first_derivative'].values,
                     title = 'Derivative mask',
                     use_y_label = True,
                     args = args)
 
    plot_single_mask(fig = fig, 
                     ax_coords = ax2_coords,
                     X = X,
                     Y = Y,
                     Z = masks['relative_sem'].values,
                     title = 'Relative SEM mask',
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax3_coords,
                     X = X,
                     Y = Y,
                     Z = masks['second_derivative'].values,
                     title = 'Second derivative mask',
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax4_coords,
                     X = X,
                     Y = Y,
                     Z = masks['shapiro_wilk'].values,
                     title = 'Shapiro-Wilk mask',
                     use_y_label = True,
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax5_coords,
                     X = X,
                     Y = Y,
                     Z = masks['cross_criterion'].values,
                     title = 'Cross-check mask',
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax6_coords,
                     X = X,
                     Y = Y,
                     Z = masks['is_positive'].values,
                     title = 'Positive signal mask',
                     args = args)

    plot_single_mask(fig = fig, 
                     ax_coords = ax7_coords,
                     X = X,
                     Y = Y,
                     Z = masks['durbin_watson'].values,
                     title = 'Durbin-Watson mask',
                     use_x_label = True,
                     use_y_label = True,
                     args = args)

    plot_single_mask(fig = fig, 
                     ax_coords = ax8_coords,
                     X = X,
                     Y = Y,
                     Z = masks['residual_extinction'].values,
                     title = 'Residual extinction mask',
                     use_x_label = True,
                     use_y_label = True,
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax9_coords,
                     X = X,
                     Y = Y,
                     Z = masks['total'].values,
                     title = 'Combined Mask',
                     use_x_label = True,
                     use_y_label = True,
                     args = args)
    
    fpath = export_plot(fig, args)
            
    return(fpath)

def plot_single_mask(fig, ax_coords, X, Y, Z, title, args, 
                     use_x_label = False, use_y_label = False):
    
    x_llim = args['x_lims'][0]
    x_ulim = args['x_lims'][1]
    y_llim = args['fit_mask_window'][0]
    y_ulim = args['fit_mask_window'][1]
    
    win = args['fit_mask_window_step']
    
    # extent = (args['fit_mask_region'][0], 
    #           args['fit_mask_region'][1], 
    #           args['fit_mask_window'][0], 
    #           args['fit_mask_window'][1])
    
    extent = (X[0], X[-1], Y[0], Y[-1])
    
    y_ticks = np.arange(y_llim, y_ulim + win, win * 10.)

    ax = fig.add_axes(ax_coords)
    # ax.pcolormesh(X, Y, Z, vmin = 0, vmax = 1)
    ax.imshow(Z, extent = extent, cmap = 'viridis', interpolation = 'nearest',
              origin = "lower", aspect = "auto")
    
    ax.set_title(title, pad = 3)
    
    if use_x_label:
        ax.set_xlabel('Window center [km]')

    if use_y_label:
        ax.set_ylabel('Window size [km]')
    
    ax.set_ylim([y_llim, y_ulim])
    ax.set_xlim([x_llim, x_ulim])
    
    ax.set_yticks(y_ticks)
    
    return(ax)