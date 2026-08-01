#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 21:31:32 2025

@author: nikos
"""

import numpy as np
from matplotlib import pyplot as plt
from ..plotting.plot_utils import export_plot
from matplotlib.colors import ListedColormap

def generate_plot(args, masks):
    
    fig = plt.figure(figsize=(12. , 6.))

    fig.suptitle(args["title"])
    
    fig_x = 0.28
    fig_y = 0.33
    
    fig_edg1_x = 0.06
    fig_edg2_x = 0.38
    fig_edg3_x = 0.70
    
    fig_edg1_y = 0.53
    fig_edg2_y = 0.08
    
    ax1_coords = [fig_edg1_x, fig_edg1_y, fig_x, fig_y]
    ax2_coords = [fig_edg2_x, fig_edg1_y, fig_x, fig_y]
    ax3_coords = [fig_edg3_x, fig_edg1_y, fig_x, fig_y]
    ax4_coords = [fig_edg1_x, fig_edg2_y, fig_x, fig_y]
    ax5_coords = [fig_edg2_x, fig_edg2_y, fig_x, fig_y]
    ax6_coords = [fig_edg3_x, fig_edg2_y, fig_x, fig_y]

    X = masks['total'].middle_point.values
    Y = masks['total'].window.values
    
    plot_single_mask(
        fig = fig, 
        ax_coords = ax1_coords,
        X = X,
        Y = Y,
        Z = masks['first_derivative'].values,
        title = 'Derivative mask',
        use_y_label = True,
        args = args
        )
   
    plot_single_mask(
        fig = fig, 
        ax_coords = ax2_coords,
        X = X,
        Y = Y,
        Z = masks['second_derivative'].values,
        title = 'Second derivative mask',
        args = args
        )
    
    plot_single_mask(
        fig = fig, 
        ax_coords = ax3_coords,
        X = X,
        Y = Y,
        Z = masks['relative_sem'].values,
        title = 'Relative SEM mask',
        args = args
        )
   
    plot_single_mask(
        fig = fig, 
        ax_coords = ax4_coords,
        X = X,
        Y = Y,
        Z = masks['shapiro_wilk'].values,
        title = 'Shapiro-Wilk mask',
        use_y_label = True,
        args = args
        )

    plot_single_mask(
        fig = fig, 
        ax_coords = ax5_coords,
        X = X,
        Y = Y,
        Z = masks['durbin_watson'].values,
        title = 'Durbin-Watson mask',
        use_x_label = True,
        use_y_label = True,
        args = args
        )
   
    plot_single_mask_with_overlay(
        fig = fig, 
        ax_coords = ax6_coords,
        X = X,
        Y = Y,
        keep_mask = masks['total_pruned'].values,
        removed_mask = masks['rejected'].values,
        title = 'Combined Mask',
        use_x_label = True,
        use_y_label = True,
        args = args
        )
    
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
        ax.set_xlabel('Window center [bins]')

    if use_y_label:
        ax.set_ylabel('Window size [bins]')
    
    ax.set_ylim([y_llim, y_ulim])
    ax.set_xlim([x_llim, x_ulim])
    
    ax.set_yticks(y_ticks)
    
    return(ax)

def plot_single_mask_with_overlay(
    fig, ax_coords, X, Y,
    keep_mask, removed_mask,
    title, args,
    use_x_label=False, use_y_label=False,
    removed_alpha=0.9,
):
    x_llim = args['x_lims'][0]
    x_ulim = args['x_lims'][1]
    y_llim = args['fit_mask_window'][0]
    y_ulim = args['fit_mask_window'][1]

    win = args['fit_mask_window_step']
    extent = (X[0], X[-1], Y[0], Y[-1])
    y_ticks = np.arange(y_llim, y_ulim + win, win * 10.)

    ax = fig.add_axes(ax_coords)

    # Base plot (keep mask) as usual
    ax.imshow(
        keep_mask.astype(float),
        extent=extent,
        cmap='viridis',
        interpolation='nearest',
        origin="lower",
        aspect="auto",
        vmin=0, vmax=1
    )

    # Overlay removed pixels in red, transparent elsewhere
    red_overlay = ListedColormap([(0, 0, 0, 0), (1, 0, 0, 1)])  # 0->transparent, 1->red
    ax.imshow(
        removed_mask.astype(int),
        extent=extent,
        cmap=red_overlay,
        interpolation='nearest',
        origin="lower",
        aspect="auto",
        vmin=0, vmax=1,
        alpha=removed_alpha
    )

    ax.set_title(title, pad=3)

    if use_x_label:
        ax.set_xlabel('Window center [bins]')
    if use_y_label:
        ax.set_ylabel('Window size [bins]')

    ax.set_ylim([y_llim, y_ulim])
    ax.set_xlim([x_llim, x_ulim])
    ax.set_yticks(y_ticks)

    return ax