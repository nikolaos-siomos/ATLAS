#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 21:31:32 2025

@author: nikos
"""

import os
import numpy as np
from matplotlib import pyplot as plt

def generate_plot(dir_out, fname, title, args, masks):
    
    [X, Y] = np.meshgrid(masks['total'].lower_limit.values, masks['total'].window)
    fig = plt.figure(figsize=(12. , 8.))

    fig.suptitle(title)
    
    fig_x = 0.44
    fig_y = 0.23
    
    fig_edg1_x = 0.06
    fig_edg2_x = 0.54
    
    fig_edg1_y = 0.65
    fig_edg2_y = 0.36
    fig_edg3_y = 0.07
    
    ax1_coords = [fig_edg1_x, fig_edg1_y, fig_x, fig_y]
    ax2_coords = [fig_edg2_x, fig_edg1_y, fig_x, fig_y]
    ax3_coords = [fig_edg1_x, fig_edg2_y, fig_x, fig_y]
    ax4_coords = [fig_edg2_x, fig_edg2_y, fig_x, fig_y]
    ax5_coords = [fig_edg1_x, fig_edg3_y, fig_x, fig_y]
    ax6_coords = [fig_edg2_x, fig_edg3_y, fig_x, fig_y]
        
    plot_single_mask(fig = fig, 
                     ax_coords = ax1_coords,
                     X = Y, 
                     Y = Y, 
                     Z = masks['first_derivative'].values,
                     x_title = 'Derivative mask',
                     use_y_label = True,
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax2_coords,
                     X = Y, 
                     Y = Y, 
                     Z = masks['relative_sem'].values,
                     x_title = 'Relative SEM mask',
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax3_coords,
                     X = Y, 
                     Y = Y, 
                     Z = masks['second_derivative'].values,
                     x_title = 'Second derivative mask',
                     use_y_label = True,
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax4_coords,
                     X = Y, 
                     Y = Y, 
                     Z = masks['shapiro_wilk'].values,
                     x_title = 'Shapiro-Wilk mask',
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax5_coords,
                     X = Y, 
                     Y = Y, 
                     Z = masks['cross_criterion'].values,
                     x_title = 'Cross-check mask',
                     use_x_label = True,
                     use_y_label = True,
                     args = args)
    
    plot_single_mask(fig = fig, 
                     ax_coords = ax6_coords,
                     X = Y, 
                     Y = Y, 
                     Z = masks['total'].values,
                     x_title = 'Combined Mask',
                     use_x_label = True,
                     args = args)
    
    fpath = export_plot(fig, args)
            
    return(fpath)

def plot_single_mask(fig, ax_coords, X, Y, Z, title, args, 
                     use_x_label = False, use_y_label = False):
    
    x_llim = 0.
    x_ulim = args['scanning_region'][1] - args['scanning_window_sizes'][1]
    y_llim = args['scanning_window_sizes'][0]
    y_ulim = args['scanning_window_sizes'][1]
    
    ax = fig.add_axes(ax_coords)
    ax.pcolormesh(X, Y, Z, vmin = 0, vmax = 1)
    
    ax.set_title(title, pad = 5)
    
    if use_x_label:
        ax.set_xlabel('Lower Limit [km]')

    if use_y_label:
        ax.set_ylabel('Window [km]')
    
    ax.set_ylim([y_llim, y_ulim])
    ax.set_xlim([x_llim, x_ulim])
    
    return(ax)

def export_plot(fig, args):
    
    dpi_val = args['dpi']

    dir_out = os.path.join(args['output_folder'],'plots'), 

    fpath = os.path.join(dir_out, f"{args['fname']}.png")
            
    fig.savefig(fpath, dpi = dpi_val)
    
    fig.clf()
    
    plt.close()
    
    return(fpath)