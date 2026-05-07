"""
@author: A. Gialitaki

Colorbars for contour plotting
"""

import matplotlib.colors as mcolors
import numpy as np

def custom_rgb(rgb, name):
    
    bins = np.linspace(0,1,rgb.shape[0])
    
    color_list = list(zip(bins,rgb))

    cmap = mcolors.LinearSegmentedColormap.from_list(name, color_list)
    
    return(cmap)

def custom_rgb_log(rgb, name):
    
    bins = np.logspace(np.log10(1E-5), np.log10(1), rgb.shape[0] - 1)
    
    bins = np.insert(bins, 0, 0)
    
    color_list = list(zip(bins, rgb))

    cmap = mcolors.LinearSegmentedColormap.from_list(name, color_list)
    
    return(cmap)

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

# auti 
c = mcolors.ColorConverter().to_rgb
cmap = make_colormap([c('lightskyblue'),c('dodgerblue'), 0.1,
                      c('dodgerblue'),c('teal'), 0.2,
                      c('teal'),c('limegreen'), 0.3,
                      c('limegreen'),c('yellow'), 0.4,
                      c('yellow'), c('orange'), 0.5,
                      c('orange'),c('orangered'), 0.6,
                      c('orangered'), c('red'), 0.7,
                      c('red'), c('firebrick'), 0.8,
                      c('firebrick'),c('darkred'), 0.9,
                      c('darkred'),c('k')])
cmap.set_bad(color='white', alpha=1)
    
cmap2 = make_colormap([c('black'), c('indigo'), 0.15,
                  c('indigo'),c('mediumblue'), 0.25,
                  c('mediumblue'),c('blue'), 0.35,
                  c('blue'),c('aqua'), 0.45,
                  c('aqua'),c('lime'), 0.55,
                  c('lime'), c('greenyellow'), 0.65,
                  c('greenyellow'),c('yellow'), 0.75,
                  c('yellow'), c('orangered'), 0.85,
                  c('orangered'), c('red'), 0.95,
                  c('red')])
    
cmap2.set_bad(color='white', alpha=1)
    
cmap3 = make_colormap([c('white'), c('mediumpurple'), 0.1,
                      c('mediumpurple'),c('slateblue'), 0.2,
                      c('slateblue'),c('blue'), 0.3,
                      c('blue'),c('dodgerblue'), 0.4,
                      c('dodgerblue'),c('lawngreen'), 0.5,
                      c('lawngreen'), c('yellow'), 0.6,
                      c('yellow'),c('orange'), 0.7,
                      c('orange'), c('orangered'), 0.8,
                      c('orangered'), c('red'), 0.9,
                      c('red')])
    
cmap3.set_bad(color='white', alpha=1)
    
cmap4 = make_colormap([c('k'), c('darkblue'), 0.15,
                  c('darkblue'),c('royalblue'), 0.25,
                  c('royalblue'),c('dodgerblue'), 0.35,
                  c('dodgerblue'),c('cornflowerblue'), 0.45,
                  c('cornflowerblue'),c('yellowgreen'), 0.55,
                  c('yellowgreen'), c('yellow'), 0.65,
                  c('yellow'),c('orange'), 0.75,
                  c('orange'), c('orangered'), 0.85,
                  c('orangered'), c('darkred'), 0.95,
                  c('darkred')])

cmap4.set_bad(color='white', alpha=1)

