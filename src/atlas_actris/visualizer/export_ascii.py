#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:00:07 2022

@author: nick
"""

import numpy as np
import os


def dark(dir_out, fname, header, x_dict, y_dict, y_err_dict):
    
    x_av = x_dict['av'] 
    x_bc = x_dict['bc'] 
    x_sm = x_dict['sm'] 
    x_rc = x_dict['rc'] 
    x_rc_ray = x_dict['rc_ray']
    x_rc_drk_ray = x_dict['rc_drk_ray']
    y_av = y_dict['av'] 
    y_bc = y_dict['bc'] 
    y_sm = y_dict['sm'] 
    y_rc = y_dict['rc'] 
    y_rc_ray = y_dict['rc_ray']
    y_rc_drk_ray = y_dict['rc_drk_ray']
    y_rc_er = y_err_dict['rc'] 
    y_rc_ray_er = y_err_dict['rc_ray']
    y_rc_drk_ray_er = y_err_dict['rc_drk_ray']
    
    body = np.vstack((
        x_av, 
        np.mean(y_av, axis = 0), 
        x_bc, 
        np.mean(y_bc, axis = 0), 
        x_sm,
        np.mean(y_sm, axis = 0),
        x_rc,
        y_rc,
        y_rc_er,
        x_rc_ray,
        y_rc_ray,
        y_rc_ray_er,
        x_rc_drk_ray,
        y_rc_drk_ray,
        y_rc_drk_ray_er,
        )).T

    mask = np.any(np.isnan(body), axis = 1)
    body = body[~mask,:]
    
    fpath = os.path.join(dir_out, 'ascii', fname)
    
    os.makedirs(os.path.join(dir_out, "ascii"), exist_ok=True)

    np.savetxt(fpath, body, header = header, comments = '', 
               delimiter = ',', fmt = '%.6e')
    
    return()

def rayleigh(dir_out, fname, header, alt, atb, rcs):
    
    body = np.vstack((alt, atb, rcs)).T

    mask = np.any(np.isnan(body), axis = 1)
    body = body[~mask,:]
    
    fpath = os.path.join(dir_out, 'ascii', fname)
    
    os.makedirs(os.path.join(dir_out, "ascii"), exist_ok=True)

    np.savetxt(fpath, body, header = header, comments = '', 
               delimiter = ',', fmt = '%.6e')
    
    return()


def telecover(dir_out, fname, header, iters, 
              alt, sectors, sectors_e):

    sectors_l = np.array([sectors[key] for key in sectors_e.keys() 
                          if isinstance(sectors[key], list) == False])
    
    extra_sector_l = np.array([sectors_e[key] for key in sectors_e.keys() 
                               if isinstance(sectors_e[key], list) == False])

    if (extra_sector_l).shape[0] > 0:
        body = np.concatenate((np.array([alt]), sectors_l, extra_sector_l), axis = 0).T
    else:
        body = np.concatenate((np.array([alt]), sectors_l), axis = 0).T
    
    mask = np.any(np.isnan(body), axis = 1)
    body = body[~mask,:]
 
    dir_ascii = os.path.join(dir_out, "ascii")
    os.makedirs(dir_ascii, exist_ok=True)

    fpath = os.path.join(dir_ascii, fname)

    np.savetxt(
        fpath,
        body,
        header=header,
        comments="",
        delimiter=",",
        fmt="%.6e",
    )
    
    return()

def polarisation_calibration(dir_out, fname, alt_cal, alt_ray,
                             r_p45, t_p45, r_m45, t_m45, ray_r, ray_t, header):
    
    if (np.abs(alt_cal - alt_ray) < 1E-3).all():
        alt = alt_cal
        
    body = np.vstack((alt, r_p45, t_p45, r_m45, t_m45, ray_r, ray_t)).T

    mask = np.any(np.isnan(body), axis = 1)
    body = body[~mask,:]
    
    fpath = os.path.join(dir_out, 'ascii', fname)
    
    np.savetxt(fpath, body, header = header, comments = '', 
               delimiter = ',', fmt = '%.6e')
    
    return()