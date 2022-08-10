"""
@author: N. Siomos & P. Paschou

Routines to construct error for analog and photon-counting signals
"""
import numpy as np
import xarray as xr

def error(sig_inp, cfg, info, tframes, isdark, drk_inp = []):

    iters  = cfg.iters
    ierror = cfg.ierror
    
    frames = np.nanmax(tframes.values)
    
    if ierror == 1:
        sig_sim = sig_inp.copy().expand_dims(iters=np.arange(0,iters+1))
        
        info['nsf'] = cfg.nsf

        sig_list = []
        for ind in info.index:
            cursor = dict(channel = ind)

            sig_ch = sig_sim.copy().loc[cursor]
            
            info_ch = info.loc[ind]
            
            an_to_pc = cfg.an_to_pc[ind]

            if info.ch_mode[ind] == 0:
               
                if isdark:
                    drk = sig_inp.copy().loc[cursor][dict(time = 0)]
                    avg_d, std_d = \
                        get_dark(drk, info_ch, frames, cfg.prt_lims.loc[ind,:])
                    noise, qmask = drk_err(sig_ch.copy(), avg_d, std_d)
                
                if not isdark:
                    # if isinstance(drk_inp,list) == False:
                    drk = drk_inp.copy().loc[cursor]
                    avg_d, std_d = \
                        get_dark(drk, info_ch, frames, cfg.prt_lims.loc[ind,:])

                    # else:
                    bgd = sig_inp.copy().loc[cursor]
                    avg_b, std_b = \
                        get_dark(bgd, info_ch, frames, cfg.prt_lims.loc[ind,:])
                    noise, groups, qmask = \
                        an_err(sig_ch.copy(), avg_d, std_d, avg_b, std_b,
                               iters, cfg.nclass, an_to_pc, cfg.an_thres,
                               info_ch, frames)

            if info.ch_mode[ind] == 1:

                noise, groups, qmask = \
                    pc_err(sig_ch.copy(), iters, cfg.nclass, info_ch, frames)
            
            sig_art = add_err(sig_ch.copy(), noise, qmask)
            
            sig_list.append(sig_art)

        # Join arrays and sort according to original channels        
        sig_sim = xr.concat(sig_list, 'channel')
        sig_sim = sig_sim.loc[{'channel': cfg.ch_map}] 
        sig_sim = sig_sim.transpose('iters','time','channel','bins')     
        sig_sim.loc[dict(iters = 0)] = sig_inp.values

    else:
        sig_sim = sig_inp.copy().expand_dims(iters=np.arange(0,1))

        
    return(sig_sim, info)

def add_err(sig, noise, qmask):
    
    sig_coords = sig.coords
    sig_dims = sig.dims
    
    sig = sig.values
    
    for i in range(noise.shape[0]):
        
        err_dist = noise[i,:]
        
        rand_ind = np.random.random_integers(0, high = len(err_dist)-1, 
                                             size = qmask[i].sum())
        
        error = err_dist[rand_ind]
        
        sig[qmask[i]] = sig[qmask[i]] + error 
        
    sig_art = xr.DataArray(data = sig, coords = sig_coords, 
                           dims = sig_dims)
    
    return(sig_art)

def quantize(sig, classes):
   
    sig_arr = np.sort(np.unique(sig))    
        
    # This ensures that the geometrical scale will not find values <= 0
    base = np.abs(sig_arr[0]) + 0.01

    min_sig = sig_arr[0] + base
            
    max_sig = sig_arr[-1] + base
        
    quants = np.geomspace(min_sig, max_sig, num = classes, endpoint = True)
        
    quants = quants - base
        
    quant_dn = quants[:-1]
    
    quant_up = quants[1:]
        
    quant = (quant_up + quant_dn) / 2.
    
    qmask = []
    
    for i in range(quant.size):
        mask = (sig >=  quant_dn[i]) & (sig <  quant_up[i])
        qmask.append(mask)

    return(quant, qmask)

def an_err(sig, avg_d, std_d, avg_b, std_b, 
           iters, classes, an_to_pc, an_thres, info, frames):
    
    sample = 100000
        
    groups, qmask = quantize(sig.copy().values - avg_d, classes = classes)
       
    noise = np.nan * np.zeros((groups.size, sample))

    noise_d = np.random.normal(scale = std_d, size = sample)
    noise_b = np.random.normal(scale = std_b, size = sample)
    
    for i in range(groups.shape[0]):
        
        factor = np.power(info.nsf,2)
        
        sampl_rate = info.sampl_rate
        
        shots = info.shots

        if groups[i] < an_thres:
            noise[i,:] = noise_b 
            
        if groups[i] >= an_thres:
            counts = groups[i] * an_to_pc * shots * frames /\
                (sampl_rate * factor)
            
            if counts <= 0:
                counts = np.min(groups[groups > 0])
            
            noise_s = np.random.poisson(lam = counts, size = sample) - counts
            
            noise_s = noise_s * sampl_rate * factor  /\
                (shots * frames * an_to_pc)
            
            noise[i,:] = noise_s + noise_d
            
    # mid_zone = (groups > 0) & (groups < an_thres)
    
    # if np.sum(mid_zone) > 0:
    #     ind_mid = np.sort(np.where((groups > 0) &\
    #                                (groups < an_thres))[0])

    #     scale = np.linspace(0.,1.,num = len(ind_mid) + 2)

    #     noise_l = noise[ind_mid[0]-1,:]

    #     noise_u = noise[ind_mid[-1]+1,:]

    #     for i in ind_mid:
    #         noise[i,:] = noise_l * (1. - scale[i -ind_mid[0] + 1]) + \
    #             noise_u * scale[i -ind_mid[0] + 1]

    return(noise, groups, qmask)

def drk_err(sig, avg, std):
    
    sample = 100000
               
    qmask = []
    
    qmask.append(sig.values == sig.values)

    noise = np.nan * np.zeros((1, sample))
        
    noise[0,:] = np.random.normal(scale = std, size = sample)
    
    noise[0,:] = noise[0,:]
    
    return(noise, qmask)

def pc_err(sig, iters, classes, info, frames):
    
    sample = 100000
    
    shots = info.shots
    
    sampl_rate = info.sampl_rate
    
    sig_cnt = sig.copy() * shots * frames / sampl_rate
    
    groups, qmask = quantize(sig_cnt.copy().values, classes = classes)
            
    noise = np.nan*np.zeros((groups.size,sample))
    
    for i in range(groups.shape[0]):
        
        counts = groups[i]
        if counts <= 0:
            counts = np.min(groups[groups > 0])

        noise[i,:] = np.random.poisson(lam = counts, size = sample)
        
        noise[i,:] =  (noise[i,:] - counts) * sampl_rate / (shots * frames)
            
    return(noise, groups, qmask)

def get_dark(sig, info, frames, prt_lims):
    
    drk_s = prt_lims.loc['start']
    drk_e = prt_lims.loc['end']
    
    # shots = info.shots
    
    cursor_d = dict(bins = slice(drk_s, drk_e))
        
    std_d = sig.loc[cursor_d].values.std()
    
    # std_d = std_d * np.sqrt(shots * frames)
    
    avg_d = sig.loc[cursor_d].values.mean()

    return(avg_d, std_d)
    
def add_param_err(param, param_err, mc_iter):
    # costruct the parameter's variability according to its statistical error
            
    if np.isscalar(param): # if the param is scalar (not array) then continue
        param = np.array([param])        
        
        if param_err >= 0.:  # variate the param only if given error is non-negative (enabling\disabling option)                                           
            ini_val = param[0]
            param = np.random.normal(param, param_err, mc_iter)
            # keep the given param value in the zero iter
            param[0] = ini_val
            
    return(param)

def extract_error(arr):
# Extract the mean value and the error (std) of an xarray (e.g. signal, dark, prod)

    err = []
    avg = []
    
    if len(arr) > 0:
        avg = arr.copy().loc[dict(iters = 0)].copy()
        iters = arr.iters.size
        if iters > 1:
            err = arr.loc[dict(iters = range(1,iters))].std(dim='iters')
        if iters == 1:
            err = arr.std(dim='iters')
 
    return(avg, err)