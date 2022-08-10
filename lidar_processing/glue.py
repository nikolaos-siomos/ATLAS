"""
@authors: N. Siomos & P. Paschou

Gluing of signals
"""
import numpy as np
from scipy.stats import linregress, pearsonr

def gluing(signal, info, gl, bc):
    
    glue_map = gl.map
    hwin = gl.hwin_len
    ref_sig = gl.ref_sig
        
    if len(glue_map)>0:
        for j in glue_map.index:
            
            gl_ind_list = []
            
            # Extracting glue channels
            ch_l  = glue_map.id_lower[j]
            ch_u  = glue_map.id_upper[j] 
            ch_gl = glue_map.id_glue[j]
            
            fail = (0.*signal[0,:,0,0].copy()).astype(bool)
            
            # Add new index in the dimension of channel for the glued signals
            signal = signal.reindex(channel = list(signal.channel.values) + [ch_gl])    
            
            # Add the info about the glued signal in the info_array 
            # It is similar to the photon channel, only type is different
            info.loc[ch_gl] = info.loc[ch_u].copy()
            info.loc[ch_gl, 'ch_mode'] = 3.
            info.loc[ch_gl, 'full_ovl_idx'] = info.loc[ch_l, 'full_ovl_idx']
            
            # Declare channels processed
            print(f'-> Glue: Adding channel {ch_gl}')
            
            # Convert float variable window to even half bin window
            ihwin = np.floor(hwin/(info.resol.loc[ch_gl])).astype(int)

            for t in range(signal.time.size):
                ind_t = dict(time = t)
                # Some shortcuts
                sig_l  = signal[ind_t].loc[dict(channel = ch_l)].values
                sig_u  = signal[ind_t].loc[dict(channel = ch_u)].values            
                bc_u = bc[ind_t].loc[dict(channel = ch_u)].values
                range_ar = signal.range.values
                
                # Searching for the best gluing point
                gl_ind = glue_region(an = sig_l.mean(axis=0), 
                                     pc = sig_u.mean(axis=0), 
                                     bc_pc = bc_u.mean(axis=0),
                                     range_ar = range_ar,
                                     ihwin = ihwin,
                                     metrics = gl)
                
                # Calculate glue signal
                if gl_ind == gl_ind:
                    for n in signal.iters.values:
                        sig_gl = glue_calc(sig_l[n,:], sig_u[n,:], 
                                           gl_ind, ihwin, ref_sig)          
                        signal[ind_t].loc[dict(channel = ch_gl, iters = n)] =\
                            sig_gl                                
                                
                else:
                # no gluing point found
                    
                    if ref_sig == 'upper':
                    # Switch to the upper signal (chosed as reference)
                        signal[ind_t].loc[dict(channel = ch_gl)] = sig_u

                    if ref_sig == 'lower':
                    # Switch to the lower signal (chosed as reference)
                        signal[ind_t].loc[dict(channel = ch_gl)] = sig_l

                    fail[ind_t] = True
                
                # Save the gluing point to report it
                gl_ind_list.append(gl_ind)

                if fail[ind_t]:
                    timefail =  signal.time.values[t]
                    print(f'        {str(timefail)[0:19]}: '+\
                          f'Cannot glue -->  switching to {ref_sig} signal!')  
            
            # Report the approximate gluing point of all time frames 
            gl_range = glue_range_check(signal, gl_ind_list)
            # save the approximate gluing point of all time frames in info dataframe
            info.loc[ch_gl, 'glue_rng'] = gl_range

    return(signal, info)

def glue_range_check(signal, ind_list):
    
    ind = np.nanmean(ind_list)

    if ind == ind:
        gl_range = signal.range.values[int(ind)]
        print(f'           Gluing at {gl_range} m range')
    else:
        print('           Could not glue for any time frame!')
        gl_range = np.nan
    return(gl_range)

def glue_map_check(info, ch_l, ch_u):
    
    if info.loc[ch_l, 'ch_mode'] != 0.:
        raise TypeError('Lower glue channel type is not analog! Program stopped \n'+\
                        'Please revise the gluing_map.ini')
        
    if info.loc[ch_u, 'ch_mode'] != 1.: 
        raise TypeError('Upper glue channel type is not an photon counting! Program stopped \n'+\
                        'Please revise the gluing_map.ini')
   
    return()     

def glue_calc(sig_l, sig_u, gl_ind, ihwin, ref_sig):

    sig_gl = np.nan*np.copy(sig_u)

    ind_l = gl_ind - ihwin
    ind_u = gl_ind + ihwin + 1
    
    # Weights to fade in - fade out in the glue region
    weights = np.linspace(1., 0., sig_l[ind_l:ind_u].size)
    
    if ref_sig == 'upper':
        # Slope with 0 intercept is the conversion (normalization) factor 
        # between the upper and lower signals -> with respect to upper signal
        slope = np.mean(sig_u[ind_l:ind_u]/sig_l[ind_l:ind_u])
        
        # Perform gluing with the upper signal as reference
        sig_gl[:ind_l] = slope*sig_l[:ind_l]
        sig_gl[ind_u:] = sig_u[ind_u:]
    
        sig_gl[ind_l:ind_u] = slope*sig_l[ind_l:ind_u]*weights + \
            sig_u[ind_l:ind_u]*(1. - weights) 

    if ref_sig == 'lower':
        # Slope with 0 intercept is the conversion (normalization) factor 
        # between the  upper and lower signals -> with respect to lower signal
        slope = np.mean(sig_l[ind_l:ind_u]/sig_u[ind_l:ind_u])

        # Perform gluing with the lower signal as reference
        sig_gl[:ind_l] = sig_l[:ind_l]
        sig_gl[ind_u:] = slope*sig_u[ind_u:]
    
        sig_gl[ind_l:ind_u] = sig_l[ind_l:ind_u]*weights + \
            slope*sig_u[ind_l:ind_u]*(1. - weights) 
    # print(f'ln 149  norm f = {slope}')
    
    return(sig_gl)#,slope)

def glue_region(an, pc, bc_pc, range_ar, ihwin, metrics):
    # search from bottom to top to find glue region
    
    ind = np.arange(0, an.size, 1)

    slope = np.nan*np.zeros(ind.size) # linear fit slope
    r = np.nan*np.zeros(ind.size) # Pearson correlation
    p = np.nan*np.zeros(ind.size) # p_value for regression

    # Get rid of regions that are impossible to glue anyway (e.g. photon saturation > 60.MHz, too low analog < 0.5mV)

    mask_gl = np.ones(an.size, dtype=bool)    

    if metrics.l_sig_min != '' and metrics.u_sig_max != '' and metrics.l_range_min == '' and metrics.u_range_max == '':
        mask_gl = (an > metrics.l_sig_min) & \
            (pc + bc_pc < metrics.u_sig_max) & \
            (an == an) & (pc == pc) 
    
    if metrics.l_sig_min == '' and metrics.u_sig_max == '' and metrics.l_range_min != '' and metrics.u_range_max != '':
        mask_gl = (range_ar > metrics.l_range_min) & \
            (range_ar < metrics.u_range_max) & \
            (an == an) & (pc == pc)
    

    if metrics.l_sig_min != '' and metrics.u_sig_max != '' and metrics.l_range_min != '' and metrics.u_range_max != '':
        mask_gl = (an > metrics.l_sig_min) & \
            (pc + bc_pc < metrics.u_sig_max) & \
            (range_ar > metrics.l_range_min) & \
            (range_ar < metrics.u_range_max) & \
            (an == an) & (pc == pc)
        
    
    ind = ind[mask_gl] 

    if len(ind) > 2.*(ihwin + 1):
        # # Normalize to get more similar signals for the correlation calculation
        # an = an/np.mean(an[mask_gl]) 
        # pc = pc/np.mean(pc[mask_gl])
        
        # Avoid gluing closer to half window towards the start and end of the profile
        # (although highly unlikely)
        # ind = ind[(ind - ihwin > 0.) & (ind[-1] - (ind + ihwin) > 0.)]
        ind = ind[ihwin:-ihwin]
            
        # Calculate correlation and regression
        for i in ind:
            s = slice(i - ihwin, i + ihwin + 1)
            #cor_coef
            r[i], _ = pearsonr(an[s], pc[s])
            # slope of the fiting in residual an - pc signals
            slope[i], _, _, p[i], _ = linregress(np.arange(i-ihwin, i+ihwin+1),
                                                 y = an[s] - pc[s])     

        # Gluing point corresponds to maximum correlation when the regression is high    
        mask_lin = (r > metrics.cor_coef)

        
        if len(r[mask_lin]) > 0:
        # Gluing point also corresponds to minimum slope of the regression between the residual: an - pc
            opt = np.nanmin(slope[mask_lin])
            gl_ind = np.where((slope == opt) & (mask_lin == True))[0][0]
            
        else:
            gl_ind = np.nan
    
    else:
        gl_ind = np.nan
        
    return(gl_ind)