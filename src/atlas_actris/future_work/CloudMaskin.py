import numpy as np

def smooth_signal(signal:np.ndarray, window_len:int):
    """ uniformly smooth the signal """
    return uniform_filter1d(signal, size=window_len, mode='nearest')

def getIdx(vals, arr):
    """ return indeces of vals for a given array arr """
    indices = np.array(((vals - arr[0]) / (arr[1] - arr[0])) + 1, dtype=int)
    return indices

def normalize_RCS(RCS, ATB, idx_region):
    """ Normalize RCS by ATB in idx_region """
    RCS_mean = np.nanmean(RCS[:, idx_region[0]:idx_region[1]])
    ATB_mean = np.nanmean(ATB[idx_region[0]:idx_region[1]])
    RCS_norm = (RCS/RCS_mean * ATB_mean) / ATB
    return RCS_norm

def cloudScreen_MSG(height:np.ndarray, RCS, slope_thres, search_region):
    """
    Cloud screen with Maximum Signal Gradient.

    INPUTS:
        - height (array[height]): Height in meters.
        - RCS (array[time, height]) Range corrected signal. [MHz*m^2]
        - slope_thres (float): Threshold of the slope to determine whether there is strong backscatter signal. [MHz*m]
        - search_region (arraylike[start, stop]):  Cloud search reagion. [m]

    OUTPUTS:
        - flagCloudFree (boolean array[time]): Indicates whether the profile is cloud free.
        - layerStatus (array[height, time]): Layer status for each bin (0: unknown, 1: cloud, 2: aerosol).
    """

    if len(search_region) != 2 or search_region[1] <= height[0]:
        raise ValueError("Not a valid search_region.")

    if search_region[0] < height[0]:
        print(f"Warning: Base of search_region is lower than {height[0]}, setting it to {height[0]}")
        search_region[0] = height[0]

    flagCloudFree = np.zeros(RCS.shape[0], dtype=bool)
    layerStatus = np.zeros_like(RCS, dtype=int)

    # Find indices corresponding to search_region
    search_indx = getIdx(search_region, height)
    
    for indx in range(RCS.shape[0]):
        if np.all(RCS[indx, 0]):
            print(f"Skipping time stamp {indx}")
            continue

        slope = np.concatenate(([0], np.diff(smooth_signal(RCS[indx, :], 10)))) / (height[1] - height[0])
        
        if not np.any(slope[search_indx[0]:search_indx[1]] >= slope_thres):
            flagCloudFree[indx] = True
            
    return flagCloudFree, layerStatus

def cloudscreen_simple(
    RCS_FR:np.ndarray,
    height:np.ndarray,
    heightSearchRegion_FR:list=[500, 14000],
    maxSigSlope4FilterCloud:float=8.0,
    # NR:
    RCS_NR:np.ndarray|None=None,
    heightSearchRegion_NR:list=[150, 5000],
    # Normalization:
    ATB:np.ndarray|None=None,
    normRegion_FR:list|None=None, #[4000, 6000]
    normRegion_NR:list|None=None,
    ):
    """ 
    Cloud Screening / Cloud Mask algorithm from PicassoPy.
    
    INPUTS:
        - RCS_FR (array[time, height]): Range corrected signal.
        - height (array[height]): Height above ground in meters
        - heightSearchRegion_FR (list[low, high]): Cloud search reagion for FR channels in height (meters above ground) the low value should represent the height of a Full Overlap. Default: [500, 7000]
        - maxSigSlope4FilterCloud (float): signal slope treshold, Default: 8.0 (only for normalized signals)
        - RCS_NR (array[time, height]): Range corrected signal. Default: None.
        - heightSearchRegion_FR (list[low, high]): Cloud search reagion for FR channels in height (meters above ground) the low value should represent the height of a Full Overlap. Default: [150, 2000]
        - ATB (array[time, height]): Attenuated Backscatter Deafualt: None
        - normRegion_FR (list[low, high]): Noramlization region for FR. The values should represent the refference height. Defualt: None
        - normRegion_NR (list[low, high]): Noramlization region for NR. The values should represent the refference height. Defualt: None
        
    OUTPUTS:
        - flagCloudFree (array[time]): a 1 dim temporal boolian array. 0 = cloudy, 1 = cloud free. 
    """
    if ATB is not None and normRegion_FR is not None:
        print("Info: Normalize FR signal")
        RCS_FR = normalize_RCS(RCS_FR, ATB, getIdx(normRegion_FR, height))
        
    flagCloudFree, layerStatus = cloudScreen_MSG(
        height=height, 
        RCS=RCS_FR,
        slope_thres=maxSigSlope4FilterCloud, 
        search_region=heightSearchRegion_FR,
        norm_region=normRegion_FR,
    )

    # and for near range if it exists
    if RCS_NR is not None:
        print("Info: Also applied to NR")
        if ATB is not None and normRegion_FR is not None:
            print("Info: Normalize NR signal")
            RCS_NR = normalize_RCS(RCS_NR, ATB, getIdx(normRegion_NR, height))
            
        flagCloudFree_NR, layerStatus_NR = cloudScreen_MSG(
            height=height,
            RCS=RCS_NR,
            slope_thres=maxSigSlope4FilterCloud,
            search_region=heightSearchRegion_NR,
            norm_region=normRegion_NR,
        )

        flagCloudFree = flagCloudFree & flagCloudFree_NR
    return flagCloudFree