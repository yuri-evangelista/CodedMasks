import numpy as np
from scipy.ndimage import shift as ndshift
from astropy.nddata import block_reduce

def shift(arr, lag):
    shifted = np.roll(arr, lag, axis=1)
    if lag >=0:
        shifted[:, 0:lag] = 0
    else:
        shifted[:, lag:] = 0
    return shifted

def fshift(arr, lagx, lagy):
    """
    Shifts a 2D array (casted to float) with fractional shifts using SciPy.
    
    Note: SciPy's shift convention is (y, x), but our images are already in the right shape
    """
    return ndshift(arr, (lagx, lagy), 'float', order=1, prefilter=True, mode='grid-constant', cval=0.0)

def erosion(arr, cut, step):
    r"""
    Function to erode mask array for vignetting.
    2D matrix erosion for simulating finite thickness effect in shadow projections.
    
    It takes a mask array and "thins" the mask elements across the columns' direction.
    The erosion is performed only on the correct side of open (1) mask elements:
        - right side if cut is negative (= thetaX negative)
        - left side if cut is positive (= thetaX positive)
    The function first erodes all integer bins (replacing 1s with 0s)
    If cut is not integer, then the function applies a fractional transparency to the last eroded bin

              \\       \\  \\
    ___________\\       \\  \\____________
               |\\       \\ |\
    ___________| \\       \\|_\___________
                  \\       \\  \ 
                   \\       \\  \
    ________________\\_______\\__\_________
               <--->           <->
               SHIFT         EROSION   
    
    """
    # number of bins to cut
    ncuts = int(cut / step)
    cutted = arr * (arr & shift(arr, ncuts)) if ncuts else arr

    # array indexes to be fractionally reduced:
    #   - the bin with the decimal values is the one
    #     to the left or right wrt the cutted bins
    erosion_value = abs(cut / step - ncuts)
    border = (cutted - shift(cutted, int(np.sign(cut)))) > 0
    return cutted - border * erosion_value

def ferosion(arr, cut, step):
    """
    2D matrix erosion for simulating finite thickness effect in shadow projections.
    It takes a mask array and "thins" the mask elements across the columns' direction.
    """
    # number of bins to cut
    ncuts = int(cut / step)# + int(np.sign(cut))
    
    arr_mask = (arr > 0) & (fshift(arr, ncuts , 0 ) > 0)
    cutted = arr * arr_mask if ncuts else arr

    # array indexes to be fractionally reduced:
    #   - the bin with the decimal values is the one
    #     to the left or right wrt the cutted bins
    erosion_value =  abs(cut / step - ncuts)
    
    cutted_mask = (
        np.array((cutted > 0), dtype=int) - np.array((fshift(cutted, int(np.sign(cut)), 0) > 0), dtype=int)
    )
    border = (cutted_mask > 0)
    return cutted * (1.0 - border * erosion_value)

def apply_vignetting(detimage, xshift, yshift, focal, ELXDIM, ELYDIM, MTHICK):
    angle_x_rad = -np.arctan(xshift*ELXDIM / focal) #note the -1
    angle_y_rad = -np.arctan(yshift*ELYDIM / focal) #note the -1

    red_factor_x = MTHICK * np.tan(angle_x_rad) 
    red_factor_y = MTHICK * np.tan(angle_y_rad) 

    #The following takes into account the fractional shift and recalculate the red_factor starting from the closer external pixel boundary
    #(see new_erosion_20251024.ipynb for more details)
    
    farest_border_x = abs(xshift - int(xshift))
    farest_border_y = abs(yshift - int(yshift))

    red_factor_x_corr = red_factor_x + np.sign(red_factor_x) * (1- farest_border_x)*ELXDIM
    red_factor_y_corr = red_factor_y + np.sign(red_factor_y) * (1- farest_border_y)*ELYDIM

    sg1 = ferosion(detimage,  red_factor_x_corr, ELXDIM )
    sg2 = ferosion(sg1.T, red_factor_y_corr, ELYDIM)
    return sg2.T

def pad_array(arr, npadx, npady):
    return np.pad(arr, ( (npadx, npadx), (npady, npady)), 'constant', constant_values=(0,0))

def float_gcd(a, b, rtol = 1e-5, atol = 1e-5):
    #Greatest common divisor for floating point numbers (with tolerances)
    t = min(abs(a), abs(b))
    while abs(b) > rtol * t + atol:
        a, b = b, a % b
    return round(a/atol)*atol

def upscale(array, fx, fy):
    for ax, factor in enumerate((fy, fx)):
        array = np.repeat(array, factor, axis=ax)
    return array

def downscale(array, fx, fy, func):
    #requires block_reduce from astropy.nddata
    #func is the downscale function (e.g. np.sum, np.mean...)
    return block_reduce(array, block_size=(fx, fy), func=func)