import numpy as np
from scipy.ndimage import shift as ndshift

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