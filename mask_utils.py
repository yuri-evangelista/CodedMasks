import numpy as np
from sympy import *
import numpy as np
import astropy.io.fits as pyfits 
from astropy.table import Table, Column
from scipy.signal import convolve
from scipy.signal import correlate
from scipy.ndimage import shift as ndshift

def ura_mura(p):
    #Checking if p is prime
    if isprime(p):
        #Preparing arrays
        A = np.zeros(p, dtype=int)
        R = np.arange(p, dtype=np.int64)**2 % p
        A[R] = 1

        #Check if URA or MURA can be generated
        URA = ( round( (p-3)/4 ) - ((p-3)/4)) == 0
        MURA = ( round( (p-1)/4 ) - ((p-1)/4)) == 0

        if URA:
            print("Generating URA array")
        elif MURA:
            print("Generating MURA array")
            A[0] = 0
        return A
    else:
        raise TypeError("p must be prime")

def bura(p, modified=False):
    #Checking if p is prime
    if not isprime(p):
        raise TypeError("Number of array elements must be prime")

    #Check if BURA can be generated
    x = np.sqrt( (p-1)/4)
    if not ( (int(x) == x) & ( (x % 2) == 1) ) :
        raise TypeError("p does not fulfill the requirement of p = 4x^2 + 1 with x odd")
        
    #Preparing arrays
    A = np.zeros(p, dtype=np.int64)
    R = np.arange(p, dtype=np.int64)**4 % p
    A[R] = 1

    if modified:
        A[0] = 0
        print("Generating M-BURA array")
    else:
        print("Generating BURA array")
    
    return A

def bura33(p):
    """
    Generates a biquadratic URA with OF ~0.33.
   
    ********************************************************
    The resulting code does not seem "perfect"!!!!!
    ********************************************************
    
    From Baumert L. D. 1971, Lecture Notes in Mathematics No. 182, Cyclic Difference Sets
    Theorem 5.18 (iii)

    ----------------------------------------------------------------------------------------------
    The biquadratic residues of primes v = 4x^2 + 9, x odd, form a difference set with ~0.33 OF
    Example primes: 13, 109, 1453, 3373, 3853, 4909, 6733

    """
    #Checking if p is prime
    if not isprime(p):
        raise TypeError("Number of array elements must be prime")

    x = np.sqrt((p - 9)/4.0)
    if not ( (int(x) == x) & ( (x % 2) == 1) ) :
        raise TypeError("p does not fulfill the requirement of p = 4x^2 + 9 with x odd")
 
    #Preparing arrays
    A = np.zeros(p, dtype=np.int64)
    R = np.arange(p, dtype=np.int64)**4 % p
    A[R] = 1
    #A[0] = 0 #not sure about that

    return A

def shift(arr, lag):
    shifted = np.roll(arr, lag, axis=1)
    if lag >=0:
        shifted[:, 0:lag] = 0
    else:
        shifted[:, lag:] = 0
    return shifted

def erosion(arr, cut, step):
    # number of bins to cut
    ncuts = int(cut / step)
    cutted = arr * (arr & shift(arr, ncuts)) if ncuts else arr

    # array indexes to be fractionally reduced:
    #   - the bin with the decimal values is the one
    #     to the left or right wrt the cutted bins
    erosion_value = abs(cut / step - ncuts)
    border = (cutted - shift(cutted, int(np.sign(cut)))) > 0
    return cutted - border * erosion_value

def erosion_old(arr, cut, step):
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
    ncuts = int(cut / step)  #find the integer numbers of bins to cut
    decimal = abs(cut/step - ncuts) #find the decimal number of bins to cut
    
    #Find arr indexes to be cut (completely)
    shifted = shift(arr, ncuts)
    if ncuts:
        eroded_int = arr * ( (arr-shifted) > 0)
    else:
        eroded_int = arr * 0
    
    #Finds arr indexes to be fractionally reduced
    if decimal:
        shifted_p1 = shift(shifted, int(np.sign(cut)))
        eroded_frac = arr * (( (shifted_p1 - arr) < 0) & ~eroded_int)
    else:
        eroded_frac = arr * 0


    #Calculate output array
    out = (arr * (eroded_int < 1)) - eroded_frac * decimal

    return out

def next_prime(n):
    p=n
    while (not isprime(p)):
        p+=1
    return(p)

def get_openfraction(mask):
	return np.sum(mask)/mask.size

def get_angular_res(m_pitch, d_pitch, m_d_distance, degrees=False):
	res = np.sqrt( np.arctan(m_pitch/m_d_distance)**2 + np.arctan(d_pitch/m_d_distance)**2)
	if degrees:
		res = np.rad2deg(res)
	return res

def get_coding_power(m_pitch, d_pitch, open_fraction):
	#from Skinner 2008
	delta = (1 - d_pitch/(3*m_pitch) * np.sqrt(4*open_fraction*(1-open_fraction)))
	return delta

def open_fraction_vs_off_axis(mask, mask_thickness, mask_x_pitch, mask_y_pitch, thetaX, thetaY, degrees=True):
	if degrees:
		thetaX = np.deg2rad(thetaX)
		thetaY = np.deg2rad(thetaY)
	cutx = mask_thickness * np.tan(thetaX) / mask_x_pitch
	cuty = mask_thickness * np.tan(thetaY) / mask_y_pitch
	#print(cutx, cuty)
	vignetted_x = erosion(mask, cutx, 1)
	vignetted_y = erosion(mask.T, cuty, 1)
	#print(get_openfraction(vignetted_x))
	#print(get_openfraction(vignetted_y))
	return get_openfraction(vignetted_x * vignetted_y.T)

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

def eff_area_vs_off_axis(mask, det, x_pitch_ups, y_pitch_ups, focal, mask_thickness, thetaX, thetaY, degrees=True):
    # Returns the area in same units of the inputs (mm2 if pitch and focal are mm, cm2 if pitch and focal are cm)
    # It does (!) take into account of vignetting 
    if degrees:
        thetaX = np.deg2rad(thetaX)
        thetaY = np.deg2rad(thetaY)

    #Apply vignetting
    cutx = mask_thickness * np.tan(thetaX) / x_pitch_ups
    cuty = mask_thickness * np.tan(thetaY) / y_pitch_ups
    vignetted_x = erosion(mask, cutx, 1)
    vignetted_y = erosion(mask.T, cuty, 1)
    vignetted = vignetted_x * vignetted_y.T

    x_shift_px, y_shift_px = int(focal * np.tan(thetaX)/x_pitch_ups), int(focal * np.tan(thetaY)/y_pitch_ups)
    #print(x_shift_px, y_shift_px)
    m = shift(vignetted, x_shift_px) #CHECK
    m = shift(m.T, y_shift_px) #CHECK
    m = m.T
    #print(get_openfraction(m))
    eff_area = np.sum( m * (  det > 0) ) * x_pitch_ups * y_pitch_ups

    return eff_area



def omega_plate_offaxis(a, b, d, A, B):
    r"""
    https://vixra.org/pdf/2001.0603v2.pdf eq.34
    This is when the line-of-sight hits the rectangular plate off center
    - d is the distance between the plate and the "observer"
    - a, b are the dimensions of the rectangular plate
    - A, B are the corner distance from the observer (0<= A <=a/2, 0<= B <=b/2)
    That means we have 4 rectangles (with the observer located on one corner of each one):
        r1 = A * (b-B)
        r2 = (a-A) * (b-B)
        r3 = A * B
        r4 = (a-A) * B

    <-----------a----------->

     _______|_______________    
    |       |               |   |
    |   1   |       2       |   |
    |       |               |   |
    |___A___|_______________|   |b
    |       |               |   |
    |   3  B|       4       |   |    
    |_______|_______________|   |
            |
    
    """

    def omega_plate(a, b, d):
        r"""
        https://vixra.org/pdf/2001.0603v2.pdf eq.27
        a, b are the dimensions of the rectangular plate
        d is the distance of the "observer" from the plate
        """
        alpha = a/(2*d)
        beta  = b/(2*d)
        return 4 * np.arctan( (alpha * beta)/ np.sqrt(1 + alpha**2 + beta**2) )

    omega1 = omega_plate(2*A,2*(b-B), d)
    omega2 = omega_plate( 2*(a-A), 2*(b-B), d)
    omega3 = omega_plate(2*A, 2*B, d)
    omega4 = omega_plate(2*(a-A), 2*B, d)

    return (omega1 + omega2 + omega3 + omega4)/4



def solid_angle(bulk, xstep, ystep, m_d_distance, nobulk=False):
    """
    bulk_shape = shape(bulk)
    xsize, ysize = (bulk_shape[0] * xstep, bulk_shape[1] * ystep)
    xcoords = np.arange(0, bulk_shape[0]) * xstep + xstep/2 - xsize/2
    ycoords = np.arange(0, bulk_shape[1]) * ystep + ystep/2 - ysize/2

    omega = np.zeros(bulk_shape)

    for ix in range(bulk_shape[0]):
        for iy in range(bulk_shape[1]):
            A = xsize/2 - abs(xcoords[ix])
            B = ysize/2 - abs(ycoords[iy])
            omega[ix, iy] = omega_plate_offaxis(xsize, ysize, m_d_distance, A, B) * bulk[ix, iy]

    return omega
    """
    bulk = np.asarray(bulk)
    xsize, ysize = bulk.shape[0] * xstep, bulk.shape[1] * ystep

    #Generates arrays of pixel centers
    xcoords = (np.arange(bulk.shape[0]) + 0.5) * xstep - xsize/2
    ycoords = (np.arange(bulk.shape[1]) + 0.5) * ystep - ysize/2

    #Calculates A, B for each pixel
    A = xsize / 2 - np.abs(xcoords[:, np.newaxis]) #shape: (Nx, 1)
    B = ysize / 2 - np.abs(ycoords[np.newaxis, :]) #shape: (1, Ny)

    if not nobulk:
        #Calculates omega for each pixel using bulk as a bitmask
        omega_values = omega_plate_offaxis(xsize, ysize, m_d_distance, A, B) * (bulk > 0)
    else:
        omega_values = omega_plate_offaxis(xsize, ysize, m_d_distance, A, B)
    return omega_values

def read_fits_events(filein, header0=False, header1=False, verbose=False):
    with pyfits.open(filein) as hdu_list:
        header_0 = hdu_list[0].header
        header_1 = hdu_list[1].header
        if verbose:
            hdu_list.info()
        data = hdu_list[1].data
        events = Table(data)

        if not header0 and not header1:
            return events
        elif header0 and not header1:
            return events, header_0
        elif not header0 and header1:
            return events, header_1
        elif header0 and header1:
            return events, header_0, header_1

def filter_source(data, ra, dec, verbose=False):
    mask = np.ones(len(data), dtype=bool)
    mask &= (np.isclose(data["RA"], ra) & np.isclose(data["DEC"], dec))
    filtered = data[mask]
    if verbose:
        print("Selected", len(filtered), "out of", len(data), "events")

    return filtered

def read_mask_bulk(fitsfile, ext, header_out=False, verbose=False):
    r"""
    Reads data from wfm_mask.fits
    Extensions are:
        0  PRIMARY       1 PrimaryHDU      28   ()      
        1  OR_MASK       1 BinTableHDU     36   676000R x 3C   [E, E, E]   
        2  MASK          1 BinTableHDU     36   676000R x 3C   [E, E, E]   
        3  RMATRIX       1 BinTableHDU     38   676000R x 3C   [E, E, E]   
        4  SENS          1 BinTableHDU     36   676000R x 3C   [E, E, E]   
    """

    with pyfits.open(fitsfile) as hdu_list:     
        if verbose:
            hdu_list.info()
        header = hdu_list[ext].header
        NELE   = header['NAXIS2']
        ELXDIM = header['ELXDIM']
        ELYDIM = header['ELYDIM']
        ELXN   = header['ELXN']
        ELYN   = header['ELYN']
        MINX   = header['MINX']
        MINY   = header['MINY']
    
    
        data=Table(hdu_list[ext].data)
        arr = np.zeros((ELXN, ELYN))
        for ele in range(NELE):
            ix = int(np.round(  ((data['X'][ele] - MINX - ELXDIM/2)/ELXDIM) ) )
            iy = int(np.round(  ((data['Y'][ele] - MINY - ELYDIM/2)/ELYDIM) ) )
            arr[ix,iy] = data['VAL'][ele]

        if header_out:
            return arr, header
        else:
            return arr

def get_detimage(data, xedges, yedges):
    detimg, _, _ = np.histogram2d( data["X"], data["Y"], bins=(xedges, yedges) )
    return detimg

def decode(detimage, rmatrix, bulk):
    cc = correlate(rmatrix, detimage, mode="full")
    balancing = correlate(rmatrix, bulk, mode="full")
    cc_bal = cc - balancing * np.sum(detimage) / np.sum(bulk)
    return cc_bal

def decode_var(detimage, rmatrix, bulk, m_d_distance, elxdim, elydim):
    #First of all we calculate the total detector counts and total active elelements
    sum_det, sum_bulk = map(np.sum, (detimage, bulk))


    #Then we "estimate" the background shape on the detector image considering only the system geometry 
    # (i.e. the solid angle seen by each active pixel)
    omega = solid_angle(bulk, elxdim, elydim, m_d_distance)

    #and we normalize it to have a sum(omega) = 1
    omega_norm = omega/np.sum(omega)

    #Now we multiplicate it for the total detector counts in order to obtain the array of the "expected values"
    Lambda = omega_norm * sum_det

    #We define
    xi = correlate(rmatrix, bulk, mode="full")
    var = correlate(np.square(rmatrix), Lambda, mode="full")
    cov = correlate(rmatrix, Lambda, mode="full")

    #Then (https://www.overleaf.com/project/6863d6e3a6ffe4d95fe704fb)
    var_bal =  var + sum_det/sum_bulk**2 * xi**2 - 2/sum_bulk * xi * cov
 
    return var_bal

def get_skysign(skyimage, varimage):
    var_clipped =  np.clip(varimage, a_min=1E-8, a_max=1E8) if np.any(varimage <= 0) else varimage
    return skyimage/np.sqrt(var_clipped)

def get_detimage_edges(xstep, ystep, nx, ny):
    xedges = np.arange(0, xstep * nx, xstep) - (xstep * nx)/2
    yedges = np.arange(0, ystep * ny, ystep) - (ystep * ny)/2

    xedges = np.append(xedges, xedges[-1] + xstep)
    yedges = np.append(yedges, yedges[-1] + ystep)

    return xedges, yedges

def get_skycoords(skyimage, xstep, ystep, m_d_distance, verbose=False, radians=False):
    s = shape(skyimage)

    x = (np.linspace(0, s[0]*xstep, num=s[0]) )
    y = (np.linspace(0, s[1]*ystep, num=s[1]) )
  
    x -= x[-1]/2
    y -= y[-1]/2

    if verbose:
        print("X bins range:", np.min(x), np.max(x))
        print("Y bins range:",np.min(y), np.max(y))

    if radians:
        return np.arctan(x/m_d_distance), np.arctan(y/m_d_distance)
    else:
        return np.rad2deg(np.arctan(x/m_d_distance)), np.rad2deg(np.arctan(y/m_d_distance))


def fshift(arr, lagx, lagy):
    """
    Shifts a 2D array (casted to float) with fractional shifts using SciPy.
    
    Note: SciPy's shift convention is (y, x), but our images are already in the right shape
    """
    return ndshift(arr, (lagx, lagy), 'float', order=1, prefilter=True, mode='grid-constant', cval=0.0)

def generate_bulk(mask_shape, ELXDIM, ELYDIM):
    #From WFM detector plane geometry
    det_ext_border_x = 157.996
    det_ext_border_y = 153.176
    det_int_border_x = 28.204
    det_int_border_y = 12.824

    #Defining bulk array
    bulk = np.ones(mask_shape)

    #Removing non-sensitive regions along the X direction
    bulk[0: int(round((bulk.shape[0] -  det_ext_border_x/ELXDIM)/2)),      :] = 0
    bulk[-int(round((bulk.shape[0] -  det_ext_border_x/ELXDIM)/2)) : ,      :] = 0
    bulk[int(round((bulk.shape[0] -  det_int_border_x/ELXDIM)/2)) - 1 : -int(round((bulk.shape[0] -  det_int_border_x/ELXDIM)/2) - 1)  ,      :] = 0
    #Removing non-sensitive regions along the Y direction
    bulk[:,      0: int(round((bulk.shape[1] -  det_ext_border_y/ELYDIM)/2))] = 0
    bulk[:,      -int(round((bulk.shape[1] -  det_ext_border_y/ELYDIM)/2)) : ] = 0
    bulk[:,      int(round((bulk.shape[1] -  det_int_border_y/ELYDIM)/2)) : -int(round((bulk.shape[1] -  det_int_border_y/ELYDIM)/2))  ] = 0

    return bulk