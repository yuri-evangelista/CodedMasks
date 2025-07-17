import numpy as np
from scipy.signal import convolve
from scipy.signal import correlate
from mask_utils.image_utils import *


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


def snr_vs_off_axis(s_counts, b_counts, mask, bulk, mask_x_pitch, mask_y_pitch, 
                    ELXDIM, ELYDIM, det_x_pitch, focal, mask_thickness, 
                    thetaX, thetaY, degrees=True, verbose=False):
    # Returns the approximate sensitivity following eq. 13 and eq. 23 of Skinner 2008
    
    if degrees:
        thetaX = np.deg2rad(thetaX)
        thetaY = np.deg2rad(thetaY)
    

    on_axis_of = np.sum(mask)/mask.size
    coding_power = get_coding_power(mask_x_pitch, det_x_pitch, on_axis_of)
    off_axis_of = open_fraction_vs_off_axis(mask, mask_thickness, ELXDIM, ELYDIM, thetaX, thetaY, degrees=False)
    
    off_axis_area =  0.01 * eff_area_vs_off_axis(mask.T.astype('int32'), bulk.T, ELXDIM, ELYDIM, focal, mask_thickness, thetaX, thetaY, degrees=False)
    

    if verbose:
        print("On-axis OF:", on_axis_of)
        print("Coding power:", coding_power)
        print("Off-axis OF:", off_axis_of)
        print("Off-axis area:", off_axis_area)

    snr = coding_power * s_counts * np.sqrt( (off_axis_area * (1-off_axis_of))/ (off_axis_of * s_counts + b_counts))


    return snr