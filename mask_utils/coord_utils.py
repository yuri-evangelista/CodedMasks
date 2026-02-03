import numpy as np

def rotation_matrices(pnt_ra_z, pnt_dec_z, pnt_ra_x, pnt_dec_x):
    theta_z = np.deg2rad(90 - pnt_dec_z)
    phi_z = np.deg2rad(pnt_ra_z)

    theta_x = np.deg2rad(90 - pnt_dec_x)
    phi_x = np.deg2rad(pnt_ra_x)

    sin_theta_x = np.sin(theta_x)
    x_axis = np.array([sin_theta_x * np.cos(phi_x), sin_theta_x * np.sin(phi_x), np.cos(theta_x)])

    sin_theta_z = np.sin(theta_z)
    z_axis = np.array([sin_theta_z * np.cos(phi_z), sin_theta_z * np.sin(phi_z), np.cos(theta_z)])

    y_axis = np.array(
        [
            z_axis[1] * x_axis[2] - z_axis[2] * x_axis[1],
            z_axis[2] * x_axis[0] - z_axis[0] * x_axis[2],
            z_axis[0] * x_axis[1] - z_axis[1] * x_axis[0],
        ]
    )

    rotmat_sky2cam = np.vstack((x_axis, y_axis, z_axis))
    rotmat_cam2sky = rotmat_sky2cam.T

    return rotmat_sky2cam, rotmat_cam2sky

def equatorial2camera(src_ra, src_dec, pnt_ra_z, pnt_dec_z, pnt_ra_x, pnt_dec_x, md_distance, angle=False):
    rotmat_sky2cam, _ = rotation_matrices(pnt_ra_z, pnt_dec_z, pnt_ra_x, pnt_dec_x)
    ra = np.deg2rad(src_ra)
    dec = np.deg2rad(src_dec)
    w = np.array(
        [
            np.cos(ra) * np.cos(dec),
            np.sin(ra) * np.cos(dec),
            np.sin(dec),
        ]
    )
    vx, vy, vz = np.matmul(rotmat_sky2cam, w)
    # the sky-shifts are computed from the versor `v` using the mask-detector distance
    shift_x = vx * md_distance / vz
    shift_y = vy * md_distance / vz
    theta_x = np.rad2deg(np.atan(shift_x/md_distance))
    theta_y = np.rad2deg(np.atan(shift_y/md_distance))
    if angle:
        return float(theta_x), float(theta_y)
    else:
        return float(shift_x), float(shift_y)

def camera2equatorial(src_x, src_y, pnt_ra_z, pnt_dec_z, pnt_ra_x, pnt_dec_x, md_distance, angle=False):

    if angle:
        shift_x = md_distance * np.tan(np.deg2rad(src_x))
        shift_y = md_distance * np.tan(np.deg2rad(src_y))
    else:
        shift_x = src_x
        shift_y = src_y
    _, rotmat_cam2sky = rotation_matrices(pnt_ra_z, pnt_dec_z, pnt_ra_x, pnt_dec_x)
    r = np.sqrt(shift_x**2 + shift_y**2 + md_distance**2)
    v = np.array([shift_x, shift_y, md_distance]) / r
    wx, wy, wz = np.matmul(rotmat_cam2sky, v)
    # the versors above are in the rectangular coordinates, we transform into angles
    dec = 0.5 * np.pi - np.arccos(wz)
    ra = np.arctan2(wy, wx)
    ra += 2 * np.pi if ra < 0 else 0.0
    dec = np.rad2deg(dec)
    ra = np.rad2deg(ra)
    return float(ra), float(dec)

'''
Examples:

CAMZRA  = 266.4 # Pointing RA of camera - Z-axis [deg]
CAMZDEC = -28.94 # Pointing DEC of camera - Z-axis [deg]
CAMXRA  = 266.4 # Pointing RA of camera - X-axis [deg]
CAMXDEC = 61.06 # Pointing DEC of camera - X-axis [deg]

src_ra, src_dec = 244.979705810546, -15.6401 #Sco X-1

#Source shifts in mm
equatorial2camera(src_ra, src_dec, CAMZRA, CAMZDEC, CAMXRA, CAMXDEC, 202.9)

#Source ThetaX, ThetaY angles in degrees
equatorial2camera(src_ra, src_dec, CAMZRA, CAMZDEC, CAMXRA, CAMXDEC, 202.9, angle=True)

#Source RA, Dec from shifts in mm
camera2equatorial(43.87696328021216, 77.98815804985986, CAMZRA, CAMZDEC, CAMXRA, CAMXDEC, 202.9)

#Source RA, Dec from ThetaX, ThetaY angles in degrees
camera2equatorial(12.20227430579432, 21.025135611459373, CAMZRA, CAMZDEC, CAMXRA, CAMXDEC, 202.9, angle=True)
'''


'''
The following should work to calculate camera pointings (CAMZRA, CAMZDEC) and (CAMXRA, CAMXDEC) for a distribution of cameras
'''
import numpy as np

def spherical_to_cartesian(ra_deg, dec_deg):
    """Convert RA/Dec to a unit cartesian vector."""
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])

def cartesian_to_spherical(vec):
    """Convert cartesian vector back to RA/Dec. Handles shape (3, N)."""
    vec = np.array(vec)
    if vec.ndim > 1:
        norm = np.linalg.norm(vec, axis=0)
    else:
        norm = np.linalg.norm(vec)
        
    x, y, z = vec / norm
    z = np.clip(z, -1.0, 1.0) # Avoid numerical errors > 1.0
    
    dec = np.arcsin(z)
    ra = np.arctan2(y, x)
    
    return np.degrees(ra) % 360, np.degrees(dec)

def get_camera_pointing(ra_inst, dec_inst, roll_deg, elevation_deg, azimuths_deg):
    """
    Calculates pointing of cameras (Vectorized).
    Returns dictionary with RA/Dec and the input Az/El for reference.
    """
    
    # --- 1. PREPARE INPUTS & BROADCASTING ---
    az_arr = np.atleast_1d(azimuths_deg)
    el_arr = np.atleast_1d(elevation_deg)

    # We use numpy broadcasting to ensure az and el arrays match in size
    # This creates two arrays of the same shape (N,), repeating values if necessary.
    az_broadcast, el_broadcast = np.broadcast_arrays(az_arr, el_arr)

    # --- 2. ESTABLISH CELESTIAL BASIS ---
    k_cel = spherical_to_cartesian(ra_inst, dec_inst)
    
    ra_rad = np.radians(ra_inst)
    dec_rad = np.radians(dec_inst)
    j_cel = np.array([
        -np.sin(dec_rad) * np.cos(ra_rad),
        -np.sin(dec_rad) * np.sin(ra_rad),
        np.cos(dec_rad)
    ])
    
    i_cel = np.array([
        -np.sin(ra_rad),
        np.cos(ra_rad),
        0.0
    ])

    # --- 3. APPLY ROLL TO INSTRUMENT BASIS ---
    phi = np.radians(roll_deg)
    
    j_inst = (np.cos(phi) * j_cel) + (np.sin(phi) * i_cel)
    i_inst = (np.cos(phi) * i_cel) - (np.sin(phi) * j_cel)
    k_inst = k_cel

    # Reshape for matrix multiplication
    i_inst = i_inst.reshape(3, 1)
    j_inst = j_inst.reshape(3, 1)
    k_inst = k_inst.reshape(3, 1)

    # --- 4. CALCULATE VECTORS ---
    # We use the explicitly broadcasted arrays here so the shapes align perfectly
    theta = np.radians(el_broadcast).reshape(1, -1)
    A     = np.radians(az_broadcast).reshape(1, -1)

    vec_z = (np.sin(theta) * k_inst) + \
            (np.cos(theta) * np.cos(A) * j_inst) + \
            (np.cos(theta) * np.sin(A) * i_inst)

    vec_x = (np.cos(A) * i_inst) - (np.sin(A) * j_inst)

    # --- 5. CONVERT TO RA/DEC ---
    z_ra, z_dec = cartesian_to_spherical(vec_z)
    x_ra, x_dec = cartesian_to_spherical(vec_x)

    # Return dictionary with everything
    return {
        'CAM_AZ': az_broadcast,   # <--- Added
        'CAM_EL': el_broadcast,   # <--- Added
        'CAMZRA': z_ra, 
        'CAMZDEC': z_dec,
        'CAMXRA': x_ra, 
        'CAMXDEC': x_dec
    }

'''
# Inputting distinct lists for both Azimuth and Elevation
# Camera 1: Az=0, El=10
# Camera 2: Az=20, El=30
# Camera 3: Az=40, El=50
az_list = [0, 20, 40]
el_list = [10, 30, 50]

results_conf = get_camera_pointing(
    ra_inst=180.0, 
    dec_inst=10.0, 
    roll_deg=0.0, 
    elevation_deg=el_list, 
    azimuths_deg=az_list
)

for i in range(len(results_conf['CAM_AZ'])):
    az = results_ring['CAM_AZ'][i]
    el = results_ring['CAM_EL'][i]
    ra = results_ring['CAMZRA'][i]
    dec = results_ring['CAMZDEC'][i]
    print(f"{az:<10.1f} | {el:<10.1f} | {ra:<15.4f} | {dec:<15.4f}")
'''