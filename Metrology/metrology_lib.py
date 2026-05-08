import numpy as np
import matplotlib.pyplot as plt
import ezdxf
import cv2
import math
import os
from skimage.transform import EuclideanTransform
from skimage.measure import ransac


def show_image(img, title, scale):
    
    # Calculate the new dimensions based on the original shape
    # img.shape[1] is width, img.shape[0] is height
    new_width = int(img.shape[1] * scale)
    new_height = int(img.shape[0] * scale)
    new_dimensions = (new_width, new_height)

    # Perform the resize using INTER_AREA
    downscaled = cv2.resize(img, new_dimensions, interpolation=cv2.INTER_AREA)
    plt.figure(figsize=(10, 8))
    plt.imshow(downscaled, cmap='gray')
    plt.title(title)
    plt.axis('off')
    plt.show()

def create_clean_binary_mask(img_gray, morph=False, otsu_fac=2):
    """
    Processes a microscope image to cleanly separate dark slits from bright background.
    Returns the original grayscale image and the cleaned binary image.
    """
    # 1. Load the image in grayscale
    #img_gray = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if img_gray is None:
        raise ValueError("Image not found or unable to load.")

    # 2. Apply CLAHE to equalize uneven microscope illumination
    # tileGridSize divides the image into 3x3 blocks for local equalization
    clahe = cv2.createCLAHE(clipLimit=1.0, tileGridSize=(3, 3))
    equalized = clahe.apply(img_gray)

    # 3. Denoise while preserving sharp geometric edges
    # The kernel size (5) must be an odd number. Increase to 7 if the image is very noisy.
    blurred = cv2.medianBlur(equalized, 7)

    # 4. Automatic Thresholding using Otsu's method
    # THRESH_BINARY_INV makes the dark slits white (255) and the bright background black (0)
    _, binary = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + otsu_fac * cv2.THRESH_OTSU)

    if morph:

        # 5. Morphological Cleanup
        kernel = np.ones((3, 3), np.uint8)
    
        # Opening: removes small isolated noise in the background
        clean_binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, kernel, iterations=1)
    
        # Closing: fills in small pinholes or imperfections inside the slits
        clean_binary = cv2.morphologyEx(clean_binary, cv2.MORPH_CLOSE, kernel, iterations=1)
    
    else:
        clean_binary = binary

    return clean_binary


'''
CONTOUR ANALYSIS
'''

#############################################################################################################################
# WIDTH
#############################################################################################################################
def robust_fit_line(contour_pts, side, distance_threshold=2.0, max_iters=10):
    """
    Iteratively fits a line using an asymmetric filter. 
    If side is 'left', it aggressively rejects points far to the left (stains).
    If side is 'right', it aggressively rejects points far to the right (stains).
    Returns the line, the clean inliers, and the rejected outliers.
    """
    if len(contour_pts) < 5:
        line = cv2.fitLine(np.array(contour_pts, dtype=np.float32), cv2.DIST_L2, 0, 0.01, 0.01)
        return line, contour_pts, np.array([])
        
    pts = np.array(contour_pts, dtype=np.float32).reshape(-1, 2)
    original_pts = pts.copy()
    
    for _ in range(max_iters):
        # Fit current line
        line = cv2.fitLine(pts, cv2.DIST_L2, 0, 0.01, 0.01)
        vx, vy, x0, y0 = line[0][0], line[1][0], line[2][0], line[3][0]
        
        # Prevent division by zero for perfectly horizontal lines (safety check)
        if abs(vy) < 1e-5:
            vy = 1e-5 if vy >= 0 else -1e-5
            
        # Calculate the expected X coordinate on the line for every point's Y coordinate
        expected_x = x0 + (vx / vy) * (pts[:, 1] - y0)
        
        # Calculate signed horizontal distance (dx). 
        # Negative dx = point is to the LEFT of the line. Positive dx = point is to the RIGHT.
        dx = pts[:, 0] - expected_x
        
        if side == 'left':
            # True edge is on the right, stains bulge to the left (negative dx).
            # Keep all points to the right, and allow only a small tolerance to the left.
            inlier_mask = dx >= -distance_threshold
        elif side == 'right':
            # True edge is on the left, stains bulge to the right (positive dx).
            # Keep all points to the left, and allow only a small tolerance to the right.
            inlier_mask = dx <= distance_threshold
        else:
            inlier_mask = np.abs(dx) <= distance_threshold
            
        if np.all(inlier_mask) or np.sum(inlier_mask) < 5:
            break
            
        pts = pts[inlier_mask]
        
    # Final fit on the clean inliers
    final_line = cv2.fitLine(pts, cv2.DIST_L2, 0, 0.01, 0.01)
    
    # Calculate final outliers for visualization
    vx, vy, x0, y0 = final_line[0][0], final_line[1][0], final_line[2][0], final_line[3][0]
    if abs(vy) < 1e-5: vy = 1e-5 if vy >= 0 else -1e-5
        
    expected_x_orig = x0 + (vx / vy) * (original_pts[:, 1] - y0)
    dx_orig = original_pts[:, 0] - expected_x_orig
    
    if side == 'left':
        outliers = original_pts[dx_orig < -distance_threshold]
    elif side == 'right':
        outliers = original_pts[dx_orig > distance_threshold]
    else:
        outliers = original_pts[np.abs(dx_orig) > distance_threshold]
    
    return final_line, pts, outliers


#############################################################################################################################
# HEIGHT
#############################################################################################################################


def robust_physical_inside_out(projections, side, expected_edge_len_px, distance_threshold=3.0, min_coverage=0.40):
    """
    Sweeps from the INSIDE out.
    Stops at the first cluster that contains an absolute number of points based on 
    the physical expected length of the edge. Completely ignores massive stain counts.
    """
    projs = np.array(projections, dtype=np.float32)
    if len(projs) == 0:
        return 0.0, []

    # Calculate absolute number of pixels required to be considered a "true solid edge"
    required_pts = max(5, int(expected_edge_len_px * min_coverage))

    # Sort from INSIDE to OUTSIDE
    if side == 'origin':
        sorted_projs = np.sort(projs)[::-1] 
    else:
        sorted_projs = np.sort(projs)       

    edge_val = None

    for y in sorted_projs:
        # Count absolute number of points within the threshold
        inliers_mask = np.abs(projs - y) <= distance_threshold
        if np.sum(inliers_mask) >= required_pts:
            edge_val = np.mean(projs[inliers_mask])
            break

    if edge_val is None:
        # Fallback to the densest peak if something goes horribly wrong
        best_y = np.median(projs)
        max_inliers = -1
        for y in projs:
            count = np.sum(np.abs(projs - y) <= distance_threshold)
            if count > max_inliers:
                max_inliers = count
                best_y = y
        edge_val = best_y

    # Refine for sub-pixel accuracy
    diffs = np.abs(projs - edge_val)
    inliers = projs[diffs <= distance_threshold]
    final_val = np.mean(inliers) if len(inliers) > 0 else edge_val

    # Everything further outward is a stain
    if side == 'origin':
        outliers = projs < (final_val - distance_threshold)
    else:
        outliers = projs > (final_val + distance_threshold)

    return final_val, outliers

#############################################################################################################################

def process_contours(image, contours, resolution_x, resolution_y):
    # Create a copy to draw the fitted lines on (Green for width, Red for height)
    output_image = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)

    print("Image dims", output_image.shape )
    print("Contours dim", len(contours))

    # This will store our final dictionaries
    measured_slits = []

    for i, cnt in enumerate(contours):
        # 1. Get bounding box to find the center and define the 5%/95% trim zones
        x, y, w, h = cv2.boundingRect(cnt)
        cx = x + w / 2.0
        cy = y + h / 2.0
    
        # Trim zones for Vertical Lines (Width)
        y_min = y + (0.05 * h)
        y_max = y + (0.95 * h)
    
        left_points, right_points = [], []
    
        # 2. Sort points into left/right edges (We no longer need to sort top/bottom here)
        for pt in cnt:
            px, py = pt[0]
        
            # Sort for Width (Left vs Right)
            if y_min < py < y_max: 
                if px < cx:
                    left_points.append([px, py])
                else:
                    right_points.append([px, py])
                
        # Convert to numpy arrays for OpenCV
        left_points = np.array(left_points, dtype=np.int32)
        right_points = np.array(right_points, dtype=np.int32)
    
        # --- FIT WIDTH (Vertical Edges) ---
        #line_left = cv2.fitLine(left_points, cv2.DIST_L2, 0, 0.001, 0.001)
        #line_right = cv2.fitLine(right_points, cv2.DIST_L2, 0, 0.001, 0.001)
        
        # cv2.DIST_HUBER ignores points that deviate too far from the consensus
        #line_left = cv2.fitLine(left_points, cv2.DIST_HUBER, 0, 0.001, 0.001)
        #line_right = cv2.fitLine(right_points, cv2.DIST_HUBER, 0, 0.001, 0.001)

        # Iteratively fit and discard outlier points (asymmetrically targeting stains)
        line_left, clean_left, outliers_left = robust_fit_line(left_points, side='left', distance_threshold=3.0)
        line_right, clean_right, outliers_right = robust_fit_line(right_points, side='right', distance_threshold=3.0)
    
        vx1, vy1, x1, y1 = line_left[0][0], line_left[1][0], line_left[2][0], line_left[3][0]
        vx2, vy2, x2, y2 = line_right[0][0], line_right[1][0], line_right[2][0], line_right[3][0]

        # Force them to be perfectly parallel
        if (vx1 * vx2 + vy1 * vy2) < 0: 
            vx2, vy2 = -vx2, -vy2 # Align vectors
        
        avg_vx_w = (vx1 + vx2) / 2.0
        avg_vy_w = (vy1 + vy2) / 2.0
        norm_w = math.sqrt(avg_vx_w**2 + avg_vy_w**2)
        avg_vx_w, avg_vy_w = avg_vx_w / norm_w, avg_vy_w / norm_w
    
        ls_width_px = abs((x2 - x1) * (-avg_vy_w) + (y2 - y1) * avg_vx_w)
        ls_width = ls_width_px * resolution_x

# --- ROI BAND FILTERING HEIGHT FIT ---
    
        # 1. Use the raw contour directly (CHAIN_APPROX_NONE is naturally dense)
        pts = cnt.reshape(-1, 2).astype(np.float32)
        
        # 2. Setup projection axes relative to the center of the slit
        height_axis = np.array([avg_vx_w, avg_vy_w])
        ortho_axis = np.array([-height_axis[1], height_axis[0]])
        center_pt = np.array([cx, cy])
        
        # 3. Project all points onto the axes
        vecs = pts - center_pt
        h_projs = np.dot(vecs, height_axis)
        w_projs = np.dot(vecs, ortho_axis)
        
# 4. Use a narrow core mask (25% of half-width)
        core_w_limit = (ls_width_px / 2.0) * 0.25
        core_mask = np.abs(w_projs) <= core_w_limit
        
        # Calculate the absolute pixel width of this tunnel
        tunnel_width_px = core_w_limit * 2.0
        
        core_pts = pts[core_mask]
        core_h_projs = h_projs[core_mask]
        
        if len(core_h_projs) > 5:
            # 5. Split into Top and Bottom
            top_mask = core_h_projs < 0
            bot_mask = core_h_projs > 0
            
            top_projs = core_h_projs[top_mask]
            bot_projs = core_h_projs[bot_mask]
            
            top_pts_filtered = core_pts[top_mask]
            bot_pts_filtered = core_pts[bot_mask]
            
            # 6. Apply the Physical Inside-Out Caliper
            # Require the edge to physically span at least 40% of the tunnel's width
            avg_proj_1, top_outliers_mask = robust_physical_inside_out(top_projs, side='origin', expected_edge_len_px=tunnel_width_px, distance_threshold=3.0, min_coverage=0.40)
            avg_proj_2, bot_outliers_mask = robust_physical_inside_out(bot_projs, side='far_end', expected_edge_len_px=tunnel_width_px, distance_threshold=3.0, min_coverage=0.40)
            
            # 7. Final calculations
            ls_height_px = abs(avg_proj_2 - avg_proj_1)
            ls_height = ls_height_px * resolution_y
            
            final_top_center = center_pt + avg_proj_1 * height_axis
            final_bot_center = center_pt + avg_proj_2 * height_axis
            
            # Extract outliers for drawing
            outliers_top = [top_pts_filtered[k] for k, is_out in enumerate(top_outliers_mask) if is_out]
            outliers_bot = [bot_pts_filtered[k] for k, is_out in enumerate(bot_outliers_mask) if is_out]
        else:
            # Absolute fallback if shape is entirely destroyed
            ls_height_px = h
            ls_height = ls_height_px * resolution_y
            final_top_center = center_pt - (h/2) * height_axis
            final_bot_center = center_pt + (h/2) * height_axis
            outliers_top, outliers_bot = [], []


        # --- STORE IN DICTIONARY ---
        if 0.200 < ls_width < 4.0 and 11.0 < ls_height < 15:
            measured_slits.append({
                "center": (cx * resolution_x, cy * resolution_y),
                "width": ls_width,
                "height": ls_height
            })
        else:
            # Print the rejected values so we know exactly why they were dropped
            print(f"Dropped slit near pixel ({cx:.0f}, {cy:.0f}) | Width: {ls_width:.3f} mm | Height: {ls_height:.3f} mm")

        # --- VISUALIZATION ---
        mult = max(h, w) / 1.9
    
        # Draw Width lines (Green: 0, 255, 0)
        cv2.line(output_image, (int(x1 - mult * avg_vx_w), int(y1 - mult * avg_vy_w)), 
                 (int(x1 + mult * avg_vx_w), int(y1 + mult * avg_vy_w)), (0, 255, 0), 1)
        cv2.line(output_image, (int(x2 - mult * avg_vx_w), int(y2 - mult * avg_vy_w)), 
                 (int(x2 + mult * avg_vx_w), int(y2 + mult * avg_vy_w)), (0, 255, 0), 1)

        # Draw rejected "stain" points in bright yellow (BGR: 0, 255, 255)
        for pt in outliers_left:
            cv2.circle(output_image, (int(pt[0]), int(pt[1])), 1, (0, 255, 255), -1)
        for pt in outliers_right:
            cv2.circle(output_image, (int(pt[0]), int(pt[1])), 1, (0, 255, 255), -1)       

        # Draw rejected height "stain" points in orange (BGR: 0, 165, 255)
        for pt in outliers_top:
            cv2.circle(output_image, (int(pt[0]), int(pt[1])), 1, (0, 165, 255), -1)
        for pt in outliers_bot:
            cv2.circle(output_image, (int(pt[0]), int(pt[1])), 1, (0, 165, 255), -1)

        # --- DRAW TIGHT HEIGHT LINES (Red) ---
        # We use half of the actual pixel width we calculated to extend the line
        # left and right from the center anchor point.
        half_w = ls_width_px / 2.0 
    
        ortho_vx, ortho_vy = -avg_vy_w, avg_vx_w
    
        # Calculate the exact start and end points for the TOP cap
        top_start = (int(final_top_center[0] - half_w * ortho_vx), 
                     int(final_top_center[1] - half_w * ortho_vy))
        top_end   = (int(final_top_center[0] + half_w * ortho_vx), 
                     int(final_top_center[1] + half_w * ortho_vy))
                 
        # Calculate the exact start and end points for the BOTTOM cap
        bot_start = (int(final_bot_center[0] - half_w * ortho_vx), 
                     int(final_bot_center[1] - half_w * ortho_vy))
        bot_end   = (int(final_bot_center[0] + half_w * ortho_vx), 
                     int(final_bot_center[1] + half_w * ortho_vy))
    
        # Draw them
        cv2.line(output_image, top_start, top_end, (0, 0, 255), 1)
        cv2.line(output_image, bot_start, bot_end, (0, 0, 255), 1)

    return output_image, measured_slits







'''
DXF ANALYSIS
'''
def process_dxf(dxf_path):
    """Parses the DXF file and extracts the nominal dimensions of the slits."""
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()
    
    dxf_slits = []
    
    # Query for all lightweight polylines in the modelspace
    for polyline in msp.query("LWPOLYLINE"):
        # The file contains closed polylines for the slits
        # ezdxf represents the vertices as a list of points
        points = polyline.get_points('xy')
        
        # A rectangular slit should have 4 or 5 points (if it explicitly closes back on the first point)
        if len(points) >= 4:
            x_coords = [p[0] for p in points]
            y_coords = [p[1] for p in points]
            
            # Calculate the bounding box of the polyline
            min_x, max_x = min(x_coords), max(x_coords)
            min_y, max_y = min(y_coords), max(y_coords)
            
            # Calculate nominal dimensions
            width = max_x - min_x
            length = max_y - min_y
            
            # Calculate nominal center points for later alignment
            center_x = (max_x + min_x) / 2.0
            center_y = (max_y + min_y) / 2.0
            
            # Ensure width is the smaller dimension (to match Phase A)
            real_width = min(width, length)
            real_length = max(width, length)
            
            # Filter out random tiny artifacts if any exist
            if real_width > 0.1 and real_length > 0.1:
                dxf_slits.append({
                    "center": (center_x, center_y),
                    "width": real_width,
                    "height": real_length
                })
            
    print(f"Found {len(dxf_slits)} nominal slits in the DXF.")
    return dxf_slits

def sort_grid(slits, y_direction='down', row_tolerance=0.5):
    """
    Sorts a list of slit dictionaries Top-to-Bottom, Left-to-Right.
    
    :param slits: List of dictionaries containing 'center': (x, y)
    :param y_direction: 'down' for Images (Y increases downwards), 'up' for DXF (Y increases upwards)
    :param row_tolerance: The maximum difference in Y to still be considered the same row (in mm)
    """
    if not slits:
        return []

    # 1. Sort by Y-axis to group into rows
    # If y_direction is 'up' (DXF), the highest Y is the top, so we reverse the sort.
    # If y_direction is 'down' (Image), the lowest Y is the top, so normal sort.
    reverse_y = (y_direction == 'up') 
    sorted_by_y = sorted(slits, key=lambda s: s['center'][1], reverse=reverse_y)

    sorted_grid = []
    current_row = [sorted_by_y[0]]
    
    # 2. Iterate through and cluster into rows based on tolerance
    for slit in sorted_by_y[1:]:
        prev_y = current_row[-1]['center'][1]
        curr_y = slit['center'][1]
        
        # If the Y difference is small, it's in the same row
        if abs(curr_y - prev_y) <= row_tolerance:
            current_row.append(slit)
        else:
            # Row is complete. Sort this row by X (Left to Right)
            current_row.sort(key=lambda s: s['center'][0])
            sorted_grid.extend(current_row)
            # Start a new row
            current_row = [slit]
            
    # 3. Sort and append the final row
    if current_row:
        current_row.sort(key=lambda s: s['center'][0])
        sorted_grid.extend(current_row)
        
    return sorted_grid


###########################################################################################
# TRANSFORMATIONS
###########################################################################################

def find_affine_transformation(physical, nominal):
    # 1. Extract the (x, y) tuples into standard N x 2 numpy arrays
    pts_phys = np.array([item["center"] for item in physical], dtype=np.float32)
    pts_nom  = np.array([item["center"] for item in nominal], dtype=np.float32)

    # 2. Find the optimal transformation matrix
    # transform_matrix will map physical points -> nominal points
    transform_matrix, inliers = cv2.estimateAffinePartial2D(pts_phys, pts_nom)

    if transform_matrix is not None:
        print("Transformation Matrix:\n", transform_matrix)
    else:
        print("Failed to find a valid transformation.")

    # 2.Transform sorted_physical
    # 1. Reshape your (N, 2) array to (N, 1, 2) for OpenCV
    pts_phys_reshaped = pts_phys.reshape(-1, 1, 2)

    # 2. Apply the 2x3 transformation matrix
    aligned_pts_reshaped = cv2.transform(pts_phys_reshaped, transform_matrix)

    # 3. Reshape back to a standard (N, 2) array
    aligned_pts = aligned_pts_reshaped.reshape(-1, 2)


    aligned_physical = []
    residual_errors = []

    for i in range(len(physical)):
        # 1. Extract the new aligned coordinates and the nominal target
        new_x, new_y = aligned_pts[i]
        nom_x, nom_y = nominal[i]["center"]
        
        # 2. Calculate the exact difference (Error)
        dx = new_x - nom_x
        dy = new_y - nom_y
        euclidean_dist = np.hypot(dx, dy)
        
        # 3. Store the new aligned data
        aligned_physical.append({
            "center": (new_x, new_y),
            "width": physical[i]["width"],
            "height": physical[i]["height"]
        })
        
        # 4. Store the errors for statistical analysis later
        residual_errors.append({
            "x_error": dx,
            "y_error": dy,
            "total_offset": euclidean_dist
        })

    # Quick readout of the worst-case scenario
    max_error = max(item["total_offset"] for item in residual_errors)
    print(f"Maximum remaining positional error: {max_error:.4f} units")

    return aligned_physical, residual_errors


def find_rigid_transformation(physical, nominal):
    #this is a rigid body fit (no scaling)

    # 1. Format your points as standard (N, 2) numpy arrays
    pts_phys = np.array([item["center"] for item in physical], dtype=np.float32)
    pts_nom  = np.array([item["center"] for item in nominal], dtype=np.float32)

    # 2. Estimate the pure Rigid Transform (Scale is strictly 1.0)
    # residual_threshold is the max pixel/unit error to be considered a "good" match
    model, inliers = ransac((pts_phys, pts_nom), EuclideanTransform, 
                            min_samples=2, residual_threshold=2.0, max_trials=100)

    if model is not None:
        # This returns a 3x3 homogeneous transformation matrix
        matrix_3x3 = model.params
        
        # You can easily extract the exact rotation angle and translation
        rotation_rads = model.rotation
        tx, ty = model.translation
        
        print(f"Rotation (degrees): {np.degrees(rotation_rads):.4f}")
        print(f"Translation: X={tx:.4f}, Y={ty:.4f}")
    else:
        print("Failed to find a valid transformation.")

    aligned_pts_rigid = model(pts_phys)

    aligned_physical = []
    residual_errors = []
    for i in range(len(physical)):
        # 1. Extract the new aligned coordinates and the nominal target
        new_x, new_y = aligned_pts_rigid[i]
        nom_x, nom_y = nominal[i]["center"]
        
        # 2. Calculate the exact difference (Error)
        dx = new_x - nom_x
        dy = new_y - nom_y
        euclidean_dist = np.hypot(dx, dy)
        
        # 3. Store the new aligned data
        aligned_physical.append({
            "center": (new_x, new_y),
            "width": physical[i]["width"],
            "height": physical[i]["height"]
        })
        
        # 4. Store the errors for statistical analysis later
        residual_errors.append({
            "x_error": dx,
            "y_error": dy,
            "total_offset": euclidean_dist
        })

    # Quick readout of the worst-case scenario
    max_error = max(item["total_offset"] for item in residual_errors)
    print(f"Maximum remaining positional error: {max_error:.4f} units")

    return aligned_physical, residual_errors



def calibrate_resolutions(measured_slits, dxf_slits, current_res_x, current_res_y):
    """
    Calculates the optimal MM_PER_PIXEL calibration by fitting a linear 
    regression between raw image pixel coordinates and nominal DXF coordinates.
    
    Parameters:
    - measured_slits: Sorted list of measured slit dictionaries.
    - dxf_slits: Sorted list of nominal DXF slit dictionaries.
    - current_res_x, current_res_y: The initial resolutions used to generate measured_slits.
    
    Returns:
    - optimal_res_x, optimal_res_y: The corrected calibration factors.
    """
    if len(measured_slits) != len(dxf_slits):
        print("Warning: Mismatch in number of measured and DXF slits.")
        
    pixel_x = []
    pixel_y = []
    nominal_x = []
    nominal_y = []
    
    for m_slit, d_slit in zip(measured_slits, dxf_slits):
        # 1. Back-calculate the raw pixel coordinates from the existing measurements
        px = m_slit["center"][0] / current_res_x
        py = m_slit["center"][1] / current_res_y
        
        # 2. Grab the nominal mm coordinates
        nx = d_slit["center"][0]
        ny = d_slit["center"][1]
        
        pixel_x.append(px)
        pixel_y.append(py)
        nominal_x.append(nx)
        nominal_y.append(ny)
        
    # Convert to numpy arrays
    pixel_x = np.array(pixel_x)
    pixel_y = np.array(pixel_y)
    nominal_x = np.array(nominal_x)
    nominal_y = np.array(nominal_y)
    
    # 3. Perform linear regression: nominal_mm = (optimal_res) * pixel + offset
    # np.polyfit(x, y, 1) returns [slope, intercept]
    fit_x = np.polyfit(pixel_x, nominal_x, 1)
    fit_y = np.polyfit(pixel_y, nominal_y, 1)
    
    optimal_res_x = abs(fit_x[0])  # Use abs() in case axes are inverted
    optimal_res_y = abs(fit_y[0])
    
    print("--- Calibration Results ---")
    print(f"Initial Res X: {current_res_x:.9f} mm/px")
    print(f"Optimal Res X: {optimal_res_x:.9f} mm/px")
    print(f"Error X:       {abs(optimal_res_x - current_res_x) / current_res_x * 100:.4f} %")
    print("-" * 25)
    print(f"Initial Res Y: {current_res_y:.9f} mm/px")
    print(f"Optimal Res Y: {optimal_res_y:.9f} mm/px")
    print(f"Error Y:       {abs(optimal_res_y - current_res_y) / current_res_y * 100:.4f} %")
    
    return optimal_res_x, optimal_res_y

###########################################################################################
# SUBIMAGE SAVING
###########################################################################################

def save_top_error_subimages(image, measured_slits, dxf_slits, res_x, res_y, output_folder, error_type='width', top_n=5, min_error_mm=0.010, margin=50):
    """
    Identifies slits with the largest errors, crops them, overlays data, 
    and saves them side-by-side with a highlighted overview of the full image.
    
    Parameters:
    - image: The original microscope image (NumPy array).
    - measured_slits: List of measured slit dictionaries (sorted).
    - dxf_slits: List of nominal DXF slit dictionaries (sorted).
    - res_x, res_y: Resolution conversion factors (mm per pixel).
    - output_folder: Destination folder for the cropped images.
    - error_type: 'width' or 'height' to define which measurement to evaluate.
    - top_n: Maximum number of top error slits to process.
    - min_error_mm: Only save subimages if the error is greater than or equal to this value (in mm).
    - margin: Pixel padding around the cropped slit.
    """
    if error_type not in ['width', 'height']:
        raise ValueError("error_type must be either 'width' or 'height'")
        
    os.makedirs(output_folder, exist_ok=True)
    
    # --- Prepare the overview image ---
    orig_h, orig_w = image.shape[:2]
    # Downscale the overview image to a manageable size (max 1000px dimension)
    scale_factor = 1000.0 / max(orig_h, orig_w)
    overview_w = int(orig_w * scale_factor)
    overview_h = int(orig_h * scale_factor)
    
    # Only resize if the image is larger than 1000x1000
    if scale_factor < 1.0:
        base_overview = cv2.resize(image, (overview_w, overview_h), interpolation=cv2.INTER_AREA)
    else:
        base_overview = image.copy()
        scale_factor = 1.0
        
    # Ensure the overview is in BGR color for red highlights
    if len(base_overview.shape) == 2:
        base_overview_color = cv2.cvtColor(base_overview, cv2.COLOR_GRAY2BGR)
    else:
        base_overview_color = base_overview.copy()
    # ---------------------------------------
    
    error_data = []
    
    # 1. Calculate error for all pairs based on the chosen switch
    for idx, (m_slit, d_slit) in enumerate(zip(measured_slits, dxf_slits)):
        if error_type == 'width':
            error_val = abs(m_slit["width"] - d_slit["width"])
        else:
            error_val = abs(m_slit["height"] - d_slit["height"])
            
        # Only keep the data if the error exceeds our minimum threshold
        if error_val >= min_error_mm:
            error_data.append({
                "index": idx,
                "error": error_val,
                "measured": m_slit
            })
        
    # 2. Sort the valid slits by error in descending order
    error_data.sort(key=lambda x: x["error"], reverse=True)
    
    actual_n = min(top_n, len(error_data))
    
    if actual_n == 0:
        print(f"No {error_type} errors found >= {min_error_mm:.4f} mm. No images saved.")
        return
    
    # 3. Extract and save the top N subimages
    for rank in range(actual_n):
        data = error_data[rank]
        m_slit = data["measured"]
        
        # Revert mm physical coordinates back to pixel coordinates
        cx_px = int(m_slit["center"][0] / res_x)
        cy_px = int(m_slit["center"][1] / res_y)
        w_px = int(m_slit["width"] / res_x)
        h_px = int(m_slit["height"] / res_y)
        
        # Calculate bounding box coordinates with the requested margin
        y1 = max(0, cy_px - (h_px // 2) - margin)
        y2 = min(image.shape[0], cy_px + (h_px // 2) + margin)
        x1 = max(0, cx_px - (w_px // 2) - margin)
        x2 = min(image.shape[1], cx_px + (w_px // 2) + margin)
        
        # Crop the region of interest (ROI)
        sub_img = image[y1:y2, x1:x2]
        
        if len(sub_img.shape) == 2:
            sub_img_color = cv2.cvtColor(sub_img, cv2.COLOR_GRAY2BGR)
        else:
            sub_img_color = sub_img.copy()
            
        # 4. Superimpose the text on the cropped image
        err_label = "Width Err" if error_type == 'width' else "Height Err"
        
        text_lines = [
            f"Rank: {rank + 1} (Orig Idx: {data['index']})",
            f"{err_label}: {data['error']:.4f} mm",
            f"Fitted W: {m_slit['width']:.4f} mm",
            f"Fitted H: {m_slit['height']:.4f} mm"
        ]
        
        for j, text in enumerate(text_lines):
            y_pos = 30 + (j * 30) 
            cv2.putText(sub_img_color, text, (15, y_pos), 
                        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 255), 2)
            
        # --- Create the highlighted overview ---
        overview_img = base_overview_color.copy()
        
        # Calculate coordinates for the highlight marker on the scaled overview
        cx_overview = int(cx_px * scale_factor)
        cy_overview = int(cy_px * scale_factor)
        
        # Draw a red circle and crosshair on the overview to highlight the slit
        marker_radius = max(15, int(max(overview_w, overview_h) * 0.02))
        cv2.circle(overview_img, (cx_overview, cy_overview), radius=marker_radius, color=(0, 0, 255), thickness=2)
        cv2.drawMarker(overview_img, (cx_overview, cy_overview), color=(0, 0, 255), 
                       markerType=cv2.MARKER_CROSS, markerSize=marker_radius, thickness=2)
        
        # --- Combine the Crop and Overview ---
        h_sub, w_sub = sub_img_color.shape[:2]
        h_over, w_over = overview_img.shape[:2]
        
        # Match heights by padding the bottom of the shorter image with black pixels
        max_h = max(h_sub, h_over)
        padded_sub = cv2.copyMakeBorder(sub_img_color, 0, max_h - h_sub, 0, 0, cv2.BORDER_CONSTANT, value=[0, 0, 0])
        padded_over = cv2.copyMakeBorder(overview_img, 0, max_h - h_over, 0, 0, cv2.BORDER_CONSTANT, value=[0, 0, 0])
        
        # Concatenate horizontally: [ Zoomed-in Crop | Full Highlighted Overview ]
        combined_img = np.hstack((padded_sub, padded_over))
        
        # 5. Save the combined layout
        filename = f"{error_type}_error_rank_{rank + 1}_idx_{data['index']}.png"
        save_path = os.path.join(output_folder, filename)
        cv2.imwrite(save_path, combined_img)
        
    print(f"Successfully generated and saved {actual_n} '{error_type}' error subimages (>= {min_error_mm:.4f} mm) to '{output_folder}'.")