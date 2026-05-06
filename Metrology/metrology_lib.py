import numpy as np
import matplotlib.pyplot as plt
import ezdxf
import cv2
import math
import os

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
#        line_left = cv2.fitLine(left_points, cv2.DIST_L2, 0, 0.001, 0.001)
#        line_right = cv2.fitLine(right_points, cv2.DIST_L2, 0, 0.001, 0.001)
        
        # cv2.DIST_HUBER ignores points that deviate too far from the consensus
        line_left = cv2.fitLine(left_points, cv2.DIST_HUBER, 0, 0.001, 0.001)
        line_right = cv2.fitLine(right_points, cv2.DIST_HUBER, 0, 0.001, 0.001)
    
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
    
        # 1. Get the bounding box to find the approximate top and bottom locations
        rect = cv2.minAreaRect(cnt)
        box = np.int32(cv2.boxPoints(rect))
    
        # 2. Identify the short sides (caps) vs long sides (walls)
        dist_0_1 = math.hypot(box[0][0] - box[1][0], box[0][1] - box[1][1])
        dist_1_2 = math.hypot(box[1][0] - box[2][0], box[1][1] - box[2][1])
    
        if dist_0_1 < dist_1_2:
            cap_mid_1 = np.mean([box[0], box[1]], axis=0)
            cap_mid_2 = np.mean([box[2], box[3]], axis=0)
        else:
            cap_mid_1 = np.mean([box[1], box[2]], axis=0)
            cap_mid_2 = np.mean([box[3], box[0]], axis=0)
        
        # 3. Use your highly accurate width vector as the projection axis
        height_axis = np.array([avg_vx_w, avg_vy_w])
        if np.dot(height_axis, cap_mid_2 - cap_mid_1) < 0:
            height_axis = -height_axis # Ensure it points from cap 1 to cap 2
        
        # 4. Filter contour points falling within a narrow band (+/- 15 pixels)
        tolerance = 20 #15 
        top_points_filtered, bot_points_filtered = [], []
    
        for pt in cnt:
            p = pt[0]
            # Check distance from cap 1
            if abs(np.dot(p - cap_mid_1, height_axis)) <= tolerance:
                top_points_filtered.append(p)
            
            # Check distance from cap 2
            if abs(np.dot(p - cap_mid_2, height_axis)) <= tolerance:
                bot_points_filtered.append(p)
            
        # 5. Fit (average) the filtered points
        if len(top_points_filtered) > 0 and len(bot_points_filtered) > 0:
            # Project all filtered points onto the shared axis to find their true distances
            avg_proj_1 = np.mean([np.dot(p - cap_mid_1, height_axis) for p in top_points_filtered])
            avg_proj_2 = np.mean([np.dot(p - cap_mid_1, height_axis) for p in bot_points_filtered])
        
            ls_height_px = abs(avg_proj_2 - avg_proj_1)
            ls_height = ls_height_px * resolution_y
        
            final_top_center = cap_mid_1 + avg_proj_1 * height_axis
            final_bot_center = cap_mid_1 + avg_proj_2 * height_axis
        else:
            # Safety fallback if the slit ends were completely missing/outside tolerance
            ls_height_px = max(dist_0_1, dist_1_2)
            ls_height = ls_height_px * resolution_y
            final_top_center, final_bot_center = cap_mid_1, cap_mid_2


        # --- STORE IN DICTIONARY ---
        if 0.200 < ls_width < 4.0 and 11.0 < ls_height < 15:
            measured_slits.append({
                "center": (cx * resolution_x, cy * resolution_y),
                "width": ls_width,
                "height": ls_height
            })

        # --- VISUALIZATION ---
        mult = max(h, w) / 1.9
    
        # Draw Width lines (Green: 0, 255, 0)
        cv2.line(output_image, (int(x1 - mult * avg_vx_w), int(y1 - mult * avg_vy_w)), 
                 (int(x1 + mult * avg_vx_w), int(y1 + mult * avg_vy_w)), (0, 255, 0), 1)
        cv2.line(output_image, (int(x2 - mult * avg_vx_w), int(y2 - mult * avg_vy_w)), 
                 (int(x2 + mult * avg_vx_w), int(y2 + mult * avg_vy_w)), (0, 255, 0), 1)

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

def save_top_width_error_subimages(image, measured_slits, dxf_slits, res_x, res_y, output_folder, top_n=5, margin=50):
    """
    Identifies slits with the largest width error, crops them, overlays data, 
    and saves them side-by-side with a highlighted overview of the full image.
    
    Parameters:
    - image: The original microscope image (NumPy array).
    - measured_slits: List of measured slit dictionaries (sorted).
    - dxf_slits: List of nominal DXF slit dictionaries (sorted).
    - res_x, res_y: Resolution conversion factors (mm per pixel).
    - output_folder: Destination folder for the cropped images.
    - top_n: Number of top error slits to process.
    - margin: Pixel padding around the cropped slit.
    """
    os.makedirs(output_folder, exist_ok=True)
    
    # --- NEW: Prepare the overview image ---
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
    
    # 1. Calculate width error for all pairs
    for idx, (m_slit, d_slit) in enumerate(zip(measured_slits, dxf_slits)):
        width_error = abs(m_slit["width"] - d_slit["width"])
        error_data.append({
            "index": idx,
            "error": width_error,
            "measured": m_slit
        })
        
    # 2. Sort the slits by width error in descending order
    error_data.sort(key=lambda x: x["error"], reverse=True)
    
    # 3. Extract and save the top N subimages
    for rank in range(min(top_n, len(error_data))):
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
        text_lines = [
            f"Rank: {rank + 1} (Orig Idx: {data['index']})",
            f"Width Err: {data['error']:.4f} mm",
            f"Fitted W: {m_slit['width']:.4f} mm",
            f"Fitted H: {m_slit['height']:.4f} mm"
        ]
        
        for j, text in enumerate(text_lines):
            y_pos = 30 + (j * 30) 
            cv2.putText(sub_img_color, text, (15, y_pos), 
                        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 255), 2)
            
        # --- NEW: Create the highlighted overview ---
        overview_img = base_overview_color.copy()
        
        # Calculate coordinates for the highlight marker on the scaled overview
        cx_overview = int(cx_px * scale_factor)
        cy_overview = int(cy_px * scale_factor)
        
        # Draw a red circle and crosshair on the overview to highlight the slit
        marker_radius = max(15, int(max(overview_w, overview_h) * 0.02))
        cv2.circle(overview_img, (cx_overview, cy_overview), radius=marker_radius, color=(0, 0, 255), thickness=2)
        cv2.drawMarker(overview_img, (cx_overview, cy_overview), color=(0, 0, 255), 
                       markerType=cv2.MARKER_CROSS, markerSize=marker_radius, thickness=2)
        
        # --- NEW: Combine the Crop and Overview ---
        h_sub, w_sub = sub_img_color.shape[:2]
        h_over, w_over = overview_img.shape[:2]
        
        # Match heights by padding the bottom of the shorter image with black pixels
        max_h = max(h_sub, h_over)
        padded_sub = cv2.copyMakeBorder(sub_img_color, 0, max_h - h_sub, 0, 0, cv2.BORDER_CONSTANT, value=[0, 0, 0])
        padded_over = cv2.copyMakeBorder(overview_img, 0, max_h - h_over, 0, 0, cv2.BORDER_CONSTANT, value=[0, 0, 0])
        
        # Concatenate horizontally: [ Zoomed-in Crop | Full Highlighted Overview ]
        combined_img = np.hstack((padded_sub, padded_over))
        
        # 5. Save the combined layout
        filename = f"error_rank_{rank + 1}_idx_{data['index']}.png"
        save_path = os.path.join(output_folder, filename)
        cv2.imwrite(save_path, combined_img)
        
    print(f"Successfully generated and saved top {top_n} combined error subimages to '{output_folder}'.")