# CodedMasks
A collection of (quick&dirty) utilities for coded mask images.
Comes with no warranty (and with a punk fashioned code...)

## Code Utils
### Tools to generate URA, MURA and other codes
[.\mask_utils\code_utils.py](#code_utils.py)
- [next_prime](#next_prime)
- [ura_mura](#ura_mura)
- [bura](#bura)
- [bura33](#bura33)
- [cura](#cura)
- [is_cyclic_difference_set](#is_cyclic_difference_set)

## Fits Utils
### Tools to read and write WFM fits
[.\mask_utils\fits_utils.py](#fits_utils.py)
- [read_fits_events](#read_fits_events)
- [read_mask_bulk](#read_mask_bulk)
- [write_mask_fits](#write_mask_fits)
- [fits_mask_to_dxf](#fits_mask_to_dxf)

## Image Utils
### Tools to manipulate images
[.\mask_utils\image_utils.py](#image_utils.py)
- [shift](#shift)
- [fshift](#fshift)
- [erosion](#erosion)
- [pad_array](#pad_array)
- [float_gcd](#float_gcd)
- [upscale](#upscale)

## Imaging Utils
### Tools to perform decoding operations and calculate system properties (effective area, solid angle, etc.)
[.\mask_utils\imaging_utils.py](#imaging_utils.py)
- [get_openfraction](#get_openfraction)
- [get_angular_res](#get_angular_res)
- [get_coding_power](#get_coding_power)
- [open_fraction_vs_off_axis](#open_fraction_vs_off_axis)
- [eff_area_vs_off_axis](#eff_area_vs_off_axis)
- [snr_vs_off_axis](#snr_vs_off_axis)
- [omega_plate_offaxis](#omega_plate_offaxis)
- [solid_angle](#solid_angle)
- [get_detimage](#get_detimage)
- [decode](#decode)
- [decode_var](#decode_var)
- [get_skysign](#get_skysign)
- [get_detimage_edges](#get_detimage_edges)
- [get_skycoords](#get_skycoords)
- [generate_bulk](#generate_bulk)

## Other Utils
### Miscellaneous tools
[.\mask_utils\other_utils.py](#other_utils.py)
- [filter_source](#filter_source)

<hr style="border:2px solid">

# Mask Utils documentation

<a id='code_utils.py'></a>
## Code Utils

This module provides functions for generating different types of coded mask arrays, specifically Uniformly Redundant Arrays (URA), Modified Uniformly Redundant Arrays (MURA), Biquadratic Uniformly Redundant Arrays (BURA), and a specific BURA variant with ~0.33 open fraction.

### <a id='next_prime'></a> `next_prime(n)`

Finds the next prime number greater than or equal to `n`.

**Parameters:**

* `n` (`int`): The starting number.

**Returns:**

* `int`: The next prime number.

### <a id='ura_mura'></a> `ura_mura(p)`

Generates a Uniformly Redundant Array (URA) or Modified Uniformly Redundant Array (MURA) based on a prime number `p`.

**Parameters:**

* `p` (`int`): A prime number.

**Returns:**

* `numpy.ndarray`: A 1D array representing the URA or MURA.

**Raises:**

* `TypeError`: If `p` is not a prime number.

### <a id='bura'></a> `bura(p, modified=False)`

Generates a Biquadratic Uniformly Redundant Array (BURA).

**Parameters:**

* `p` (`int`): A prime number that fulfills the requirement $p = 4x^2 + 1$ with $x$ being an odd integer.
* `modified` (`bool`, optional): If `True`, generates a Modified BURA (M-BURA) where the first element is set to 0. Defaults to `False`.

**Returns:**

* `numpy.ndarray`: A 1D array representing the BURA or M-BURA.

**Raises:**

* `TypeError`: If `p` is not prime or does not fulfill the condition $p = 4x^2 + 1$ with $x$ odd.

### <a id='bura33'></a> `bura33(p)`

Generates a biquadratic URA with an open fraction of approximately 0.33.

**Note:** The resulting code might not be "perfect".
**Reference:** Based on Baumert L. D. 1971, Lecture Notes in Mathematics No. 182, Cyclic Difference Sets, Theorem 5.18 (iii).
The biquadratic residues of primes $v = 4x^2 + 9$, where $x$ is odd, form a difference set with ~0.33 open fraction.
**Example primes:** 13, 109, 1453, 3373, 3853, 4909, 6733.

**Parameters:**

* `p` (`int`): A prime number that fulfills the requirement $p = 4x^2 + 9$ with $x$ being an odd integer.

**Returns:**

* `numpy.ndarray`: A 1D array representing the biquadratic URA.

**Raises:**

* `TypeError`: If `p` is not prime or does not fulfill the condition $p = 4x^2 + 9$ with $x$ odd.

### <a id='ura_mura'></a> `ura_mura(p)`

Generates a Uniformly Redundant Array (URA) or Modified Uniformly Redundant Array (MURA) based on a prime number `p`.

**Parameters:**

* `p` (`int`): A prime number.

**Returns:**

* `numpy.ndarray`: A 1D array representing the URA or MURA.

### <a id='is_cyclic_difference_set'></a> `is_cyclic_difference_set(s, v, k, lambda_val)`

Tests if a given set is a cyclic (v, k, lambda) difference set.

**Parameters:**

* `s` (`numpy.ndarray`): The set to be tested.
* `v` (`int`): The modulus (size of the cyclic group).
* `k` (`int`): The expected size of the set D..
* `lambda_val` (`int`): The expected number of times each non-zero residue appears as a difference.
  
**Returns:**

* `(bool, str)`: Indicating whether the set is a difference set and a message.
<a id='fits_utils.py'></a>
## FITS Utils

This module provides functions for reading and writing FITS files, particularly for mask data, and converting FITS mask data to DXF format.

### <a id='read_fits_events'></a> `read_fits_events(filein, header0=False, header1=False, verbose=False)`

Reads event data from a FITS file.

**Parameters:**

* `filein` (`str`): The path to the input FITS file.
* `header0` (`bool`, optional): If `True`, returns the primary header. Defaults to `False`.
* `header1` (`bool`, optional): If `True`, returns the header of the first extension (events). Defaults to `False`.
* `verbose` (`bool`, optional): If `True`, prints FITS file information using `hdu_list.info()`. Defaults to `False`.

**Returns:**

* `astropy.table.Table`: The event data.
* `astropy.io.fits.Header`, optional: The primary header if `header0` is `True`.
* `astropy.io.fits.Header`, optional: The header of the first extension if `header1` is `True`.

### <a id='read_mask_bulk'></a> `read_mask_bulk(fitsfile, ext, header_out=False, verbose=False)`

Reads mask data from a WFM mask FITS file (`wfm_mask.fits`).

**Extensions are:**
* `0`: PRIMARY
* `1`: OR_MASK
* `2`: MASK
* `3`: RMATRIX
* `4`: SENS

**Parameters:**

* `fitsfile` (`str`): The path to the WFM mask FITS file.
* `ext` (`int` or `str`): The extension number or name to read.
* `header_out` (`bool`, optional): If `True`, returns the header of the specified extension. Defaults to `False`.
* `verbose` (`bool`, optional): If `True`, prints FITS file information using `hdu_list.info()`. Defaults to `False`.

**Returns:**

* `numpy.ndarray`: A 2D array representing the mask data.
* `astropy.io.fits.Header`, optional: The header of the specified extension if `header_out` is `True`.

### <a id='write_mask_fits'></a> `write_mask_fits(fitsfile, mask, rmatrix, bulk, props)`

Writes a FITS file with mask-related extensions.

**The FITS file will have the following extensions:**
* `1`: OR_MASK
* `2`: MASK
* `3`: RMATRIX
* `4`: SENS

**Parameters:**

* `fitsfile` (`str`): The path where the FITS file will be written.
* `mask` (`numpy.ndarray`): The mask array.
* `rmatrix` (`numpy.ndarray`): The R-matrix array.
* `bulk` (`numpy.ndarray`): The bulk array.
* `props` (`dict`): A dictionary containing properties to be written into the FITS headers, such as 'ELXN', 'ELYDIM', 'ELXDIM', 'MXDIM', 'MYDIM'.

### <a id='fits_mask_to_dxf'></a> `fits_mask_to_dxf(fitsin, dxfout)`

Writes a DXF file of the mask with all open elements as polylines.

**Parameters:**

* `fitsin` (`str`): The path to the input FITS mask file.
* `dxfout` (`str`): The path where the DXF file will be written.

<a id='image_utils.py'></a>
## Image Utils

This module provides various utility functions for image manipulation, including shifting, erosion, padding, and upscaling.

### <a id='shift'></a> `shift(arr, lag)`

Shifts a 2D array horizontally (along `axis=1`) by a given `lag`. Elements shifted out are replaced with zeros.

**Parameters:**

* `arr` (`numpy.ndarray`): The input 2D array.
* `lag` (`int`): The number of positions to shift. Positive values shift to the right, negative values shift to the left.

**Returns:**

* `numpy.ndarray`: The shifted array.

### <a id='fshift'></a> `fshift(arr, lagx, lagy)`

Shifts a 2D array with fractional shifts using SciPy's `ndimage.shift`.

**Parameters:**

* `arr` (`numpy.ndarray`): The input 2D array.
* `lagx` (`float`): The fractional shift along the x-axis.
* `lagy` (`float`): The fractional shift along the y-axis.

**Returns:**

* `numpy.ndarray`: The fractionally shifted array.

### <a id='erosion'></a> `erosion(arr, cut, step)`

Performs 2D matrix erosion on a mask array to simulate finite thickness effects in shadow projections. It "thins" the mask elements across the columns' direction. The erosion is performed only on the correct side of open (1) mask elements: right side if `cut` is negative (thetaX negative), and left side if `cut` is positive (thetaX positive). The function first erodes all integer bins (replacing 1s with 0s). If `cut` is not an integer, a fractional transparency is applied to the last eroded bin.

**Parameters:**

* `arr` (`numpy.ndarray`): The input mask array.
* `cut` (`float`): The amount of erosion to apply.
* `step` (`float`): The size of one bin in the erosion process.

**Returns:**

* `numpy.ndarray`: The eroded mask array.

### <a id='pad_array'></a> `pad_array(arr, npadx, npady)`

Pads a 2D array with zeros.

**Parameters:**

* `arr` (`numpy.ndarray`): The input 2D array.
* `npadx` (`int`): The number of zeros to pad on each side along the x-axis.
* `npady` (`int`): The number of zeros to pad on each side along the y-axis.

**Returns:**

* `numpy.ndarray`: The padded array.

### <a id='float_gcd'></a> `float_gcd(a, b, rtol=1e-05, atol=1e-05)`

Calculates the greatest common divisor (GCD) for floating-point numbers with specified relative and absolute tolerances.

**Parameters:**

* `a` (`float`): The first floating-point number.
* `b` (`float`): The second floating-point number.
* `rtol` (`float`, optional): Relative tolerance. Defaults to `1e-5`.
* `atol` (`float`, optional): Absolute tolerance. Defaults to `1e-5`.

**Returns:**

* `float`: The greatest common divisor.

### <a id='upscale'></a> `upscale(array, fx, fy)`

Upscales a 2D array by repeating its elements.

**Parameters:**

* `array` (`numpy.ndarray`): The input 2D array.
* `fx` (`int`): The upscaling factor for the x-axis.
* `fy` (`int`): The upscaling factor for the y-axis.

**Returns:**

* `numpy.ndarray`: The upscaled array.

<a id='imaging_utils.py'></a>
## Imaging Utils

This module provides functions for calculating imaging properties such as open fraction, angular resolution, coding power, and effective area, as well as functions for decoding images and handling detector geometry.

### <a id='get_openfraction'></a> `get_openfraction(mask)`

Calculates the open fraction of a mask.

**Parameters:**

* `mask` (`numpy.ndarray`): The mask array.

**Returns:**

* `float`: The open fraction (sum of mask elements divided by total size).

### <a id='get_angular_res'></a> `get_angular_res(m_pitch, d_pitch, m_d_distance, degrees=False)`

Calculates the angular resolution of the imaging system.

**Parameters:**

* `m_pitch` (`float`): Mask pitch.
* `d_pitch` (`float`): Detector pitch.
* `m_d_distance` (`float`): Mask-detector distance.
* `degrees` (`bool`, optional): If `True`, returns the resolution in degrees. Defaults to `False`.

**Returns:**

* `float`: The angular resolution.

### <a id='get_coding_power'></a> `get_coding_power(m_pitch, d_pitch, open_fraction)`

Calculates the coding power based on Skinner 2008.

**Parameters:**

* `m_pitch` (`float`): Mask pitch.
* `d_pitch` (`float`): Detector pitch.
* `open_fraction` (`float`): The open fraction of the mask.

**Returns:**

* `float`: The coding power.

### <a id='open_fraction_vs_off_axis'></a> `open_fraction_vs_off_axis(mask, mask_thickness, mask_x_pitch, mask_y_pitch, thetaX, thetaY, degrees=True)`

Calculates the open fraction of the mask as a function of off-axis angles, considering vignetting.

**Parameters:**

* `mask` (`numpy.ndarray`): The mask array.
* `mask_thickness` (`float`): The thickness of the mask.
* `mask_x_pitch` (`float`): The mask pitch along the x-axis.
* `mask_y_pitch` (`float`): The mask pitch along the y-axis.
* `thetaX` (`float`): The off-axis angle along the x-axis.
* `thetaY` (`float`): The off-axis angle along the y-axis.
* `degrees` (`bool`, optional): If `True`, `thetaX` and `thetaY` are in degrees. Defaults to `True`.

**Returns:**

* `float`: The vignetted open fraction.

### <a id='eff_area_vs_off_axis'></a> `eff_area_vs_off_axis(mask, det, x_pitch_ups, y_pitch_ups, focal, mask_thickness, thetaX, thetaY, degrees=True)`

Calculates the effective area of the system as a function of off-axis angles, accounting for vignetting.

**Parameters:**

* `mask` (`numpy.ndarray`): The mask array.
* `det` (`numpy.ndarray`): The detector array (a binary mask indicating sensitive regions).
* `x_pitch_ups` (`float`): Upscaled x-pitch.
* `y_pitch_ups` (`float`): Upscaled y-pitch.
* `focal` (`float`): Focal length.
* `mask_thickness` (`float`): The thickness of the mask.
* `thetaX` (`float`): The off-axis angle along the x-axis.
* `thetaY` (`float`): The off-axis angle along the y-axis.
* `degrees` (`bool`, optional): If `True`, `thetaX` and `thetaY` are in degrees. Defaults to `True`.

**Returns:**

* `float`: The effective area in the same units as the input pitches and focal length (e.g., mm$^2$ or cm$^2$).

### <a id='snr_vs_off_axis'></a> `snr_vs_off_axis(s_counts, b_counts, mask, bulk, mask_x_pitch, mask_y_pitch, ELXDIM, ELYDIM, det_x_pitch, focal, mask_thickness, thetaX, thetaY, degrees=True, verbose=False)`

Estimate the system sensitivity as a function of the off-axis angle for a given camera geometry and mask code. Makes use of and from G. K. Skinner - "The sensitivity of coded mask telescopes" (2008)

**Parameters:**
  * `s_counts` (`float` or `int`):
    The expected number of source counts (photons from the source).
  * `b_counts` (`float` or `int`):
    The expected number of background counts (noise photons).
  * `mask` (`numpy.ndarray`):
    A 2D NumPy array representing the coded aperture mask. Typically, this is a binary mask (e.g., 0s and 1s).
  * `bulk` (`numpy.ndarray`):
    A 2D NumPy array representing the detector bulk. Its specific structure depends on how the `eff_area_vs_off_axis` function is implemented, but it generally relates to the detector's properties.
  * `mask_x_pitch` (`float`):
    The physical size of a single element (pixel) of the mask along the x-axis, in appropriate units (e.g., mm).
  * `mask_y_pitch` (`float`):
    The physical size of a single element (pixel) of the mask along the y-axis, in appropriate units (e.g., mm).
  * `ELXDIM` (`int`):
    The number of elements (pixels) in the detector along the x-dimension.
  * `ELYDIM` (`int`):
    The number of elements (pixels) in the detector along the y-dimension.
  * `det_x_pitch` (`float`):
    The physical size of a single element (pixel) of the detector along the x-axis, in appropriate units (e.g., mm).
  * `focal` (`float`):
    The focal length of the system (distance between the mask and the detector), in appropriate units (e.g., mm).
  * `mask_thickness` (`float`):
    The physical thickness of the mask material, in appropriate units (e.g., mm). This is important for understanding off-axis attenuation effects.
  * `thetaX` (`float`):
    The off-axis angle in the x-direction. By default, it's assumed to be in degrees.
  * `thetaY` (`float`):
    The off-axis angle in the y-direction. By default, it's assumed to be in degrees.
  * `degrees` (`bool`, optional):
    If `True` (default), `thetaX` and `thetaY` are interpreted as degrees and converted to radians internally. If `False`, they are assumed to be in radians.
  * `verbose` (`bool`, optional):
    If `True`, the function will print intermediate calculation values (on-axis open fraction, coding power, off-axis open fraction, and off-axis area) to the console. Defaults to `False`.

**Returns:**
  * **`float`**:
    The calculated approximate **Signal-to-Noise Ratio (SNR)** for the given off-axis angles and system parameters.


### <a id='omega_plate_offaxis'></a> `omega_plate_offaxis(a, b, d, A, B)`

Calculates the solid angle subtended by a rectangular plate when the line-of-sight hits it off-center. This function is based on equation 34 from [https://vixra.org/pdf/2001.0603v2.pdf](https://vixra.org/pdf/2001.0603v2.pdf).

* `d`: is the distance between the plate and the "observer"
* `a`, `b`: are the dimensions of the rectangular plate
* `A`, `B`: are the corner distances from the observer (where $0 \le A \le a/2$ and $0 \le B \le b/2$)

The calculation considers four sub-rectangles:
* `r1 = A * (b-B)`
* `r2 = (a-A) * (b-B)`
* `r3 = A * B`
* `r4 = (a-A) * B`

**Parameters:**

* `a` (`float`): Dimension of the plate along the x-axis.
* `b` (`float`): Dimension of the plate along the y-axis.
* `d` (`float`): Distance from the plate to the observer.
* `A` (`float`): Corner distance from the observer along the x-axis.
* `B` (`float`): Corner distance from the observer along the y-axis.

**Returns:**

* `float`: The solid angle.

### <a id='solid_angle'></a> `solid_angle(bulk, xstep, ystep, m_d_distance, nobulk=False)`

Calculates the solid angle for each pixel of a bulk array.

**Parameters:**

* `bulk` (`numpy.ndarray`): The bulk array (detector sensitive region).
* `xstep` (`float`): The size of each pixel along the x-axis.
* `ystep` (`float`): The size of each pixel along the y-axis.
* `m_d_distance` (`float`): Mask-detector distance.
* `nobulk` (`bool`, optional): If `True`, calculates solid angle for all pixels regardless of bulk value. Defaults to `False`.

**Returns:**

* `numpy.ndarray`: An array of solid angle values for each pixel.

### <a id='get_detimage'></a> `get_detimage(data, xedges, yedges)`

Generates a 2D histogram (detector image) from event data.

**Parameters:**

* `data` (`numpy.ndarray` or `astropy.table.Table`): The input event data, expected to have 'X' and 'Y' columns.
* `xedges` (`numpy.ndarray`): The bin edges for the x-axis.
* `yedges` (`numpy.ndarray`): The bin edges for the y-axis.

**Returns:**

* `numpy.ndarray`: The 2D detector image.

### <a id='decode'></a> `decode(detimage, rmatrix, bulk)`

Decodes a detector image using the R-matrix and bulk array.

**Parameters:**

* `detimage` (`numpy.ndarray`): The detector image.
* `rmatrix` (`numpy.ndarray`): The R-matrix.
* `bulk` (`numpy.ndarray`): The bulk array.

**Returns:**

* `numpy.ndarray`: The decoded sky image.

### <a id='decode_var'></a> `decode_var(detimage, rmatrix, bulk, m_d_distance, elxdim, elydim)`

Calculates the variance image for decoding, based on the statistical approach outlined in [this document](https://github.com/yuri-evangelista/CodedMasks/blob/main/Coded_aperture_variance_202507.pdf).

**Parameters:**

* `detimage` (`numpy.ndarray`): The detector image.
* `rmatrix` (`numpy.ndarray`): The R-matrix.
* `bulk` (`numpy.ndarray`): The bulk array.
* `m_d_distance` (`float`): Mask-detector distance.
* `elxdim` (`float`): Elemental dimension along X.
* `elydim` (`float`): Elemental dimension along Y.

**Returns:**

* `numpy.ndarray`: The variance image.

### <a id='get_skysign'></a> `get_skysign(skyimage, varimage)`

Calculates the sky significance image from a sky image and its corresponding variance image.

**Parameters:**

* `skyimage` (`numpy.ndarray`): The decoded sky image.
* `varimage` (`numpy.ndarray`): The variance image.

**Returns:**

* `numpy.ndarray`: The sky significance image.

### <a id='get_detimage_edges'></a> `get_detimage_edges(xstep, ystep, nx, ny)`

Generates the bin edges for a detector image.

**Parameters:**

* `xstep` (`float`): The step size for the x-axis.
* `ystep` (`float`): The step size for the y-axis.
* `nx` (`int`): The number of bins along the x-axis.
* `ny` (`int`): The number of bins along the y-axis.

**Returns:**

* `tuple`: A tuple containing two `numpy.ndarray` objects: `xedges` and `yedges`.

### <a id='get_skycoords'></a> `get_skycoords(skyimage, xstep, ystep, m_d_distance, verbose=False, radians=False)`

Calculates the sky coordinates (angles) corresponding to the sky image pixels.

**Parameters:**

* `skyimage` (`numpy.ndarray`): The sky image.
* `xstep` (`float`): The step size for the x-axis.
* `ystep` (`float`): The step size for the y-axis.
* `m_d_distance` (`float`): Mask-detector distance.
* `verbose` (`bool`, optional): If `True`, prints the X and Y bin ranges. Defaults to `False`.
* `radians` (`bool`, optional): If `True`, returns angles in radians. Defaults to `False` (returns degrees).

**Returns:**

* `tuple`: A tuple containing two `numpy.ndarray` objects: x-coordinates and y-coordinates.

### <a id='generate_bulk'></a> `generate_bulk(mask_shape, elxdim, elydim)`

Generates a bulk array representing the sensitive regions of the WFM detector plane.

**Parameters:**

* `mask_shape` (`tuple`): The desired shape of the bulk array (rows, columns).
* `elxdim` (`float`): Elemental dimension along X.
* `elydim` (`float`): Elemental dimension along Y.

**Returns:**

* `numpy.ndarray`: The generated bulk array where sensitive regions are marked with 1s and non-sensitive regions with 0s.

<a id='other_utils.py'></a>
## Other Utils

This module contains utility functions for general data manipulation.

### <a id='filter_source'></a> `filter_source(data, ra, dec, verbose=False)`

Filters a dataset based on Right Ascension (RA) and Declination (DEC) values.

**Parameters:**

* `data` (`numpy.ndarray` or `astropy.table.Table`): The input data, expected to have 'RA' and 'DEC' columns.
* `ra` (`float`): The Right Ascension value to filter by.
* `dec` (`float`): The Declination value to filter by.
* `verbose` (`bool`, optional): If `True`, prints the number of events selected. Defaults to `False`.

**Returns:**

* `numpy.ndarray` or `astropy.table.Table`: The filtered data containing only events matching the specified RA and DEC.
