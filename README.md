# CodedMasks
A collection of (quick&dirty) utilities for coded mask images.
Comes with no warranty (and with a punk fashioned code...)

## Code Utils
### Tools to generate URA, MURA and other codes
.\mask_utils\code_utils.py
- next_prime
- ura_mura
- bura
- bura33
## Fits Utils
### Tools to read and write WFM fits
.\mask_utils\fits_utils.py
- read_mask_bulk
- read_mask_bulk
- write_mask_fits
- fits_mask_to_dxf
## Image Utils
### Tools to manipulate images
.\mask_utils\image_utils.py
- shift
- fshift
- erosion
- pad_array
- float_gcd
- upscale
## Imaging Utils
### Tools to perform decoding operations and calculate system properties (effective area, solid angle, etc.)
.\mask_utils\imaging_utils.py
- get_openfraction
- get_angular_res
- get_coding_power
- open_fraction_vs_off_axis
- eff_area_vs_off_axis
- omega_plate_offaxis
- solid_angle
- get_detimage
- decode
- decode_var
- get_skysign
- get_detimage_edges
- get_skycoords
- generate_bulk
## Other Utils
### Miscellaneous tools
.\mask_utils\other_utils.py
- filter_source
