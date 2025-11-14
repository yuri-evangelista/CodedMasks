import numpy as np
from regions import Regions, PixCoord, TextPixelRegion, RegionVisual
from regions import CircleSkyRegion
from astropy.coordinates import Angle
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm

def filter_source(data, ra, dec, verbose=False):
    mask = np.ones(len(data), dtype=bool)
    mask &= (np.isclose(data["RA"], ra) & np.isclose(data["DEC"], dec))
    filtered = data[mask]
    if verbose:
        print("Selected", len(filtered), "out of", len(data), "events")

    return filtered

def filter_energy(data, emin, emax, verbose=False):
    mask = np.ones(len(data), dtype=bool)
    mask &= ((data["ENERGY"] >= emin) & (data["ENERGY"] < emax))
    filtered = data[mask]
    if verbose:
        print("Selected", len(filtered), "out of", len(data), "events")

    return filtered

def plot_lb(img_f, img_b, wcs, fileout, spacing=15, regions=None):
    '''
    Plots composed image in l,b.
    Example usage:
    img_opt = np.pow(np.clip(img, a_min=3.5, a_max=35), 0.25)
    img_bkg = np.pow(np.clip(img, a_min=0, a_max=75), 0.25) #this to re-add some low significance noise
    regions = Regions.read("F:/CodedMasks/Simulations/galctr_rxte_sax_cat_cut_8kcts.reg", format='ds9')
    plot_lb(img_opt, img_bkg, wcs,  "GC_composed_sky_image_sign_preIROS.png", spacing=10, regions=regions))
    '''
    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs), figsize=(7, 7))

    # Plot the image data
    ax.imshow(img_f, cmap=cm.gnuplot2, interpolation='none', origin='lower')
    ax.imshow(img_b, cmap=cm.gnuplot2, interpolation='none', origin='lower', alpha=0.175)

    # Explicitly remove the Equatorial (RA, Dec) tick labels and axes labels.
    # This targets the *primary* WCS axes.
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticklabel_visible(False)
    ax.coords[0].set_axislabel('')
    ax.coords[1].set_axislabel('')

    # Also remove the major and minor tick *marks* themselves.
    ax.coords[0].set_ticks_visible(False)
    ax.coords[1].set_ticks_visible(False)

    # Add the Galactic Coordinate Grid (l, b)
    overlay = ax.get_coords_overlay('galactic')

    # Set the desired tick separation to 15 degrees
    #spacing = 15 * u.degree

    # Apply the spacing to Galactic Longitude (l) - overlay[0]
    overlay[0].set_major_formatter('d') # Use decimal degrees for formatting
    overlay[0].set_ticks(spacing=spacing * u.degree)
    overlay[0].set_axislabel('Galactic Longitude ($l$)', fontsize=16)
    overlay[0].set_ticklabel(fontsize=14)

    # Apply the spacing to Galactic Latitude (b) - overlay[1]
    overlay[1].set_major_formatter('d') # Use decimal degrees for formatting
    overlay[1].set_ticks(spacing=spacing * u.degree)
    overlay[1].set_axislabel('Galactic Latitude ($b$)', fontsize=16)
    overlay[1].set_ticklabel(fontsize=14)

    # The line below is also often the one that plots the RA/Dec grid lines,
    # which you should make sure is NOT present:
    # ax.grid(True, color='white', ls='solid') 
    # This overlay will now be the only one displaying ticks and labels.

    overlay.grid(color='gray', ls='dotted')
    overlay[0].set_axislabel('Galactic Longitude ($l$)')
    overlay[1].set_axislabel('Galactic Latitude ($b$)')

    if regions != None:
        visual = RegionVisual({'textangle': 60, 'color': 'white', 'fontsize':10})
        delta = PixCoord(x=0, y=30)
        for i, region in enumerate(regions):
            #patch = region.to_patch()
            #ax.add_patch(patch)
            text = region.meta['text']
            pixel_region = region.to_pixel(wcs)
            pixel_region.plot(ax=ax)
            reg = TextPixelRegion(center=pixel_region.center + delta, text=text, visual=visual)
            reg.plot(ax=ax)
           

    plt.savefig(fileout, format="png", bbox_inches="tight", dpi=300)

    plt.show()