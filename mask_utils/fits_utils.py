import numpy as np
import astropy.io.fits as pyfits 
from astropy.table import Table, Column
import ezdxf
from ezdxf.addons import r12writer
from scipy.ndimage import label
import pprint 

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

def write_mask_fits(FITSFILE, MASK, RMATRIX, BULK, PROPS):
    r"""
    Writes a FITS file with the following extensions:
        1 - OR_MASK
        2 - MASK
        3 - RMATRIX
        4 - SENS    
    """
    pprint.pp(PROPS, sort_dicts=False)

    row_count = PROPS['ELXN'] * PROPS['ELYN']
    print("\nTotal row count", row_count)

    #Generate arrays
    el = np.arange(row_count)
    idx = np.unravel_index(el, MASK.shape)
    X = idx[0] *  PROPS['ELXDIM'] +  PROPS['ELXDIM']/2 -  PROPS['MXDIM']/2 
    Y = idx[1] *  PROPS['ELYDIM'] +  PROPS['ELYDIM']/2 -  PROPS['MYDIM']/2 
    MASK_VAL = MASK[idx].astype('int')
    RMATRIX_VAL = RMATRIX[idx]
    BULK_VAL =  BULK[idx]

    #Generating primary HDU
    primary_hdu = pyfits.PrimaryHDU()

    #Generating columns
    x_col = pyfits.Column(name='X', format='E', array=X)
    y_col = pyfits.Column(name='Y', format='E', array=Y)
    mask_val_col = pyfits.Column(name='VAL', format='E', array=MASK_VAL)
    rmatrix_val_col = pyfits.Column(name='VAL', format='E', array=RMATRIX_VAL)
    bulk_val_col = pyfits.Column(name='VAL', format='E', array=BULK_VAL)

    #Grouping columns in hdus
    mask_cols = pyfits.ColDefs([x_col, y_col, mask_val_col])
    rmatrix_cols = pyfits.ColDefs([x_col, y_col, rmatrix_val_col])
    bulk_cols = pyfits.ColDefs([x_col, y_col, bulk_val_col])

    #Generating HDUs from columns
    or_mask_hdu = pyfits.BinTableHDU.from_columns(mask_cols, name="OR_MASK")
    mask_hdu = pyfits.BinTableHDU.from_columns(mask_cols, name="MASK")
    rmatrix_hdu = pyfits.BinTableHDU.from_columns(rmatrix_cols, name="RMATRIX")
    bulk_hdu = pyfits.BinTableHDU.from_columns(bulk_cols, name="SENS")

    #Updating headers
    or_mask_hdu.header.update(PROPS)
    mask_hdu.header.update(PROPS)
    rmatrix_hdu.header.update(PROPS)
    bulk_hdu.header.update(PROPS)

    #Putting together HDUs in HDU list
    out_hdu_list = pyfits.HDUList([primary_hdu, or_mask_hdu, mask_hdu, rmatrix_hdu, bulk_hdu])

    #Writing to file
    out_hdu_list.writeto(FITSFILE, overwrite=True)

def fits_mask_to_dxf(fitsin, dxfout):
    r"""
    Write a DXF file of the mask with all the open elements as polyline
    """
    def get_label_index(arr, _label):
        # Returns the first and last indexes of a labeled region
        # e.g. 
        # (x1, y1), (x2, y2) = get_label_index(l, 2)
        idx = np.where(arr == _label)
        return (idx[0][0], idx[1][0]), (idx[0][-1], idx[1][-1])

    #reading mask file as an image
    mask, header = read_mask_bulk(fitsin, 'MASK', header_out=True, verbose=False)
    ELXDIM = header['ELXDIM']
    ELYDIM = header['ELYDIM']
    MXDIM  = header['MXDIM']
    MYDIM  = header['MYDIM']

    #Finding open regions
    l, num_features = label(mask)

    #Initializing DXF file
    doc = ezdxf.new(dxfversion='R2010')
    msp = doc.modelspace()

    for feature in range(1, num_features):
        #Getting open region corner indexes
        (i1, j1), (i2, j2) = get_label_index(l, feature)
        #Transforming indexes in spatial coordinates
        x1, y1 = i1 * ELXDIM -  MXDIM/2,  j1 * ELYDIM -  MYDIM/2
        x2, y2 = (i2+1) * ELXDIM -  MXDIM/2, (j2+1) * ELYDIM - MYDIM/2

        #Drawing rectangle as closed polyline
        msp.add_lwpolyline([(x1, y1),(x2, y1),(x2, y2),(x1, y2),(x1, y1)],close=True)

    #Writing DXF file
    doc.saveas(dxfout)
    print(dxfout, 'file created')