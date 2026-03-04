##!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Routine originally written by David Long for getting Mg II profile properties using quartiles
## Modified to also work for other lines

import warnings
warnings.filterwarnings('ignore')

import numpy as np
from tqdm import tqdm
from astropy.io import fits
from iris_fitting import extract_irisL2data
from astropy.wcs import WCS
import pdb

# Calculate the quartiles of a spectral line profile
def estimate_quartiles(wavelength, profile):
    """
    Estimate the quartiles of the spectral line profile.
    
    Parameters:
    profile : array-like
        The intensity profile of the spectral line.
    wavelength : array-like
        The corresponding wavelengths for the profile.
    
    Returns:
    tuple : (q1_wvl, q2_wvl, q3_wvl)
        The wavelengths corresponding to the first, second (median), and third quartiles.
    """

    profile[profile < 0] = 0  # Set negative values to zero

    if np.count_nonzero(profile)>0:
        cs = np.nancumsum(profile)
        cdf = cs/cs[-1]
        integ_int = cs[-1]

        xrng = np.linspace(wavelength[0],wavelength[-1],2000)
        yrng = np.interp(xrng,wavelength,cdf)
        q1_wvl_interp = xrng[np.argmin(np.abs(yrng - 0.25))]
        q2_wvl_interp = xrng[np.argmin(np.abs(yrng - 0.5))]
        q3_wvl_interp = xrng[np.argmin(np.abs(yrng - 0.75))]

    else:
        integ_int = np.nan
        q1_wvl_interp = np.nan
        q2_wvl_interp = np.nan
        q3_wvl_interp = np.nan
    
    return q1_wvl_interp, q2_wvl_interp, q3_wvl_interp, integ_int


# Main routine to get the quartiles of the spectral line
def get_quartiles(file, line_wvl, iris_window, fulldisk=False):

    A_to_nm = 10  # convert wavelength to nm

    match line_wvl:
        case 'Mg II k':
            est_wvl = 279.63509493
            spectral_window = 0.1
        case 'Mg II h':
            est_wvl = 280.35297192
            spectral_window = 0.1
        case 'C II 1334':
            est_wvl = 133.453
            spectral_window = 0.03
        case 'C II 1335':
            est_wvl = 133.571
            spectral_window = 0.03
        case 'Si IV 1393':
            est_wvl = 139.376
            spectral_window = 0.03
        case 'Si IV 1394':
            est_wvl = 139.376
            spectral_window = 0.03
        case 'Si IV 1403':
            est_wvl = 140.277
            spectral_window = 0.03

    if fulldisk:
        # Open the IRIS FITS file
        sp = fits.open(file, memmap=False)

        ind = 0
        header = sp[ind].header
        data = sp[ind].data.T
        # Get the wavelength
        nwave = data.shape[2]
        wavelength = (np.array(range(header['NAXIS3'])) * header['CDELT3'] + 
                      (header['CRVAL3']-((header['NAXIS3']/2)*header['CDELT3'])))/A_to_nm
    else:
        sp = fits.open(file, memmap=False)
        sub_header = sp[1].header

        if sub_header['CDELT3'] == 0:
            ind = np.where(extract_irisL2data.show_lines(file)==iris_window)[0][0]
            extension = ind+1
            data = sp[extension].data
            header = sp[extension].header
            nwave = data.shape[2]
        # Get the wavelength
            wavelength = (np.array(range(header['NAXIS1'])) * header['CDELT1'] + header['CRVAL1'])/A_to_nm
        else:
            rast = extract_irisL2data.load(file,window_info=[iris_window],verbose=False)
            loc = np.where(extract_irisL2data.show_lines(file)==iris_window)[0][0]
            extension = loc+1
            head = extract_irisL2data.only_header(file,extension=extension)
            wcs = WCS(head)
            data = rast.raster[iris_window].data
            m_to_nm = 1e9  # convert wavelength to nm
            nwave = data.shape[2]
            wavelength = wcs.all_pix2world(np.arange(nwave), [0.], [0.], 0)[0] * m_to_nm
    
    k_min = est_wvl - spectral_window
    k_max = est_wvl + spectral_window

    wvl_crop = wavelength[(wavelength >= k_min) & (wavelength <= k_max)]
    datacube_crop = data[:,:,(wavelength >= k_min) & (wavelength <= k_max)]

    y_size = datacube_crop.shape[0]
    x_size = datacube_crop.shape[1]

    profile = [(datacube_crop[y, x, :]) for y in range(y_size) for x in range(x_size)]

    results = []
    for ind in tqdm(np.arange(0, len(profile)), desc="Calculating quartiles"):
        res_array = estimate_quartiles(wvl_crop, profile[ind])
        results.append(res_array)

    # Reshape the results back into the (y, x, 4) shape (numpy style)
    res = np.array([result for result in results], dtype='object').reshape(y_size, x_size, 4)
    quartiles = res.astype(np.float32)

    return quartiles*A_to_nm