#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#
import pdb
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from scipy.constants import c
import matplotlib.pyplot as plt

# Define event and output location
event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

# Load FITS file
fits_file = glob.glob(os.path.join(output_loc, "iris_l2_20230503_072923_4204700135_raster_t000_r00000.fits"))[0]

# Compute Doppler shifts from spectral cube
hdul = fits.open(fits_file)
header = hdul[0].header
print('Window. Name      : wave start - wave end\n')
for i in range(header['NWIN']):
    win = str(i + 1)
    print('{0}. {1:15}: {2:.2f} - {3:.2f} Å'
          ''.format(win, header['TDESC' + win], header['TWMIN' + win], header['TWMAX' + win]))

data = hdul[5].data
wcs = WCS(hdul[1].header)

# GET DIMENSIONS
nw, ny, nt = data.shape
print(data.shape)

# Solar y axis
crval_y = header['CRVAL2']
cdelt_y = header['CDELT2']
crpix_y = header['CRPIX2']

y_arcsec = crval_y + (np.arange(ny) - (crpix_y - 1)) * cdelt_y

# Time axis
crval_t = header['CRVAL3']
cdelt_t = header['CDELT3']
crpix_t = header['CRPIX3']

time_sec = crval_t + (np.arange(nt) - (crpix_t - 1)) * cdelt_t
# Bulid wavelength axis
crval = header['CRVAL1']
cdelt = header['CDELT1']
crpix = header['CRPIX1']

y_chosen = 240
t_chosen = 7000

wavelength = crval + (np.arange(nw) - (crpix -1)) * cdelt
# Compute Doppler velocity using centroid
# Choose rest wavelength (Si IV for now)
rest_wavlen = 1403 # Angstrom

# Make Doppler map
y_plot = 100
t_plot = 50

spectrum = data[:, y_plot, t_plot]

#pdb.set_trace()
# Plot spectral profiles
plt.figure(figsize=(8, 6))
plt.plot(wavelength, spectrum,
         label=f"y={y_arcsec[y_plot]:.1f}\"  "
               f"t={time_sec[t_plot]:.1f}s  ")
plt.xlim(1402, 1404)
plt.xlabel("Wavelength (Å)")
plt.ylabel("Intensity")
plt.legend()
plt.show()

