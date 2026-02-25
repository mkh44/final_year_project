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
from scipy.constants import c
import matplotlib.pyplot as plt

# Define event and output location
event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

# Load FITS file
fits_file = glob.glob(os.path.join(output_loc, "*.fits"))[0]

# Compute Doppler shifts from spectral cube
hdul = fits.open(fits_file)
cube = hdul[4].data.astype(np.float32)
header = hdul[4].header


# GET DIMENSIONS
nw, ny, nt = cube.shape

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

spectrum = cube[:, y_plot, t_plot]

print(wavelength.min(), wavelength.max())
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

