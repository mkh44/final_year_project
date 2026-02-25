#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#

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

# Bulid wavelength axis
crval = header['CRVAL1']
cdelt = header['CDELT1']
crpix = header['CRPIX1']

nw = cube.shape[0]

wavelength = crval + (np.arange(nw) - (crpix -1)) * cdelt

nw, ny, nt = cube.shape

# Compute Doppler velocity using centroid
# Choose rest wavelength (Si IV for now)
rest_wavlen = 1393.27 # Angstrom
window = 2 # Angstrom
wl_min = rest_wavlen - window/2
wl_max = rest_wavlen + window/2

# Refine wavelength
wl_mask = (wavelength >= wl_min) & (wavelength <= wl_max)
wl_short = wavelength[wl_mask]

# Make Doppler map
y_plot = cube.shape[1] // 2  # middle row
t_plot = cube.shape[2] // 2  # middle time

spectrum = cube[:, y_plot, t_plot]
centroid = np.sum(wavelength * spectrum) / np.sum(spectrum)
velocity_pixel = c * (centroid - rest_wavlen) / rest_wavlen

# Only consider positive emission for centroid
spectrum_pos = np.clip(spectrum, a_min=0, a_max=None)  # set negatives to zero

spectrum_short = spectrum_pos[wl_mask]

# Compute Dopper Centroid
if np.sum(spectrum_short) > 0:
    centroid = np.sum(wl_short * spectrum_short) / np.sum(spectrum_short)
    velocity_pixel = c * (centroid - rest_wavlen) / rest_wavlen / 1000  # km/s
else:
    centroid = np.nan
    velocity_pixel = np.nan


# Plot spectral profiles
plt.figure(figsize=(8, 6))
plt.plot(wl_short, spectrum_short,
         label=f"y={y_plot}, t={t_plot}, v={velocity_pixel:.1f} km/s")
plt.xlabel("Wavelength (Å)")
plt.ylabel("Intensity")
plt.legend()

plt.show()