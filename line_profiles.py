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
window = 1
wl_min = rest_wavlen - window/2
wl_max = rest_wavlen + window/2

# Refine wavelength
wl_mask = (wavelength >= wl_min) & (wavelength <= wl_max)
wl_short = wavelength[wl_mask]

# Compute centroid
velocity = np.zeros((cube.shape[1], cube.shape[2]))

# Make Doppler map
for y in range(cube.shape[1]):
    for t in range(cube.shape[2]):

        spectrum = cube[:, y, t]

        if np.all(np.isnan(spectrum)):
            velocity[y, t] = np.nan
            continue

        centroid = np.sum(wavelength * spectrum) / np.sum(spectrum)
        velocity[y, t] = c * (centroid - rest_wavlen) / rest_wavlen


# Find fastest redshift pixels
flat_indices = np.argsort(velocity.flatten())[::-1]

# Get top N redshift pixels
N = 2

indices = np.unravel_index(flat_indices[:N], velocity.shape)

y_indices = indices[0]
t_indices = indices[1]

# Plot spectral profiles
for i, (y, t) in enumerate((y_indices, t_indices)):
    spectrum = cube[:, y, t]
    plt.figure(figsize=(8, 6))
    plt.plot(wl_short, spectrum[wl_mask],
             label=f"y={y}, t={t}, v={velocity[y, t]:.1f} km/s")
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Intensity")
    plt.legend()

    plt.show()