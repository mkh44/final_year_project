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
from fit_iris_lines import fit_raster
from astropy import units as u

# Define event and output location
event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

# Load FITS file
fits_file = glob.glob(os.path.join(output_loc, "*.fits"))

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

cube.shape = (nw, ny, nt)

# Compute Doppler velocity using centroid
# Choose rest wavelength (Si IV for now)
rest_wavlen = 1393.27 # Angstrom

# Compute centroid


