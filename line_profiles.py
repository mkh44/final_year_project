#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from fit_iris_lines import fit_raster

# Define event and output location
event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

# Load FITS file
fits_file = glob.glob(os.path.join(output_loc, "*.fits"))[0]

# Compute Doppler shifts from spectral cube
hdul = fits.open(fits_file)


