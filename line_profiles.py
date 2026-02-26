#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#
import pdb
import glob
import os
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

# print('Window. Name      : wave start - wave end\n')
# for i in range(header['NWIN']):
#     win = str(i + 1)
#     print('{0}. {1:15}: {2:.2f} - {3:.2f} Å'
#           ''.format(win, header['TDESC' + win], header['TWMIN' + win], header['TWMAX' + win]))
header = hdul[5].header
#print(hdul[0].header)

crval = header['CRVAL1']
cdelt = header['CDELT1']
crpix = header['CRPIX1']

start_t = hdul[0].header['STARTOBS']
print(start_t)

#wcs = WCS(hdul[5].header)
nwave = hdul[5].data.shape[2]
print(nwave)
wavelength = crval + (np.arange(nwave) -(crpix-1)) * cdelt
print(wavelength)
data_arr = hdul[5].data
print(data_arr.shape)
# wavelen to velocity (dopp shift strong 1st line)
#only interested in >0 so set anything <0 equal to zero
#pdb.set_trace()
# Plot spectral profiles
plt.figure(figsize=(8, 6))
plt.plot(wavelength, data_arr[11850, 100])
plt.ylim(0, 100)
plt.xlabel("Wavelength (Å)")
plt.ylabel("Intensity")
plt.legend()
plt.show()

