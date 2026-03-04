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
from matplotlib.pyplot import twiny
from scipy.constants import c
import matplotlib.pyplot as plt
from get_quartiles import get_quartiles
from fit_iris_lines import get_line_references
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

# Load FITS file
fits_file = glob.glob(os.path.join(output_loc, "iris_l2_20230503_072923_4204700135_raster_t000_r00000.fits"))[0]

lines = [5, 1, 9]
time_x = 7000
height = 100
y_lim = 210
x_lim = 200

hdul = fits.open(fits_file)

def line_profile(lines, time_x):
    hdr = hdul[0].header

    # Get line names
    si_title = hdr['TDESC' + str(lines[0])]
    cii_title = hdr['TDESC' + str(lines[1])]
    mg_title = hdr['TDESC' + str(lines[2])]

    # Get Si IV properties
    si_header = hdul[lines[0]].header
    si_crval = si_header['CRVAL1']
    si_cdelt = si_header['CDELT1']
    si_crpix = si_header['CRPIX1']
    si_wave = hdul[lines[0]].data.shape[2]
    si_wavelength = si_crval + (np.arange(si_wave) - (si_crpix - 1)) * si_cdelt
    si_data_arr = hdul[lines[0]].data

    # Get Cii properties
    cii_header = hdul[lines[1]].header
    cii_crval = cii_header['CRVAL1']
    cii_cdelt = cii_header['CDELT1']
    cii_crpix = cii_header['CRPIX1']
    cii_wave = hdul[lines[1]].data.shape[2]
    cii_wavelength = cii_crval + (np.arange(cii_wave) - (cii_crpix - 1)) * cii_cdelt
    cii_data_arr = hdul[lines[1]].data

    # Get mg properties
    mg_header = hdul[lines[2]].header
    mg_crval = mg_header['CRVAL1']
    mg_cdelt = mg_header['CDELT1']
    mg_crpix = mg_header['CRPIX1']
    mg_wave = hdul[lines[2]].data.shape[2]
    mg_wavelength = mg_crval + (np.arange(mg_wave) - (mg_crpix - 1)) * mg_cdelt
    mg_data_arr = hdul[lines[2]].data

    # Doppler velocities
    lambda_si = 1402.8
    si_v_dopp = ((si_wavelength - lambda_si) / lambda_si) * (c / 1e3)

    lambda_cii = 1335.7
    cii_v_dopp = ((cii_wavelength - lambda_cii) / lambda_cii) * (c / 1e3)

    lambda_mg = 2795.5
    mg_v_dopp = ((mg_wavelength - lambda_mg) / lambda_mg) * (c / 1e3)

    # Plotting
    plt.figure(figsize=(8, 6))
    fig, ax = plt.subplots(3, 1)
    ax[0].plot(si_v_dopp, si_data_arr[time_x, height], color='k')
    ax[0].set_xlabel(' ')
    ax[0].set_ylim(0, y_lim)
    ax0 = ax[0].twinx()
    ax0.set_yticks([])
    ax0.set_yticklabels([])
    ax0.set_ylabel(f'{si_title}')
    si_vmax = np.max(np.abs(si_v_dopp))
    ax[0].set_xlim(-x_lim, x_lim)
    ax[0].xaxis.set_minor_locator(MultipleLocator(10))

    ax[1].plot(cii_v_dopp, cii_data_arr[time_x, height], color='k')
    ax[1].set_xlabel(' ')
    ax[1].set_ylabel("Intensity")
    ax[1].set_ylim(0, y_lim)
    ax1 = ax[1].twinx()
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_ylabel(f'{cii_title}')
    cii_vmax = np.max(np.abs(cii_v_dopp))
    ax[1].set_xlim(-x_lim, x_lim)
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))

    ax[2].plot(mg_v_dopp, mg_data_arr[time_x, height], color='k')
    ax[2].set_ylabel(' ')
    ax[2].set_xlabel('Doppler Velocity (km/s)')
    ax[2].set_ylim(0, y_lim)
    ax2 = ax[2].twinx()
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax2.set_ylabel(f'{mg_title}')
    mg_vmax = np.max(np.abs(mg_v_dopp))
    ax[2].set_xlim(-x_lim, x_lim)
    ax[2].xaxis.set_minor_locator(MultipleLocator(10))

    plt.suptitle(f'Time: {time_x} s')
    plt.show()


line_profile(lines, time_x)


