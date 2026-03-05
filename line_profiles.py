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

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\line_profiles"
input_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

# Load FITS file
fits_file = glob.glob(os.path.join(input_loc, "iris_l2_20230503_072923_4204700135_raster_t000_r00000.fits"))[0]

lines = [5, 1, 9]
time_x = 7000
height = 100
y_lim = 75
x_lim = 150

hdul = fits.open(fits_file)
hdr = hdul[0].header

def get_wavelength(header, line_index):
    crval = header['CRVAL1']
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    n_wave = hdul[lines[line_index]].data.shape[2]
    wavelength = crval + (np.arange(n_wave) - (crpix - 1)) * cdelt
    return wavelength

def line_profile(lines, time_x):

    # Get line names
    si_title = hdr['TDESC' + str(lines[0])]
    cii_title = hdr['TDESC' + str(lines[1])]
    mg_title = hdr['TDESC' + str(lines[2])]

    # Get Si IV properties
    si_header = hdul[lines[0]].header
    si_wavelength = get_wavelength(si_header,0)
    si_data_arr = hdul[lines[0]].data

    # Get Cii properties
    cii_header = hdul[lines[1]].header
    cii_wavelength = get_wavelength(cii_header, 1)
    cii_data_arr = hdul[lines[1]].data

    # Get mg properties
    mg_header = hdul[lines[2]].header
    mg_wavelength = get_wavelength(mg_header, 2)
    mg_data_arr = hdul[lines[2]].data

    # Doppler velocities
    lambda_si = 1402.8
    si_v_dopp = ((si_wavelength - lambda_si) / lambda_si) * (c / 1e3)

    lambda_cii = 1335.7
    cii_v_dopp = ((cii_wavelength - lambda_cii) / lambda_cii) * (c / 1e3)

    lambda_mg = 2795.5
    mg_v_dopp = ((mg_wavelength - lambda_mg) / lambda_mg) * (c / 1e3)

# Plotting Doppler velocity

    fig, ax = plt.subplots(3, 1, sharex = True)


    # Masks for colour difference on plot
    si_redshift = si_v_dopp >= 0
    si_blueshift = si_v_dopp <= 0

    cii_redshift = cii_v_dopp >= 0
    cii_blueshift = cii_v_dopp <= 0

    mg_redshift = mg_v_dopp >= 0
    mg_blueshift = mg_v_dopp <= 0

    #Si IV 1403 plot
    ax[0].plot(si_v_dopp, si_data_arr[time_x, height], color='k')
    ax[0].plot(si_v_dopp[si_redshift], si_data_arr[time_x, height][si_redshift], color='red')
    ax[0].plot(si_v_dopp[si_blueshift], si_data_arr[time_x, height][si_blueshift], color='blue')

    # Set Si IV axis labels abd tickmarks
    ax[0].set_xlabel(' ')
    ax[0].set_ylim(0, y_lim)
    ax[0].set_xlim(-x_lim, x_lim)
    ax[0].xaxis.set_minor_locator(MultipleLocator(10))


    # Second Si IV axis for titles
    ax0 = ax[0].twinx()
    ax0.set_yticks([])
    ax0.set_ylim(0, y_lim)
    ax0.set_yticklabels([])
    ax0.set_ylabel(f'{si_title}')

    plt.title(f'Time: {time_x} s', loc='right')

    # Cii plot
    ax[1].plot(cii_v_dopp, cii_data_arr[time_x, height], color='k')
    ax[1].plot(cii_v_dopp[cii_redshift], cii_data_arr[time_x, height][cii_redshift], color='red')
    ax[1].plot(cii_v_dopp[cii_blueshift], cii_data_arr[time_x, height][cii_blueshift], color='blue')

    # Set Cii axis limits and labels
    ax[1].set_xlim(-x_lim, x_lim)
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))
    ax[1].set_ylim(0, y_lim)
    ax[1].set_xlabel(' ')
    ax[1].set_ylabel("Intensity")

    # Seconds Cii axis for titles
    ax1 = ax[1].twinx()
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_ylabel(f'{cii_title}')

    # Mg plot
    ax[2].plot(mg_v_dopp, mg_data_arr[time_x, height], color='k')
    ax[2].plot(mg_v_dopp[mg_redshift], mg_data_arr[time_x, height][mg_redshift], color='red')
    ax[2].plot(mg_v_dopp[mg_blueshift], mg_data_arr[time_x, height][mg_blueshift], color='blue')

    # Set Mg axis limits and labels
    ax[2].set_ylabel(' ')
    ax[2].set_xlabel('Doppler Velocity (km/s)')
    ax[2].set_ylim(0, y_lim)
    ax[2].set_xlim(-x_lim, x_lim)
    ax[2].xaxis.set_minor_locator(MultipleLocator(10))

    # Second Mg axis for titles
    ax2 = ax[2].twinx()
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax2.set_ylabel(f'{mg_title}')

    # Dotted line at x=0
    ax[0].axvline(0, color='k', linestyle='dashed', linewidth=1)
    ax[1].axvline(0, color='k', linestyle='dashed', linewidth=1)
    ax[2].axvline(0, color='k', linestyle='dashed', linewidth=1)

    # Saving plot and displaying
    save_path = os.path.join(output_loc, f"doppler_profiles_{time_x}s.png")
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# Plotting wavelength
    fig, ax = plt.subplots(3, 1)
    ax[0].plot(si_wavelength, si_data_arr[time_x, height], color='k')
    ax[0].set_xlabel(' ')
    ax[0].set_ylim(0, y_lim)
    ax0 = ax[0].twinx()
    ax0.set_yticks([])
    ax0.set_yticklabels([])
    ax0.set_ylabel(f'{si_title}')
    ax[0].xaxis.set_minor_locator(MultipleLocator(10))
    plt.title(f'Time: {time_x} s', loc='right')

    ax[1].plot(cii_wavelength, cii_data_arr[time_x, height], color='k')
    ax[1].set_xlabel(' ')
    ax[1].set_ylabel("Intensity")
    ax[1].set_ylim(0, y_lim)
    ax1 = ax[1].twinx()
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_ylabel(f'{cii_title}')
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))

    ax[2].plot(mg_wavelength, mg_data_arr[time_x, height], color='k')
    ax[2].set_ylabel(' ')
    ax[2].set_xlabel('Wavelength (Å)')
    ax[2].set_ylim(0, y_lim)
    ax2 = ax[2].twinx()
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax2.set_ylabel(f'{mg_title}')
    ax[2].xaxis.set_minor_locator(MultipleLocator(10))


    save_path = os.path.join(output_loc, f"wavelength_profiles_{time_x}s.png")
    #plt.savefig(save_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)

line_profile(lines, time_x)


