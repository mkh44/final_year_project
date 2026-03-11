#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#
import pdb
import glob
import os

import datetime as dt
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import twiny
from scipy.constants import c
import matplotlib.pyplot as plt
import matplotlib as mpl
import asdf
from get_quartiles import get_quartiles

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\line_profiles"
input_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"
si_file = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\IRIS_fitting_Si_IV_1403_20230503_072923.asdf"
cii_file = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\IRIS_fitting_C_II_1336_20230503_072923.asdf"
mg_file = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\IRIS_fitting_smooth_MgII_20230503_072923.asdf"

# Load FITS file
fits_file = glob.glob(os.path.join(input_loc, "iris_l2_20230503_072923_4204700135_raster_t000_r00000.fits"))[0]


lines = [5, 1, 9]
time_x = 12000
height = 100
y_lim = 210
x_lim = 250

hdul = fits.open(fits_file)
hdr = hdul[0].header

def get_wavelength(header, line_index):
    crval = header['CRVAL1']
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    n_wave = hdul[lines[line_index]].data.shape[2]
    wavelength = crval + (np.arange(n_wave) - (crpix - 1)) * cdelt
    return wavelength

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

def plot_line_profile(lines, time_x):
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


#plot_line_profile(lines, time_x)

# Defining line titles
def get_int_map(file):
    with asdf.open(file) as af:
        keys = af.tree.keys()
        # Si IV / C II files
        if 'q_int_map' in keys:
            return af.tree['q_int_map']

        # Mg II files
        elif 'mgii_k_integ_int' in keys:
            return af.tree['mgii_k_integ_int']
        elif 'mgii_h_integ_int' in keys:
            return af.tree['mgii_h_integ_int']
        else:
            raise KeyError("No intensity map found in ASDF file")


# Plotting Intensity quartiles reference
def plot_iris_sns_quartile_fits(si_title, mg_title, event, main_header):

    #Get int maps
    si_q_int_map = get_int_map(si_file)
    #cii_q_int_map = get_int_map(cii_file)
    mg_q_int_map = get_int_map(mg_file)

    cadence = main_header['STEPT_AV']

    fig, ax = plt.subplots(3, 1, sharex=True)

    alpha = 1

 # Si IV
    si_t_array = np.arange(si_q_int_map.data.shape[1]) * cadence
    si_slit_pos = si_q_int_map.meta['crval2'] + si_q_int_map.meta['cdelt2'] * (
        np.arange(si_q_int_map.data.shape[0]) - si_q_int_map.meta['crpix2'])

    si_upr_bnd = np.nanpercentile(si_q_int_map.data, 100 - alpha)

    im0 = ax[0].imshow(si_q_int_map.data, origin='lower', cmap='Reds_r', aspect='auto',
        extent=[si_t_array.min(), si_t_array.max(), si_slit_pos.min(), si_slit_pos.max()],
        norm=colors.Normalize(vmin=0, vmax=si_upr_bnd))

    ax[0].set_ylabel(' ')

# second Si axis for label
    ax_0 = ax[0].twinx()
    ax_0.set_ylabel(si_title)
    ax_0.set_yticks([])


    # # C ii
    # cii_t_array = np.arange(cii_q_int_map.data.shape[1]) * cadence
    # cii_slit_pos = cii_q_int_map.meta['crval2'] + cii_q_int_map.meta['cdelt2'] * (
    #         np.arange(cii_q_int_map.data.shape[0]) - cii_q_int_map.meta['crpix2'])
    #
    # cii_upr_bnd = np.nanpercentile(cii_q_int_map.data, 100 - alpha)
    #
    # im1 = ax[1].imshow(cii_q_int_map.data, origin='lower', cmap='Reds_r', aspect='auto',
    #     extent=[cii_t_array.min(), cii_t_array.max(), cii_slit_pos.min(),
    #     cii_slit_pos.max()], norm=colors.Normalize(vmin=0, vmax=cii_upr_bnd))
    #
    # ax[1].set_ylabel("Solar Y")
    # ax[1].set_title(' ')

    # second cii axis for label
    # ax_1 = ax[1].twinx()
    # ax_1.set_ylabel(cii_title)
    # ax_1.set_yticks([])

# Mg ii
    mg_t_array = np.arange(mg_q_int_map.data.shape[1]) * cadence
    mg_slit_pos = mg_q_int_map.meta['crval2'] + mg_q_int_map.meta['cdelt2'] * (
            np.arange(mg_q_int_map.data.shape[0]) - mg_q_int_map.meta['crpix2'])

    mg_upr_bnd = np.nanpercentile(mg_q_int_map.data, 100 - alpha)

    im2 = ax[2].imshow(mg_q_int_map.data, origin='lower', cmap='Reds_r', aspect='auto',
        extent=[mg_t_array.min(), mg_t_array.max(), mg_slit_pos.min(), mg_slit_pos.max()],
        norm=colors.Normalize(vmin=0, vmax=mg_upr_bnd))

    ax[2].set_ylabel(" ")
    ax[2].set_xlabel("Time (s)")
    ax[2].set_title(' ')

    # second mg axis for label
    ax_2 = ax[2].twinx()
    ax_2.set_ylabel(mg_title)
    ax_2.set_yticks([])

    # Colorbar
    fig.colorbar(im2, ax=ax, label="Integrated Intensity")

    plt.suptitle(f"Intensity Quartiles")

    save_path = os.path.join(output_loc, f"quartile_maps_{event}.png")
    plt.savefig(save_path, bbox_inches="tight")

    plt.show()

plot_iris_sns_quartile_fits(si_title, mg_title, event, hdr)

