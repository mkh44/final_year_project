#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#
import pdb
import glob
import os
from matplotlib import ticker, gridspec
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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, FixedLocator)


event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\line_profiles"
input_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"
si_file = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\IRIS_fitting_Si_IV_1403_20230503_072923.asdf"
cii_file = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\IRIS_fitting_C_II_1334_20230503_072923.asdf"
mg_file = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\IRIS_fitting_smooth_MgII_20230503_072923.asdf"

# Load FITS file
fits_file = glob.glob(os.path.join(input_loc, "iris_l2_20230503_072923_4204700135_raster_t000_r00000.fits"))[0]


lines = [5, 1, 9]
time_seconds = [9795, 9820, 9850, 9900]
position_solar_y = 266

# [6900, 6950, 7000, 7050]
# [9600, 9700, 9750, 9800]
# [11500, 11750, 11900, 12050]



x_lim = 450
zoom = 50
fs = 16 # font size

def time_to_index(time_array, target_time):
    return np.argmin(np.abs(time_array - target_time))

def get_data_arr(line):
    return hdul[lines[line]].data

def get_wavelength(header, line_index):
    crval = header['CRVAL1']
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    n_wave = hdul[lines[line_index]].data.shape[2]
    wavelength = crval + (np.arange(n_wave) - (crpix - 1)) * cdelt
    return wavelength


# Defining line titles
def get_int_map(file):
    with asdf.open(file) as af:
        keys = af.tree.keys()

        # Si IV / C II
        if 'q_int_map' in keys:
            return af.tree['q_int_map']

        # Mg II files
        elif 'mgii_k_integ_int' in keys:
            return af.tree['mgii_k_integ_int']
        elif 'mgii_h_integ_int' in keys:
            return af.tree['mgii_h_integ_int']
        else:
            raise KeyError("No intensity map found in ASDF file")

def get_line_profiles(line, expected_wl):
# Getting pixel array
    data_arr = get_data_arr(line)
    title = hdul[lines[line]].header

# Getting time and height in index
    cadence = hdr['STEPT_AV']
    time_array = np.arange(data_arr.shape[0]) * cadence

    wavelength = get_wavelength(title, line)

    v_dopp = ((wavelength - expected_wl) / expected_wl) * (c / 1e3)

    # spatial index
    slit_pos = None
    if line == 0:
        slit_pos = si_slit_pos
    elif line == 1:
        slit_pos = cii_slit_pos
    elif line == 2:
        slit_pos = mg_slit_pos

    pos_idx = np.argmin(np.abs(slit_pos - position_solar_y))

# time index
    t_indices = [time_to_index(time_array, t) for t in time_seconds]

    profiles = []
    for t_idx in t_indices:
        profiles.append(data_arr[t_idx, pos_idx, :])

    return v_dopp, profiles

# # Defining limits
#
#     # y (intensity)
#     y_max = max(si_data_arr[time_idx, pos_idx].max(),
#                 cii_data_arr[time_idx, pos_idx].max(),
#                 mg_data_arr[time_idx, pos_idx].max())
#     y_lim = y_max + y_max*0.05
#
#
# # Plotting Doppler velocity
#     fig, ax = plt.subplots(3, 1, sharex = True, figsize=(6, 6))
#
#
#     # Masks for colour difference on plot
#     si_redshift = si_v_dopp >= 0
#     si_blueshift = si_v_dopp <= 0
#
#     cii_redshift = cii_v_dopp >= 0
#     cii_blueshift = cii_v_dopp <= 0
#
#     mg_redshift = mg_v_dopp >= 0
#     mg_blueshift = mg_v_dopp <= 0
#
#     #Si IV 1403 plot
#     ax[0].plot(si_v_dopp, si_data_arr[time_idx, pos_idx], color='k')
#     ax[0].plot(si_v_dopp[si_redshift], si_data_arr[time_idx, pos_idx][si_redshift], color='red')
#     ax[0].plot(si_v_dopp[si_blueshift], si_data_arr[time_idx, pos_idx][si_blueshift], color='blue')
#
#     # Set Si IV axis labels abd tickmarks
#     ax[0].set_xlabel(' ')
#     ax[0].set_ylim(0, y_lim)
#     ax[0].set_xlim(-x_lim, x_lim)
#     ax[0].xaxis.set_minor_locator(MultipleLocator(10))
#
#     # Second Si IV axis for titles
#     ax0 = ax[0].twinx()
#     ax0.set_yticks([])
#     ax0.set_ylim(0, y_lim)
#     ax0.set_yticklabels([])
#     ax0.set_ylabel(f'{si_title}', fontsize=fs)
#
#     plt.title(f'Time: {actual_time:.1f} s, Height: {position_solar_y} arcsec', loc='right', fontsize=fs)
#
#     # Cii plot
#     ax[1].plot(cii_v_dopp, cii_data_arr[time_idx, pos_idx], color='k')
#     ax[1].plot(cii_v_dopp[cii_redshift], cii_data_arr[time_idx, pos_idx][cii_redshift], color='red')
#     ax[1].plot(cii_v_dopp[cii_blueshift], cii_data_arr[time_idx, pos_idx][cii_blueshift], color='blue')
#
#     # Set Cii axis limits and labels
#     ax[1].set_xlim(-x_lim, x_lim)
#     ax[1].xaxis.set_minor_locator(MultipleLocator(10))
#     ax[1].set_ylim(0, y_lim)
#     ax[1].set_xlabel(' ')
#     ax[1].set_ylabel("Intensity", fontsize=fs)
#
#     # Seconds Cii axis for titles
#     ax1 = ax[1].twinx()
#     ax1.set_yticks([])
#     ax1.set_yticklabels([])
#     ax1.set_ylabel(f'{cii_title}', fontsize=fs)
#
#     # Mg plot
#     ax[2].plot(mg_v_dopp, mg_data_arr[time_idx, pos_idx], color='k')
#     ax[2].plot(mg_v_dopp[mg_redshift], mg_data_arr[time_idx, pos_idx][mg_redshift], color='red')
#     ax[2].plot(mg_v_dopp[mg_blueshift], mg_data_arr[time_idx, pos_idx][mg_blueshift], color='blue')
#
#     # Set Mg axis limits and labels
#     ax[2].set_ylabel(' ')
#     ax[2].set_xlabel('Doppler Velocity (km/s)', fontsize=fs)
#     ax[2].set_ylim(0, y_lim)
#     ax[2].set_xlim(-x_lim, x_lim)
#     ax[2].xaxis.set_minor_locator(MultipleLocator(10))
#
#     # Second Mg axis for titles
#     ax2 = ax[2].twinx()
#     ax2.set_yticks([])
#     ax2.set_yticklabels([])
#     ax2.set_ylabel(f'{mg_title}', fontsize=fs)
#
#     # Dotted line at x=0
#     ax[0].axvline(0, color='k', linestyle='dashed', linewidth=1)
#     ax[1].axvline(0, color='k', linestyle='dashed', linewidth=1)
#     ax[2].axvline(0, color='k', linestyle='dashed', linewidth=1)
#
#     # Saving plot and displaying
#     save_path = os.path.join(output_loc, f"doppler_profiles_{actual_time}s_{position_solar_y}.png")
#     plt.savefig(save_path, bbox_inches='tight')
#     plt.show()
#     plt.close(fig)


# Plotting Intensity quartiles reference
#def plot_iris_sns_quartile_fits(si_title, cii_title, mg_title, event, main_header, time):


    # Defining common extent
    global_max = np.nanmax([
        np.nanpercentile(si_q_int_map.data, 100 - alpha),
        np.nanpercentile(cii_q_int_map.data, 100 - alpha),
        np.nanpercentile(mg_q_int_map.data, 100 - alpha)])

    global_min = 0
    norm = colors.Normalize(vmin=global_min, vmax=global_max)



# PLOTTING
    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(10, 5))

# Si plotting
    im0 = ax[0].imshow(si_q_int_map.data, origin='lower', cmap='Reds_r', aspect='auto',
        extent=[si_t_array.min(), si_t_array.max(), si_slit_pos.min(),
        si_slit_pos.max()], norm=norm)
    ax[0].set_ylabel(' ')

    # second Si axis for label
    ax_0 = ax[0].twinx()
    ax_0.set_ylabel(si_title, fontsize=fs)
    ax_0.set_yticks([])

    plt.title(f'Height: {position_solar_y}', loc='right', fontsize=fs)


# C ii plotting
    im1 = ax[1].imshow(cii_q_int_map.data, origin='lower', cmap='Reds_r', aspect='auto',
        extent=[cii_t_array.min(), cii_t_array.max(), cii_slit_pos.min(),
        cii_slit_pos.max()], norm=norm)

    ax[1].set_ylabel("Solar Y", fontsize=fs)
    ax[1].set_title(' ')

    #second cii axis for label
    ax_1 = ax[1].twinx()
    ax_1.set_ylabel(cii_title, fontsize=fs)
    ax_1.set_yticks([])


# Mg ii plotting
    im2 = ax[2].imshow(mg_q_int_map.data, origin='lower', cmap='Reds_r', aspect='auto',
        extent=[mg_t_array.min(), mg_t_array.max(), mg_slit_pos.min(), mg_slit_pos.max()],
                       norm=norm)

    ax[2].set_ylabel(" ")
    ax[2].set_xlabel("Time (s)", fontsize=fs)
    ax[2].set_title(' ')

    # second mg axis for label
    ax_2 = ax[2].twinx()
    ax_2.set_ylabel(mg_title, fontsize=fs)
    ax_2.set_yticks([])

    # Defining profile pixel

    target_y = position_solar_y
    max_target_t = max(time)
    min_target_t = min(time)

    for time in time_seconds:
        for axis in ax:
            axis.scatter(time, target_y, marker='x', s=250, c='k', lw=3)
    for axis in ax:
        axis.set_xlim(min_target_t - zoom, max_target_t + zoom)
        axis.set_ylim(target_y - 20, target_y + 20)

# Colorbar
    cbar = fig.colorbar(im2, ax=ax, orientation='vertical', pad=0.12, fraction=0.05)
    cbar.set_label("Integrated Intensity", fontsize=fs)
    tick_locator = ticker.LinearLocator(numticks=4)
    cbar.locator = tick_locator
    cbar.update_ticks()

    save_path = os.path.join(output_loc, f"quartile_maps_{event}_{time_seconds}_{position_solar_y}.png")
    plt.savefig(save_path, bbox_inches="tight")
    plt.show()



#------------- Universal data -------------

# Get headers
hdul = fits.open(fits_file)
hdr = hdul[0].header

# Get line names
si_title = hdr['TDESC' + str(lines[0])]
cii_title = hdr['TDESC' + str(lines[1])]
mg_title = hdr['TDESC' + str(lines[2])]


# Get Si IV properties
si_header = hdul[lines[0]].header
si_wavelength = get_wavelength(si_header,0)


# Get Cii properties
cii_header = hdul[lines[1]].header
cii_wavelength = get_wavelength(cii_header, 1)


# Get mg properties
mg_header = hdul[lines[2]].header
mg_wavelength = get_wavelength(mg_header, 2)


# Doppler velocities
lambda_si = 1402.8
si_v_dopp = ((si_wavelength - lambda_si) / lambda_si) * (c / 1e3)

lambda_cii = 1335.7
cii_v_dopp = ((cii_wavelength - lambda_cii) / lambda_cii) * (c / 1e3)

lambda_mg = 2795.5
mg_v_dopp = ((mg_wavelength - lambda_mg) / lambda_mg) * (c / 1e3)


# Get int maps
si_q_int_map = get_int_map(si_file)
cii_q_int_map = get_int_map(cii_file)
mg_q_int_map = get_int_map(mg_file)

cadence = hdr['STEPT_AV']
alpha = 1

# Si time and position
si_t_array = np.arange(si_q_int_map.data.shape[1]) * cadence
si_slit_pos = si_q_int_map.meta['crval2'] + si_q_int_map.meta['cdelt2'] * (
            np.arange(si_q_int_map.data.shape[0]) - si_q_int_map.meta['crpix2'])

# Cii time and position
cii_t_array = np.arange(cii_q_int_map.data.shape[1]) * cadence
cii_slit_pos = cii_q_int_map.meta['crval2'] + cii_q_int_map.meta['cdelt2'] * (
            np.arange(cii_q_int_map.data.shape[0]) - cii_q_int_map.meta['crpix2'])

# Mg time and postion
mg_t_array = np.arange(mg_q_int_map.data.shape[1]) * cadence
mg_slit_pos = mg_q_int_map.meta['crval2'] + mg_q_int_map.meta['cdelt2'] * (
            np.arange(mg_q_int_map.data.shape[0]) - mg_q_int_map.meta['crpix2'])


def plot_combined_fig():
    fig = plt.figure(figsize=(20, 10))
    fig.text(0.5, 0.96, 'Doppler Velocity (km/s)', ha='center', va='center', fontsize=fs)
    fig.text(0.15, 0.96, f'{position_solar_y} arcsec', ha='center', va='center', fontsize=fs-4)
    gs = gridspec.GridSpec(9, 4, figure=fig, height_ratios=[
        1, 1, 0.3,
        1, 1, 0.3,
        1, 1, 0.3,
    ], hspace=0.14, wspace=0.05)
    plt.subplots_adjust(top=0.92)

    line_info = [
        ('Si IV', 0, si_q_int_map, si_t_array, si_slit_pos, 1402.8),
        ('C II', 1, cii_q_int_map, cii_t_array, cii_slit_pos, 1335.7),
        ('Mg II', 2, mg_q_int_map, mg_t_array, mg_slit_pos, 2795.5)
    ]



    for i, (title, line_idx, q_map, t_arr, slit_pos, rest_wave) in enumerate(line_info):
        data = get_data_arr(line_idx)
        header = hdul[lines[line_idx]].header
        wavelength = get_wavelength(header, line_idx)

        v_dopp, profiles = get_line_profiles(line_idx, rest_wave)


        # Time indices
        cadence = hdr['STEPT_AV']
        time_array = np.arange(data.shape[0]) * cadence
        t_indices = [time_to_index(time_array, t) for t in time_seconds]

        pos_idx = np.argmin(np.abs(slit_pos - position_solar_y))
        max_intensity = data[t_indices, pos_idx].max()

    # TOP ROW
        for j, t_idx in enumerate(t_indices):
            ax = fig.add_subplot(gs[i * 3, j])

            intensity = data[t_idx, pos_idx]

            intensity_norm = intensity / max_intensity

            red = v_dopp >= 0
            blue = v_dopp < 0

            ax.plot(v_dopp, intensity_norm, color='k')
            ax.plot(v_dopp[red], intensity_norm[red], color='red')
            ax.plot(v_dopp[blue], intensity_norm[blue], color='blue')

            ax.axvline(0, linestyle='--', color='k', linewidth=1)
            ax.set_xlim(-x_lim, x_lim)
            ax.set_ylim(0, 1.1)

            if j == 3:
                ax2 = ax.twinx()
                ax2.set_ylabel(f"{title}", fontsize=fs)
                ax2.set_yticks([])


            if j != 0:
                ax.set_yticklabels([' '])

            if j == 0:
                ax.set_ylabel('Intensity', fontsize=fs)

            ax.text(
                0.98, 0.95, f"t = {time_seconds[j]} s",
                transform=ax.transAxes,
                ha='right', va='top',
                fontsize=fs - 3)

            ax.yaxis.set_major_locator(FixedLocator([0.0, 0.50, 1.0]))
            if i == 0:
                ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

            if i != 0:
                ax.set_xticklabels([' '])
                ax.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)



    # BOTTOM ROW QUARTILES
        ax_q = fig.add_subplot(gs[i * 3 + 1, :])

        # Defining common extent
        global_max = np.nanmax([
            np.nanpercentile(si_q_int_map.data, 100 - alpha),
            np.nanpercentile(cii_q_int_map.data, 100 - alpha),
            np.nanpercentile(mg_q_int_map.data, 100 - alpha)])

        global_min = 0
        norm = colors.Normalize(vmin=global_min, vmax=global_max)

        im = ax_q.imshow(
            q_map.data,
            origin='lower',
            aspect='auto',
            cmap='Reds_r',
            extent=[t_arr.min(), t_arr.max(), slit_pos.min(), slit_pos.max()], norm=norm
        )


        # Mark selected points
        for t in time_seconds:

            ax_q.scatter(t, position_solar_y, marker='x', c='k', s=300, lw=4)

        ax_q.set_xlim(min(time_seconds) - zoom, max(time_seconds) + zoom)
        ax_q.set_ylim(position_solar_y - 20, position_solar_y +20)

        ax_q.set_ylabel('Solar Y', fontsize=fs)


        if i == 2:
            ax_q.set_xlabel("Time (s)", fontsize=fs)
    fig.align_ylabels()
    save_path = os.path.join(output_loc, f"line_profiles_{time_seconds}_{position_solar_y}.png")
    plt.savefig(save_path, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    plot_combined_fig()


