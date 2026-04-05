
import pdb
import glob
import os

from astropy.wcs.docstrings import crval
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
from datetime import datetime, timedelta
from get_quartiles import get_quartiles
import matplotlib.dates as mdates
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
time_seconds = [12000]
position_solar_y = 253

# [6900, 6950, 7000, 7050]
# [9600, 9700, 9750, 9800]
# [11500, 11750, 11900, 12050]


z = 10
x_lim = 450

zoom = 50
fs = 22 # font size

def time_to_index(time_array, target_time):
    return np.argmin(np.abs(time_array - target_time))

def seconds_to_realtime(time_seconds, start_time):
    return start_time + timedelta(seconds=time_seconds)

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

def get_velocity_map(file):
    with asdf.open(file) as af:
        if 'q_dopp_map' in af.tree:
            return af.tree['q_dopp_map']

        else:
            raise KeyError("No velocity map found")

def get_line_profiles(line, expected_wl):
# Getting pixel array
    data_arr = get_data_arr(line)
    title = hdul[lines[line]].header

# time index
    t_indices = [time_to_index(time_array, t) for t in time_seconds]
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

n = len(time_seconds)

def plot_combined_fig():
    fig = plt.figure(figsize=(20,14))
    fig.text(0.5, 0.96, 'Doppler Velocity (km/s)', ha='center', va='center', fontsize=fs)
    fig.text(0.15, 0.96, f'{position_solar_y} arcsec', ha='center', va='center', fontsize=fs)
    gs = gridspec.GridSpec(9, n, figure=fig, height_ratios=[
        1, 1, 0.08,
        1, 1, 0.08,
        1, 1, 0.08,
    ], hspace=0.22, wspace=0.05)
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


            if j == 0:
                ax.set_ylabel('Intensity', fontsize=fs-3)
                ax.tick_params(labelsize=fs-5)
            else:
                ax.set_yticklabels([' '])

            if i ==2:
                ax.text(
                0.05, 0.95, f"{time_labels[t_idx]}",
                transform=ax.transAxes,
                ha='left', va='top',
                fontsize=fs - 3)
            else:
                ax.text(
                    0.98, 0.95, f"{time_labels[t_idx]}",
                    transform=ax.transAxes,
                    ha='right', va='top',
                    fontsize=fs - 3)

            ax.yaxis.set_major_locator(FixedLocator([0.0, 0.50, 1.0]))

            if i == 0:
                ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
                ax.tick_params(labelsize=fs - 5)

            else:
                ax.set_xticklabels([' '])
                ax.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)

        # Twin axis for line names on right
            if j == n-1:
                ax2 = ax.twinx()
                ax2.set_ylabel(f"{title}", fontsize=fs-4)
                ax2.set_yticks([])

    # BOTTOM ROW QUARTILES
        ax_q = fig.add_subplot(gs[i * 3 + 1, :])

        # Defining common extent
        global_max = np.nanmax([
            np.nanpercentile(si_q_int_map.data, 99.9),
            np.nanpercentile(cii_q_int_map.data, 99,),
            np.nanpercentile(mg_q_int_map.data, 99.9)])

        global_min = 0
        norm = colors.Normalize(vmin=global_min, vmax=global_max)

        t_arr_clock = np.array([start_time + timedelta(seconds=t) for t in t_arr])


        im = ax_q.imshow(
            q_map.data,
            origin='lower',
            aspect='auto',
            cmap='Reds_r',
            extent=[t_arr_clock.min(), t_arr_clock.max(), slit_pos.min(), slit_pos.max()], norm=norm
        )


        # Mark selected points
        selected_times = [clock_times[idx] for idx in t_indices]
        for t in selected_times:
            ax_q.scatter(t, position_solar_y, marker='x', c='k', s=300, lw=4)

        # Set axis limits
        min_x = seconds_to_realtime((min(time_seconds) - zoom), start_time)
        max_x = seconds_to_realtime((max(time_seconds) + zoom), start_time)
        ax_q.set_xlim(min_x, max_x)
        ax_q.set_ylim(position_solar_y - 20, 280)

        ax_q.set_ylabel('Solar y', fontsize=fs-3)
        ax_q.tick_params(labelsize=fs-5)
        # ax_q.xaxis_date()
        ax_q.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        # fig.autofmt_xdate()

        if i == 2:
            ax_q.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
            ax_q.tick_params(labelsize=fs-5)
        else:
            ax_q.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)

    fig.text(0.5, 0.09, "Time", fontsize=fs)
    fig.align_ylabels()
    save_path = os.path.join(output_loc, f"line_profiles_{time_seconds}_{position_solar_y}.png")
    plt.savefig(save_path, bbox_inches="tight")
    plt.show()



#------------- Universal data -------------

# Get headers
hdul = fits.open(fits_file)
hdr = hdul[0].header

# Get start time
start_time = datetime.fromisoformat(hdr['STARTOBS'])
time_array = np.arange(hdul[1].data.shape[0]) * hdr['STEPT_AV']

# Clock times
clock_times = np.array([start_time + timedelta(seconds=t) for t in time_array])
time_labels = [ct.strftime('%H:%M:%S') for ct in clock_times]

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


def plot_wavelength_with_velocity(time_seconds):
    fig = plt.figure(figsize=(16, 12))
    #fig.text(0.5, 0.96, 'Doppler Velocity (km/s)', ha='center', va='center', fontsize=fs, color='#8c0010')
    fig.text(0.15, 0.06, f'{position_solar_y} arcsec', ha='center', va='center', fontsize=fs)
    # fig.text(0.48, 0.05, 'Wavelength (Å)', fontsize=fs, color='blue')
    gs = gridspec.GridSpec(3, n, figure=fig, hspace=0.22, wspace=0.05)
    plt.subplots_adjust(top=0.92)


    line_info = [
        ('Si IV', 0, si_q_int_map, si_t_array, si_slit_pos, 1402.8),
        ('C II', 1, cii_q_int_map, cii_t_array, cii_slit_pos, 1335.7),
        ('Mg II', 2, mg_q_int_map, mg_t_array, mg_slit_pos, 2795.5)
    ]

    # Handles for legend
    dopp_handle=None
    wl_handle=None
    restwl_handle=None

    # For each time in timeseconds
    for i, (title, line_idx, q_map, t_arr, slit_pos, rest_wave) in enumerate(line_info):
        data = get_data_arr(line_idx)
        header = hdul[lines[line_idx]].header
        wavelength = get_wavelength(header, line_idx)

        v_dopp, profiles = get_line_profiles(line_idx, rest_wave)

        # Time indices
        t_indices = [time_to_index(time_array, t) for t in time_seconds]

        pos_idx = np.argmin(np.abs(slit_pos - position_solar_y))
        max_intensity = data[t_indices, pos_idx].max()


        for j, t_idx in enumerate(t_indices):
            ax = fig.add_subplot(gs[i, j])

            intensity = data[t_idx, pos_idx]
            intensity_norm = intensity / max_intensity

            red = v_dopp >= 0
            blue = v_dopp < 0

            # Doppler velocity plots
            dopp_line, = ax.plot(
                v_dopp, intensity_norm,
                color='red',
                linestyle='dashdot',
                linewidth=3,
                label='Doppler Velocity (km/s)')
            # ax.plot(v_dopp[red], intensity_norm[red], color='red', lw=3)
            # ax.plot(v_dopp[blue], intensity_norm[blue], color='blue', lw=3)

            # Central line for Dopp vel
            ax.axvline(0, linestyle='--', color='k', linewidth=1.5)
            ax.set_xlim(-x_lim, x_lim)
            ax.set_ylim(0, 1.1)

            # Twin axis for wavelengths
            axw = ax.twiny()
            wl_line, = axw.plot(
                wavelength, intensity_norm,
                color='blue',
                label='Wavelength (Å)')
            axw.tick_params(axis='x', labelcolor='blue', labelsize = fs-3)
            axw.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

            # Central line for rest wavelength (NTS: should match the dopp vel line exactly)
            rest_wl_line = axw.axvline(
                rest_wave,
                linestyle='--',
                color='k',
                linewidth=1.5,
                label='Rest Wavelength/Zero Doppler shift')

            # Axis limits for wavelength
            if title == 'Mg II':
                axw.set_xlim((rest_wave - 4), rest_wave + 4)
            else:
                axw.set_xlim((rest_wave - 2), rest_wave + 2)

            # Y axis title
            if j == 0 and i == 1:
                ax.set_ylabel('Intensity', fontsize=fs)
                ax.tick_params(labelsize=fs - 5)
            elif j != 0:
                ax.set_yticklabels([' '])
            elif j == 0:
                ax.tick_params(labelsize=fs - 5)

            # Y axis tick spacing
            ax.yaxis.set_major_locator(FixedLocator([0.0, 0.50, 1.0]))

            # Time stamps
            ax.text(
                0.05, 0.95, f"{time_labels[t_idx]}",
                transform=ax.transAxes,
                ha='left', va='top',
                fontsize=fs - 3)

            #Wl axis label
            if i == 2 and j == 0:
                axw.xaxis.set_label_position('bottom')
                axw.set_xlabel('Wavelength (Å)', fontsize=fs, color='blue')

            # Dopp axis label
            if i == 0 and j == 0:
                ax.xaxis.set_label_position('top')
                ax.set_xlabel('Doppler Velocity (km/s)', fontsize=fs, color='#8c0010')

            # Dopp tick labels
            if i == 0:
                ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
                ax.tick_params(axis='x', labelsize=fs - 5, labelcolor='#8c0010')

            else:
                ax.set_xticklabels([' '])
                ax.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)

            # Twin axis for line names on right
            if j == n - 1:
                ax2 = ax.twinx()
                ax2.set_ylabel(f"{title}", fontsize=fs)
                ax2.set_yticks([])

            # Legend
            if i == 0 and j == n -1:
                ax.legend(
                    handles=[dopp_line, wl_line, rest_wl_line],
                    loc='upper right',
                    fontsize=fs - 5,
                    frameon=False
                )
            # Final legend handles
            # if dopp_handle is None:
            #     dopp_handle = dopp_line
            # if wl_handle is None:
            #     wl_handle = wl_line
            # if restwl_handle is None:
            #     restwl_handle = rest_wl_line

    fig.align_ylabels()

    # Legend

    # fig.legend(
    #     handles=[dopp_handle, wl_handle, restwl_handle],
    #     loc='upper center',
    #     ncols=1,
    #     fontsize=fs - 3
    # )
    plt.tight_layout()
    plt.show()


# if __name__ == "__main__":
#     plot_combined_fig()

plot_wavelength_with_velocity(time_seconds)
