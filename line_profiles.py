
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
time_seconds = [6950, 7000, 7050, 7075, 7090, 7100, 7110, 7120]
position_solar_y = 268

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

def plot_wl(lines, time):
    # Plotting wavelength
    clock_time = start_time + timedelta(seconds=time)
    time_label = clock_time.strftime('%H:%M:%S')

    si_data_arr = get_data_arr(0)
    cii_data_arr = get_data_arr(1)
    mg_data_arr = get_data_arr(2)

    t_idx = time_to_index(time_array, time)

    pos_idx = np.argmin(np.abs(si_slit_pos - position_solar_y))
    max_intensity = max(si_data_arr[t_idx, pos_idx].max(), cii_data_arr[t_idx, pos_idx].max(), mg_data_arr[t_idx, pos_idx].max())

    si_intensity = si_data_arr[t_idx, pos_idx]
    si_intensity_norm = si_intensity / max_intensity

    cii_intensity = cii_data_arr[t_idx, pos_idx]
    cii_intensity_norm = cii_intensity / max_intensity

    mg_intensity = mg_data_arr[t_idx, pos_idx]
    mg_intensity_norm = mg_intensity / max_intensity


    fig, ax = plt.subplots(3, 1)
    ax[0].plot(si_wavelength, si_intensity_norm, color='#8c0010')
    ax[0].set_xlabel(' ')
    ax[0].set_ylim(0, 1.1)
    ax0 = ax[0].twinx()
    ax0.set_yticks([])
    ax0.set_yticklabels([])
    ax0.set_ylabel(f'{si_title}')
    ax[0].xaxis.set_minor_locator(MultipleLocator(10))
    ax[0].axvline(lambda_si, color='k', linestyle='--')
    ax[0].set_xlim(lambda_si - 3, lambda_si + 3)
    plt.title(f'{time_label}, {position_solar_y}\"', loc='right')

    ax[1].plot(cii_wavelength, cii_intensity_norm, color='#8c0010')
    ax[1].set_xlabel(' ')
    ax[1].set_ylabel("Intensity")
    ax[1].set_ylim(0, 1.1)
    ax1 = ax[1].twinx()
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_ylabel(f'{cii_title}')
    ax[1].axvline(lambda_cii, color='k', linestyle='--')
    ax[1].set_xlim(lambda_cii - 3, lambda_cii + 3)
    ax[1].xaxis.set_minor_locator(MultipleLocator(10))

    ax[2].plot(mg_wavelength, mg_intensity_norm, color='#8c0010')
    ax[2].set_ylabel(' ')
    ax[2].set_xlabel('Wavelength (Å)')
    ax[2].set_ylim(0, 1.1)
    ax2 = ax[2].twinx()
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax2.set_ylabel(f'{mg_title}')
    ax[2].xaxis.set_minor_locator(MultipleLocator(10))
    ax[2].set_xlim((lambda_mg - z), lambda_mg + z)
    ax[2].axvline(lambda_mg, color='k', linestyle='--')
    plt.show()

if __name__ == "__main__":
    plot_combined_fig()

#plot_wl(lines, [11500, 11550, 11575, 11600, 11625, 11650, 11700, 11750])
