#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%
import glob
import os
import pdb
import sunpy.map
#from iris_fitting.fit_iris_lines import fit_raster
#from iris_fitting import iris_get_mg_features_lv2 as get_mg
#from iris_fitting import get_mgii_quartiles
import asdf
from astropy.io import fits
import datetime as dt
from datetime import timedelta
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

from sunpy.net import fido_factory, attrs as a
#from iris_fitting import extract_irisL2data
from astropy import units as u
import tarfile
from scipy.constants import speed_of_light


# Plot the output of the fitting routine if fitting sit-and-stare data
def plot_iris_sns_fits(int_map,dopp_map,width_map,vnt_map,asym_map,iris_window,event,main_header,asym_rng,dopp_rng,max_wid,max_vnt):

    fig = plt.figure(constrained_layout=True, figsize=(10, 10))
    plt.rcParams['font.size'] = '10'

    if vnt_map is not None:
        nrows = 5
    else:
        nrows = 4
    gs = GridSpec(nrows=nrows, ncols=1, hspace=0.05, wspace=0.05)
    gs.update(left=0.08, right=0.95, bottom=0.06, top=0.95)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

    # Absolute time of file
    obs_start = dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f')
    abs_start = obs_start + dt.timedelta(seconds=float(t_start))
    abs_end = obs_start + dt.timedelta(seconds=float(t_end))
    abs_start_str = abs_start.strftime('%H:%M:%S')
    abs_end_str = abs_end.strftime('%H:%M:%S')



# Set the plotting parameters
    cadence = main_header['STEPT_AV']
    t_array = np.arange(0, int_map.data.shape[1])*cadence
    t_mask = (t_array >= t_start) & (t_array <= t_end)
    t_inds = np.where(t_mask)[0]
   
# Slice raster maps in time
    int_data   = int_map.data[:, t_start:t_end]
    dopp_data  = dopp_map.data[:, t_start:t_end]
    width_data = width_map.data[:, t_start:t_end]
    asym_data  = asym_map.data[:, t_start:t_end]
    #Only create vnt_data if vnt_map exists
    if vnt_map is not None:
        vnt_data = vnt_map.data[:, t_start:t_end]
    else:
        vnt_data = None


    t_plot = t_array[t_start:t_end]

    # Slit position
    slit_pos = int_map.meta['crval2'] + int_map.meta['cdelt2'] * (np.arange(int_map.data.shape[0]) - int_map.meta['crpix2'])

    os.makedirs(os.path.join(output_loc), exist_ok=True)
# Intensity map
    ax1 = fig.add_subplot(gs[0,0], label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(int_data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    
    plt.imshow(int_data, norm=norm, cmap = mpl.colormaps['Reds_r'], axes=ax1, extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax1.set_ylabel(" ")
    ax1.set_xlabel(" ")
    ax1.set_xticklabels([])

    plt.colorbar(location='right', label=r'a) Intensity', shrink=0.6, ax = ax1, ticks=[0,upr_bnd])

# Asymmetry map
    ax2 = fig.add_subplot(gs[1,0], label='b)')
    norm = colors.Normalize(vmin = -asym_rng, vmax = asym_rng)
    plt.imshow(asym_data, norm=norm, cmap = mpl.colormaps['seismic'], axes=ax2, extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax2.set_ylabel(" ")
    ax2.set_xlabel(" ")
    ax2.set_xticklabels([])

    plt.colorbar(location='right', label=r'b) RB Asym.', shrink=0.6, ax = ax2)

# Doppler map
    ax3 = fig.add_subplot(gs[2,0], label='c)')
    norm = colors.Normalize(vmin = -dopp_rng, vmax = dopp_rng)
    plt.imshow(dopp_data, norm=norm, cmap = mpl.colormaps['coolwarm'], axes=ax3, extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax3.set_ylabel(" ")
    ax3.set_xlabel(" ")
    ax3.set_xticklabels([])

    plt.colorbar(location='right', label=r'c) v$_{dopp}$ ($km~s^{-1}$)', shrink=0.6, ax = ax3)

# Line width
    ax4 = fig.add_subplot(gs[3,0], label='d)')
    norm = colors.Normalize(vmin = 0, vmax = max_wid)
    plt.imshow(width_data, norm=norm, cmap = mpl.colormaps['cubehelix'], axes=ax4, extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax4.set_ylabel(" ")
    ax4.set_xlabel(" ")
    if vnt_data is not None:
        ax4.set_xticklabels([])

    plt.colorbar(location='right', label=r'd) Width ($\AA$)', shrink=0.6, ax = ax4)
    
# Nonthermal velocity
    if vnt_data is not None:
        ax5 = fig.add_subplot(gs[4,0], label='e)')
        norm = colors.Normalize(vmin = 0, vmax = max_vnt)
        plt.imshow(vnt_data, norm=norm, cmap = mpl.colormaps['inferno'], axes=ax5, extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
        plt.colorbar(location='right', label=r'e) v$_{nt}$ ($km~s^{-1}$)', shrink=0.6, ax = ax5)
        ax5.set_ylabel(" ")
        ax5.set_xlabel(" ")

    fig.supxlabel("Time from raster start (s)")
    fig.supylabel('Solar Y (arcsec)')
    plt.suptitle(iris_window + r'$\AA$; ' + abs_start_str + '-' + abs_end_str)
    plt.savefig(os.path.join(output_loc, event) + f"/IRIS_zoomed_plot_"+iris_window.replace(' ', '_')+'_'+file_time+'_'+f"{t_start:.0f}-{t_end:.0f}s_from_raster_start.png", bbox_inches='tight')
    plt.close(fig)


def plot_y_for_time_series(int_map, dopp_map, width_map, vnt_map, asym_map,
                      iris_window, event, main_header,
                      y_value_arcsec=None, y_index=None):
    #Cadence
    cadence = main_header.get('STEPT_AV')
    nt = int_map.data.shape[1]
    t_array = np.arange(nt) * cadence

    # Slit pos
    slit_pos = int_map.meta['crval2'] + int_map.meta['cdelt2'] * (
    np.arange(int_map.data.shape[0]) - int_map.meta['crpix2'])

    # Choose Y
    y_index = np.argmin(np.abs(slit_pos - y_value_arcsec))

    y_selected = slit_pos[y_index]

    #Get time series
    int_ts = int_map.data[y_index, :]
    dopp_ts = dopp_map.data[y_index, :]
    width_ts = width_map.data[y_index, :]
    asym_ts = asym_map.data[y_index, :]
    vnt_ts = vnt_map.data[y_index, :] if vnt_map is not None else None

    # Check frequency
    fft = np.fft.rfft(int_ts - np.nanmean(int_ts))
    freq = np.fft.rfftfreq(len(int_ts), cadence)

    valid = freq > 0

    plt.plot(1/freq[valid], np.abs(fft[valid]))
    plt.xlabel("Period (s)")
    plt.ylabel("Amplitude")
    plt.show()


    #Plotting
    fig, ax = plt.subplots(4 if vnt_map is None else 5, 1, figsize=(10, 8), sharex=True)
    ax[0].plot(t_array, int_ts)
    ax[0].set_ylabel("Intensity")

    ax[1].plot(t_array, asym_ts)
    ax[1].set_ylabel("RB Asym")

    ax[2].plot(t_array, dopp_ts)
    ax[2].set_ylabel("v_dopp (km/s)")

    ax[3].plot(t_array, width_ts)
    ax[3].set_ylabel("Width (Å)")

    if vnt_ts is not None:
        ax[4].plot(t_array, vnt_ts)
        ax[4].set_ylabel("v_nt (km/s)")

    ax[-1].set_xlabel("Time (s)")

    plt.suptitle(f"{iris_window} - Y = {y_selected:.2f} arcsec")
    plt.tight_layout()

    os.makedirs(os.path.join(output_loc, y_loc), exist_ok=True)

    plt.savefig(os.path.join(y_loc, "y_plots") + iris_window.replace(' ', '_') + f"_time_series_Y_{y_selected:.1f}.png")
    plt.close(fig)

# EVENT
event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive - Dublin City University\PHA4\Final_Year_Project\outputs"
y_loc = os.path.join(output_loc, "y_plots")
os.makedirs(y_loc, exist_ok=True)

# Define time frame
t_start = 10000    # seconds from raster start
t_end   = 13000    # seconds from raster start

asdf_files = glob.glob(os.path.join(output_loc, "*.asdf"))

if not asdf_files:
    raise FileNotFoundError(f"No ASDF files found in {output_loc}")



# Get map names
#----------------------
for asdf_file in asdf_files:
    iris_window_underscore = os.path.basename(asdf_file).replace('IRIS_fitting_', '').replace(event, '').replace('.asdf', '')
    iris_window = iris_window_underscore.replace('_', ' ').strip()


    print(iris_window)
    with asdf.open(asdf_file) as af:
        possible_int_keys = [key for key in af.tree.keys() if 'int' in key.lower()]
        if not possible_int_keys:
            raise KeyError(f"No intensity map found in {asdf_file}")
        int_key = possible_int_keys[0]
        int_map = af.tree[int_key]

        possible_dopp_keys = [key for key in af.tree.keys() if 'dopp' in key.lower()]
        if not possible_dopp_keys:
            raise KeyError(f"No Doppler map found in {asdf_file}")
        dopp_key = possible_dopp_keys[0]
        dopp_map = af.tree[dopp_key]

        possible_width_keys = [key for key in af.tree.keys() if 'width' in key.lower()]
        if not possible_width_keys:
            raise KeyError(f"No width map found in {asdf_file}")
        width_key = possible_width_keys[0]
        width_map = af.tree[width_key]

        possible_vnt_keys = [key for key in af.tree.keys() if 'vnt' in key.lower()]
        if possible_vnt_keys:
            vnt_key = possible_vnt_keys[0]
            vnt_map = af.tree[vnt_key]
        else:
            vnt_map = None
            print(f'Warning: No VNT map found in {asdf_file}. Skipping vnt plot.')


        possible_asym_keys = [key for key in af.tree.keys() if 'asym' in key.lower()]
        if not possible_asym_keys:
            raise KeyError(f"No asym map found in {asdf_file}")
        asym_key = possible_asym_keys[0]
        asym_map = af.tree[asym_key]

    # Find ranges

    # Define per-window plotting ranges
    plot_ranges = {
        "Si IV 1394": {"dopp_rng": 10, "max_wid": 0.1, "asym_rng": 1, "max_vnt": 30},
        "C II 1334": {"dopp_rng": 10, "max_wid": 0.1, "asym_rng": 1, "max_vnt": 30},
        "C II 1335": {"dopp_rng": 10, "max_wid": 0.1, "asym_rng": 1, "max_vnt": 30},
        "smooth MgII": {"dopp_rng": 10, "max_wid": 0.8, "asym_rng": 1, "max_vnt": 30},}

    # Default if a window is not in the dictionary
    default_ranges = {"dopp_rng": 10, "max_wid": 0.1, "asym_rng": 1, "max_vnt": 30}

    # Get cadence from .fits file as asdf does not contain it
    iris_fits = glob.glob(os.path.join(r"C:\Users\molly\Downloads\iris", "*.fits"))
    main_header = fits.getheader(iris_fits[0], 0)

    # Ensure output directory exists
    ranges = plot_ranges.get(iris_window, default_ranges)
    dopp_rng = ranges["dopp_rng"]
    max_wid = ranges["max_wid"]
    asym_rng = ranges["asym_rng"]
    max_vnt = ranges["max_vnt"]

    print(dopp_rng, max_wid, max_vnt, asym_rng)

    #plot_iris_sns_fits(int_map, dopp_map, width_map, vnt_map, asym_map, iris_window, event, main_header, asym_rng, dopp_rng, max_wid, max_vnt)
    plot_y_for_time_series(int_map, dopp_map, width_map, vnt_map, asym_map, iris_window, event, main_header, y_value_arcsec=235)



