from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from matplotlib import colors
from matplotlib.gridspec import GridSpec
import os
from datetime import timedelta
from datetime import datetime as dt
import matplotlib.patheffects as pe
import astropy.units as u

from scipy.ndimage import rotate

fs = 20

def time_to_index(time_array, target_time):
    return np.argmin(np.abs(time_array - target_time))

sji_filepath = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\iris_l2_20230503_072923_4204700135_SJI_2796_t000.fits"

time_seconds = [7050, 7090, 7110]
position_solar_y = 268


output = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\sji"


hdul = fits.open(sji_filepath)
hdul.info()
data = hdul[0].data
hdr = hdul[0].header

#times real
cadence = hdr.get('CDELT3')
start_time_str = hdr['STARTOBS']
start_time = dt.fromisoformat(start_time_str)
clock_times = [start_time + timedelta(seconds=t) for t in time_seconds]
time_labels = [ct.strftime('%H:%M:%S') for ct in clock_times]


# headers
cdelt2 = hdr['CDELT2']
crval2 = hdr['CRVAL2']
crpix2 = hdr['CRPIX2']


y_pix = (position_solar_y - crval2) / cdelt2 + (crpix2 - 1)

time_sji = np.arange(data.shape[0]) *cadence

sji_indices = [time_to_index(time_sji, t) for t in time_seconds]


wcs = WCS(hdul[0].header)
frame = 100
frames = [data[idx, :, :] for idx in sji_indices]
frame_idx = sji_indices[0]

vmin = 0
vmax = np.nanpercentile(np.stack(frames), 98.8)



if len(time_seconds) <= 4:
    rows = 1
else:
    rows = 2

fig = plt.figure(figsize=(18, 10))
gs = GridSpec(rows, 4, figure=fig, wspace=0.07, hspace=0.07)



for i, (t, idx) in enumerate(zip(time_seconds, sji_indices)):

    wcs_2d = wcs.slice([idx, slice(None), slice(None)])


    ax = fig.add_subplot(gs[i // 4, i%4], projection=wcs_2d)

    img = data[idx, :, :]

    im = ax.imshow(img, origin='lower', cmap='magma_r', vmin=vmin, vmax=vmax)

    # lines to show pixel
    ax.axhline(y_pix, color='white', linestyle='--', linewidth=2)
    ax.axvline(crpix2 - 9, color='white', linewidth=3)


    if rows == 2:
        fig.text(0.45, 0.95, "Solar x (arcsec)", fontsize=fs)
        fig.text(0.06, 0.45, 'Solar y (arcsec)', rotation=90, fontsize=fs)

        if i == 0 or i== 4:
            ax.coords[1].set_ticklabel_visible(True)
            ax.coords[1].tick_params(axis='y', labelsize=fs)

            ax.coords[1].set_ticks(spacing=20 * u.arcsec)
            ax.coords[1].set_ticks_position('l')
            ax.coords[1].set_ticklabel_position('l')
        else:
            ax.coords[1].set_ticklabel_visible(False)
            ax.coords[1].set_ticks(spacing=20 * u.arcsec)
            ax.coords[1].set_ticks_position('l')
            ax.coords[1].set_ticklabel_position('l')
        if i < 4:
            ax.coords[0].set_ticklabel_visible(True)
            ax.coords[0].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, axis='x', labelsize=fs)
            ax.coords[0].set_axislabel(' ')
        else:
            ax.coords[0].set_ticklabel_visible(False)

    else:
        fig.text(0.45, 0.75, "Solar x (arcsec)", fontsize=fs)
        fig.text(0.06, 0.4,'Solar y (arcsec)', rotation= 90, fontsize=fs)
        if i == 0:
            ax.coords[1].set_ticklabel_visible(True)
            ax.coords[1].tick_params(axis='y', labelsize=fs)
            ax.coords[1].set_ticks(spacing=20 * u.arcsec)
            ax.coords[1].set_ticks_position('l')
            ax.coords[1].set_ticklabel_position('l')
        else:
            ax.coords[1].set_ticklabel_visible(False)
            ax.coords[1].set_ticks(spacing=20 * u.arcsec)
            ax.coords[1].set_ticks_position('l')
            ax.coords[1].set_ticklabel_position('l')
        if i <= 4:
            ax.coords[0].tick_params(top=True, labeltop=True, bottom=True, labelbottom=False, axis='x', labelsize=fs)
            ax.coords[0].set_axislabel(' ')
        else:
            ax.coords[0].set_ticklabel_visible(False)


    ax.set_ylabel(' ')
    ax.set_xlabel(" ")
    ax.text(3, 2, f"{time_labels[i]}", color='k', fontsize= fs)
    ax.text(3.5, 16, f"{position_solar_y}\"", color='k', fontsize=fs)



save_paths = os.path.join(output, f"SJI_{time_seconds}_{position_solar_y}.png")
plt.savefig(save_paths, bbox_inches="tight")

plt.show()

