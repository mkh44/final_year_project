from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from matplotlib import colors
from matplotlib.gridspec import GridSpec
import os
from datetime import timedelta
from datetime import datetime as dt



def time_to_index(time_array, target_time):
    return np.argmin(np.abs(time_array - target_time))

sji_filepath = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\iris_l2_20230503_072923_4204700135_SJI_2796_t000.fits"

time_seconds = [11500, 11550, 11575, 11600, 11625, 11650, 11700, 11750]
position_solar_y = 270


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
vmax = np.nanpercentile(np.stack(frames), 99)

n = len(time_seconds)



fig = plt.figure(figsize=(18, 10))
gs = GridSpec(2, 4, figure=fig, wspace=0.07, hspace=0.07)



for i, (t, idx) in enumerate(zip(time_seconds, sji_indices)):

    wcs_2d = wcs.slice([idx, slice(None), slice(None)])


    ax = fig.add_subplot(gs[i // 4, i%4], projection=wcs_2d)

    img = data[idx, :, :]

    im = ax.imshow(img, origin='lower', cmap='magma', vmin=vmin, vmax=vmax)
    ax.axhline(y_pix, color='white', linestyle='--', linewidth=2)
    ax.axvline(crpix2 - 9, color='white', linewidth=3)



    if i == 0 or i== 4:
        ax.coords[1].set_ticklabel_visible(True)
    else:
        ax.coords[1].set_ticklabel_visible(False)

    ax.set_ylabel(' ')
    ax.set_xlabel(" ")
    ax.text(0.1, 2, f"{time_labels[i]}", color='white', fontsize= 12)
    ax.text(0.1, 15, f"{position_solar_y}\"", color='white', fontsize=12)

fig.text(0.48, 0.06,"Solar x (arcsec)")
fig.text(0.09, 0.48, 'Solar y (arcsec)', rotation=90)


save_paths = os.path.join(output, f"SJI_{time_seconds}_{position_solar_y}.png")
plt.savefig(save_paths, bbox_inches="tight")

plt.show()

