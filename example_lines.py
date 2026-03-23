from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from matplotlib.gridspec import GridSpec


def time_to_index(time_array, target_time):
    return np.argmin(np.abs(time_array - target_time))

sji_filepath = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\iris_l2_20230503_072923_4204700135_SJI_2796_t000.fits"

time_seconds = [9700, 9750, 9800, 9850]
position_solar_y = 267


hdul = fits.open(sji_filepath)
hdul.info()
data = hdul[0].data
hdr = hdul[0].header
cadence = hdr.get('CDELT3')

cdelt2 = hdr['CDELT2']
crval2 = hdr['CRVAL2']
crpix2 = hdr['CRPIX2']


y_pix = (position_solar_y - crval2) / cdelt2 + (crpix2 - 1)

time_sji = np.arange(data.shape[0]) *cadence

sji_indices = [time_to_index(time_sji, t) for t in time_seconds]


wcs = WCS(hdul[0].header)
frame = 100

frame_idx = sji_indices[0]  # or loop through indices


wcs_2d = wcs.slice([frame_idx, slice(None), slice(None)])


# fig, axes = plt.subplots(1, len(time_seconds), figsize=(16, 4), subplot_kw={'projection': wcs_2d}, constrained_layout=True)
#
# for i, (t, idx) in enumerate(zip(time_seconds, sji_indices)):
#     img = data[idx, :, :]
#
#     vmin = np.percentile(img, 0.5)
#     vmax = np.percentile(img, 99.5)
#
#     axes[i].imshow(img, origin='lower', cmap='magma', vmin=vmin, vmax=vmax)
#
#
#     axes[i].axhline(y_pix, color='k', linestyle='--', linewidth=2)
#     axes[i].axvline(crpix2 - 9, color='k', linewidth=2)
#     axes[i].set_xlabel("Solar X (arcsec)")
#
#     axes[i].set_title(f"t = {t}s")
#
#     if i == 0:
#         axes[i].set_ylabel("Solar Y (arcsec)")
#     else:
#         axes[i].set_ylabel(' ')
n = len(time_seconds)


fig = plt.figure(figsize=(16, 4))
gs = GridSpec(1, n, figure=fig, wspace=0.01)

for i, (t, idx) in enumerate(zip(time_seconds, sji_indices)):

    wcs_2d = wcs.slice([idx, slice(None), slice(None)])


    ax = fig.add_subplot(gs[0, i], projection=wcs_2d)

    img = data[idx, :, :]
    vmin = np.percentile(img, 0.5)
    vmax = np.percentile(img, 99.5)

    im = ax.imshow(img, origin='lower', cmap='magma', vmin=vmin, vmax=vmax)
    ax.axhline(y_pix, color='black', linestyle='--', linewidth=2)


    if i == 0:
        ax[i].set_ylabel("Solar Y (arcsec)")
    else:
        axes[i].set_ylabel(' ')

    ax.set_title(f"t = {t}s")
    ax.set_xlabel("Solar X (arcsec)")
    ax.set_ylabel("Solar Y (arcsec)")



plt.show()


