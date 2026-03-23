from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from matplotlib import colors
from matplotlib.gridspec import GridSpec
import os



def time_to_index(time_array, target_time):
    return np.argmin(np.abs(time_array - target_time))

sji_filepath = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\iris_l2_20230503_072923_4204700135_SJI_2796_t000.fits"

time_seconds = [11800, 11900, 12000, 12100]
position_solar_y = 253
output = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\sji"

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
frames = [data[idx, :, :] for idx in sji_indices]
frame_idx = sji_indices[0]

vmin = 0
vmax = np.nanpercentile(frames, 99)


n = len(time_seconds)

fig = plt.figure(figsize=(16, 4))
gs = GridSpec(1, n, figure=fig, wspace=0.01)

for i, (t, idx) in enumerate(zip(time_seconds, sji_indices)):

    wcs_2d = wcs.slice([idx, slice(None), slice(None)])


    ax = fig.add_subplot(gs[0, i], projection=wcs_2d)

    img = data[idx, :, :]

    im = ax.imshow(img, origin='lower', cmap='magma', vmin=vmin, vmax=vmax)
    ax.axhline(y_pix, color='black', linestyle='--', linewidth=2)
    ax.axvline(crpix2 - 9, color='black', linewidth=3)

    if i == 0:
        ax.set_ylabel("Solar Y (arcsec)")
    else:
        ax.set_ylabel(' ')
        ax.coords[1].set_ticklabel_visible(False)



    ax.set_title(f"{t}s, {position_solar_y}\"")
    ax.set_xlabel("Solar X (arcsec)")

save_paths = os.path.join(output, f"SJI_{time_seconds}_{position_solar_y}.png")
plt.savefig(save_paths, bbox_inches="tight")

plt.show()


