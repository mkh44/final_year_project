
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import glob
import os
from astropy.wcs import WCS


# Set up some default matplotlib options
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['image.origin'] = 'lower'
plt.rcParams['image.cmap'] = 'viridis'

input_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"
file = glob.glob(os.path.join(input_loc, "iris_l2_20230503_072923_4204700135_raster_t000_r00000.fits"))[0]

f = fits.open(file)


data = f[1].data
hd = f[1].header

for i in range(1, 10):
    print(i, f[i].header.get('TWAVE1'), f[i].header.get('TDESC'))
plt.imshow(f[1].data[0], vmin=0, vmax=10)

plt.show()


