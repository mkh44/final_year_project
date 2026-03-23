from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

sji_filepath = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs\iris_l2_20230503_072923_4204700135_SJI_2796_t000.fits"

hdul = fits.open(sji_filepath)
hdul.info()
data = hdul[0].data

plt.imshow(data[0], origin='lower', cmap='gray')
plt.colorbar(label="Intensity")
plt.title("SJI Frame 0")
plt.show()
