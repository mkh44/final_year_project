#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#

import matplotlib.pyplot as plt
import numpy as np
from fit_iris_lines import fit_raster
from zoomed_plots import output_loc


def plot_pixel_spectrum(file, iris_window, xpix, ypix):
    # Load raster
    a = fit_raster(file, iris_window, fulldisk=False)
    data, main_header, header, wavelength = a.open_iris_file(iris_window)

    # Extract spectrum at pixel
    spectrum = data[ypix, xpix, :]

    # Plot
    plt.figure(figsize=(6, 4))
    plt.plot(wavelength, spectrum, color='black')
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Intensity")
    plt.title(f"{iris_window} | Pixel (x={xpix}, y={ypix})")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_loc, ) + iris_window.replace(' ', '_') + .png")
    plt.close(fig)

event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"
plot_pixel_spectrum(file, "Si IV 1403", xpix=100, ypix=50)