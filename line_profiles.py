#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%#

import matplotlib.pyplot as plt
import numpy as np
from fit_iris_lines import fit_raster

def plot_pixel_spectrum(file, iris_window, xpix, ypix):
