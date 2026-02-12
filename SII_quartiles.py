# %% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%
import glob
import os
from fit_iris_lines import fit_raster
from import iris_get_mg_features_lv2 as get_mg
from iris_fitting import get_quartiles
import asdf
from astropy.io import fits
import datetime as dt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from iris_fitting import extract_irisL2data
from scipy.constants import speed_of_light
import pdb