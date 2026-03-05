
import glob
import os
import numpy as np
from astropy.io import fits
import matplotlib as mpl
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from scipy.constants import c
import matplotlib.pyplot as plt
import asdf
import datetime as dt

event = '20230503_072923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

asdf_files = glob.glob(os.path.join(output_loc, "*.asdf"))

if not asdf_files:
    raise FileNotFoundError(f"No ASDF files found in {output_loc}")

af = asdf.open(asdf_files[3])
mg_integ_int = af.tree['mg_integ_int']
def plot_mgii_quartiles(mg_integ_int, aspect_ratio, iris_window, event):
    fig = plt.figure(constrained_layout=True, figsize=(9, 9))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=3, ncols=3, hspace=0.1, wspace=0.01)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.96)

    dopp_rng = 20
    max_wid = 1.0
    asym_rng = 0.5
    min_rat = 0.7
    max_rat = 1.7

    # k integrated intensity map
    ax00 = fig.add_subplot(gs[0, 0], projection=mg_integ_int, label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(mg_integ_int.data, 100 - alpha)

    norm = colors.Normalize(vmin=0, vmax=upr_bnd)
    mg_integ_int.plot_settings['norm'] = norm
    mg_integ_int.plot(axes=ax00, cmap=mpl.colormaps['Reds_r'], title='', aspect=aspect_ratio)
    ax00.set_ylabel("Solar Y (arcsec)")
    ax00.set_xlabel(" ")
    x = ax00.coords[0]
    x.set_ticklabel_visible(False)

    plt.colorbar(location='top', label=r'a) Mg II k Intensity', shrink=0.6, ax=ax00)

plot_mgii_quartiles()
