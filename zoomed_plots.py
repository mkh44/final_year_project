#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%
import glob
import os
from iris_fitting.fit_iris_lines import fit_raster
from iris_fitting import iris_get_mg_features_lv2 as get_mg
from iris_fitting import get_mgii_quartiles
import asdf
from astropy.io import fits
import datetime as dt
from datetime import timedelta
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from iris_fitting import extract_irisL2data
from sunpy.net import Fido, attrs as a
from astropy import units as u
import tarfile
from scipy.constants import speed_of_light


c_1334 = glob.glob('/mnt/nas/ug/hurlem24/iris_data/iris_output/iris_output/20230503_072923/IRIS_fitting_C_II_1334_20230503_072923.asdf')

with asdf.open(c_1334[0]) as af:
    int_1334 = af.tree['int_map']
    dopp_1334 = af.tree['dopp_map']
    width_1334 = af.tree['width_map']
    vnt_1334 = af.tree['vnt_map']
    asym_1334 = af.tree['asym_map']


# Define time frame
t_start = 6500    # seconds from raster start
t_end   = 7000    # seconds from raster start


# Plot the output of the fitting routine if fitting sit-and-stare data
def plot_iris_sns_fits(int_map,dopp_map,width_map,vnt_map,asym_map,iris_window,event,main_header):

    fig = plt.figure(constrained_layout=True, figsize=(10, 10))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=5, ncols=1, hspace=0.05, wspace=0.05)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.95)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

    dopp_rng = 10
    max_wid = 0.1
    asym_rng = 1
    max_vnt = 30

# Slice raster maps in time
    int_data   = int_data[:, t_inds]
    dopp_data  = dopp_map.data[:, t_inds]
    width_data = width_map.data[:, t_inds]
    vnt_data   = vnt_map.data[:, t_inds]
    asym_data  = asym_map.data[:, t_inds]

    t_plot = t_array[t_inds]

# Set the plotting parameters
    cadence = main_header['STEPT_AV']
    t_array = np.arange(0, int_data.shape[1])*cadence
    t_mask = (t_array >= t_start) & (t_array <= t_end)
    t_inds = np.where(t_mask)[0]
    slit_pos = int_map.meta['crval2'] + int_map.meta['cdelt2'] * (np.arange(int_data.shape[0]) - int_map.meta['crpix2'])

# Intensity map
    ax1 = fig.add_subplot(gs[0,0], label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(int_data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    
    plt.imshow(int_data, norm=norm, cmap = mpl.colormaps['Reds_r'], axes=ax1, 
               extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax1.set_ylabel(" ")
    ax1.set_xlabel(" ")
    ax1.set_xticklabels([])

    plt.colorbar(location='right', label=r'a) Intensity', shrink=0.6, ax = ax1, ticks=[0,upr_bnd])

# Asymmetry map
    ax2 = fig.add_subplot(gs[1,0], label='b)')
    norm = colors.Normalize(vmin = -asym_rng, vmax = asym_rng)
    plt.imshow(asym_map.data, norm=norm, cmap = mpl.colormaps['seismic'], axes=ax2, 
               extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax2.set_ylabel(" ")
    ax2.set_xlabel(" ")
    ax2.set_xticklabels([])

    plt.colorbar(location='right', label=r'b) RB Asym.', shrink=0.6, ax = ax2)

# Doppler map
    ax3 = fig.add_subplot(gs[2,0], label='c)')
    norm = colors.Normalize(vmin = -dopp_rng, vmax = dopp_rng)
    plt.imshow(dopp_map.data, norm=norm, cmap = mpl.colormaps['coolwarm'], axes=ax3, 
               extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax3.set_ylabel("Solar Y (arcsec)")
    ax3.set_xlabel(" ")
    ax3.set_xticklabels([])

    plt.colorbar(location='right', label=r'c) v$_{dopp}$ ($km~s^{-1}$)', shrink=0.6, ax = ax3)

# Line width
    ax4 = fig.add_subplot(gs[3,0], label='d)')
    norm = colors.Normalize(vmin = 0, vmax = max_wid)
    plt.imshow(width_map.data, norm=norm, cmap = mpl.colormaps['cubehelix'], axes=ax4, 
               extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax4.set_ylabel(" ")
    ax4.set_xlabel(" ")
    ax4.set_xticklabels([])

    plt.colorbar(location='right', label=r'd) Width ($\AA$)', shrink=0.6, ax = ax4)
    
# Nonthermal velocity
    ax5 = fig.add_subplot(gs[4,0], label='e)')
    norm = colors.Normalize(vmin = 0, vmax = max_vnt)
    plt.imshow(vnt_map.data, norm=norm, cmap = mpl.colormaps['inferno'], axes=ax5, 
               extent=[t_plot.min(), t_plot.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax5.set_ylabel(" ")
    ax5.set_xlabel("Time from raster start (s)")

    plt.colorbar(location='right', label=r'e) v$_{nt}$ ($km~s^{-1}$)', shrink=0.6, ax = ax5)

    plt.suptitle(iris_window+r'$\AA$; '+plot_time)
    plt.savefig(output_loc+event+'/IRIS_analysis_'+iris_window.replace(' ', '_')+'_'+file_time+'.png', bbox_inches='tight')
    plt.close(fig)


            