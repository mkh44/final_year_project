#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%
import glob
import os
import matplotlib as mpl
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
import pdb

if os.uname().sysname == 'Darwin':
    IRIS_data_loc = '/mnt/nas/ug/hurlem24/iris_data/iris_input_data/iris_l2_20503_072923_4204700135_raster.tar.gz'
    output_loc = '/mnt/nas/ug/hurlem24/iris_data/iris_output/iris_output/'
else:
    IRIS_data_loc = '/mnt/nas/ug/hurlem24/iris_data/iris_input_data/iris_l2_20503_072923_4204700135_raster.tar.gz'
    output_loc = 'C:/Users/molly/OneDrive - Dublin City University/PHA4/Finalar_Project'


# Plot the output of the fitting routine
def plot_iris_fits(int_map,dopp_map,width_map,vnt_map,asym_map,aspect_ratio,iris_window,event):

    fig = plt.figure(constrained_layout=True, figsize=(10, 6))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=1, ncols=5, hspace=0.1, wspace=0.1)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.99)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

# Intensity map
    ax1 = fig.add_subplot(gs[0,0], projection=int_map, label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(int_map.data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    int_map.plot_settings['norm'] = norm
    int_map.plot(axes=ax1, cmap = mpl.colormaps['Reds_r'], title='', aspect=aspect_ratio)
    ax1.set_ylabel("Solar Y (arcsec)")
    ax1.set_xlabel("Solar X (arcsec)")
    x = ax1.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    plt.colorbar(location='top', label=r'a) Intensity', shrink=0.6, ax = ax1, ticks=[0,upr_bnd])

# Asymmetry map
    ax2 = fig.add_subplot(gs[0,1], projection=asym_map, label='b)')
    norm = colors.Normalize(vmin = -5, vmax = 5)
    asym_map.plot_settings['norm'] = norm
    asym_map.plot(axes=ax2, cmap = mpl.colormaps['seismic'], title='', aspect=aspect_ratio)
    ax2.tick_params('x', rotation=30)
    ax2.set_xlabel("Solar X (arcsec)")
    ax2.set_ylabel(" ")
    y = ax2.coords[1]
    y.set_ticklabel_visible(False)
    x = ax2.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    
    plt.colorbar(location='top', label=r'b) RB Asym.', shrink=0.6, ax = ax2)

# Doppler map
    ax3 = fig.add_subplot(gs[0,2], projection=dopp_map, label='c)')
    norm = colors.Normalize(vmin = -20, vmax = 20)
    dopp_map.plot_settings['norm'] = norm
    dopp_map.plot(axes=ax3, cmap = mpl.colormaps['coolwarm'], title='', aspect=aspect_ratio)
    ax3.tick_params('x', rotation=30)
    ax3.set_xlabel("Solar X (arcsec)")
    ax3.set_ylabel(" ")
    y = ax3.coords[1]
    y.set_ticklabel_visible(False)
    x = ax3.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    
    plt.colorbar(location='top', label=r'c) v$_{dopp}$ ($km~s^{-1}$)', shrink=0.6, ax = ax3)

# Line width
    ax4 = fig.add_subplot(gs[0,3], projection=width_map, label='d)')
    norm = colors.Normalize(vmin = 0, vmax = 0.2)
    width_map.plot_settings['norm'] = norm
    width_map.plot(axes=ax4, cmap = mpl.colormaps['cubehelix_r'], title='', aspect=aspect_ratio)
    ax4.tick_params('x', rotation=30)
    ax4.set_xlabel("Solar X (arcsec)")
    ax4.set_ylabel(" ")
    y = ax4.coords[1]
    y.set_ticklabel_visible(False)
    x = ax4.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    
    plt.colorbar(location='top', label=r'd) Width ($\AA$)', shrink=0.6, ax = ax4)
    
# Nonthermal velocity
    ax5 = fig.add_subplot(gs[0,4], projection=vnt_map, label='e)')
    norm = colors.Normalize(vmin = 0, vmax = 50)
    vnt_map.plot_settings['norm'] = norm
    vnt_map.plot(axes=ax5, cmap = mpl.colormaps['cubehelix_r'], title='', aspect=aspect_ratio)
    ax5.tick_params('x', rotation=30)
    ax5.set_xlabel("Solar X (arcsec)")
    ax5.set_ylabel(" ")
    y = ax5.coords[1]
    y.set_ticklabel_visible(False)
    x = ax5.coords[0]
    x.set_ticklabel(exclude_overlapping=True)
    
    plt.colorbar(location='top', label=r'e) v$_{nt}$ ($km~s^{-1}$)', shrink=0.6, ax = ax5)

    plt.suptitle(iris_window+r'$\AA$; '+plot_time)
    plt.xticks(rotation=30, ha='right')
    plt.savefig(output_loc+event+'/IRIS_analysis_'+iris_window.replace(' ', '_')+'_'+file_time+'.png', bbox_inches='tight')
    plt.close(fig)

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

# Set the plotting parameters
    cadence = main_header['STEPT_AV']
    t_array = np.arange(0, int_map.data.shape[1])*cadence
    slit_pos = int_map.meta['crval2'] + int_map.meta['cdelt2'] * (np.arange(int_map.data.shape[0]) - int_map.meta['crpix2'])

# Intensity map
    ax1 = fig.add_subplot(gs[0,0], label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(int_map.data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    
    plt.imshow(int_map.data, norm=norm, cmap = mpl.colormaps['Reds_r'], axes=ax1, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax1.set_ylabel(" ")
    ax1.set_xlabel(" ")
    ax1.set_xticklabels([])

    plt.colorbar(location='right', label=r'a) Intensity', shrink=0.6, ax = ax1, ticks=[0,upr_bnd])

# Asymmetry map
    ax2 = fig.add_subplot(gs[1,0], label='b)')
    norm = colors.Normalize(vmin = -asym_rng, vmax = asym_rng)
    plt.imshow(asym_map.data, norm=norm, cmap = mpl.colormaps['seismic'], axes=ax2, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax2.set_ylabel(" ")
    ax2.set_xlabel(" ")
    ax2.set_xticklabels([])

    plt.colorbar(location='right', label=r'b) RB Asym.', shrink=0.6, ax = ax2)

# Doppler map
    ax3 = fig.add_subplot(gs[2,0], label='c)')
    norm = colors.Normalize(vmin = -dopp_rng, vmax = dopp_rng)
    plt.imshow(dopp_map.data, norm=norm, cmap = mpl.colormaps['coolwarm'], axes=ax3, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax3.set_ylabel("Solar Y (arcsec)")
    ax3.set_xlabel(" ")
    ax3.set_xticklabels([])

    plt.colorbar(location='right', label=r'c) v$_{dopp}$ ($km~s^{-1}$)', shrink=0.6, ax = ax3)

# Line width
    ax4 = fig.add_subplot(gs[3,0], label='d)')
    norm = colors.Normalize(vmin = 0, vmax = max_wid)
    plt.imshow(width_map.data, norm=norm, cmap = mpl.colormaps['cubehelix'], axes=ax4, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax4.set_ylabel(" ")
    ax4.set_xlabel(" ")
    ax4.set_xticklabels([])

    plt.colorbar(location='right', label=r'd) Width ($\AA$)', shrink=0.6, ax = ax4)
    
# Nonthermal velocity
    ax5 = fig.add_subplot(gs[4,0], label='e)')
    norm = colors.Normalize(vmin = 0, vmax = max_vnt)
    plt.imshow(vnt_map.data, norm=norm, cmap = mpl.colormaps['inferno'], axes=ax5, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax5.set_ylabel(" ")
    ax5.set_xlabel("Time from raster start (s)")

    plt.colorbar(location='right', label=r'e) v$_{nt}$ ($km~s^{-1}$)', shrink=0.6, ax = ax5)

    plt.suptitle(iris_window+r'$\AA$; '+plot_time)
    plt.savefig(output_loc+event+'/IRIS_analysis_'+iris_window.replace(' ', '_')+'_'+file_time+'.png', bbox_inches='tight')
    plt.close(fig)
    
# Plot the output of the Mg II fitting routine
def plot_mgii_fits(dv_k3_map,dv_h3_map,k2_sep_map,h2_sep_map,aspect_ratio,iris_window,event):

    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=1, ncols=4, hspace=0.1, wspace=0.01)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.99)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(dv_k3_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(dv_k3_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

# dv k3 map
    ax1 = fig.add_subplot(gs[0,0], projection=dv_k3_map, label='a)')
    norm = colors.Normalize(vmin = -20, vmax = 20)
    dv_k3_map.plot_settings['norm'] = norm
    dv_k3_map.plot(axes=ax1, cmap = mpl.colormaps['coolwarm'], title='', aspect=aspect_ratio)
    ax1.tick_params('x', labelrotation=30)
    ax1.set_ylabel("Solar Y (arcsec)")
    ax1.set_xlabel("Solar X (arcsec)")

    plt.colorbar(location='top', label=r'a) $\Delta v_{k3}$ (km/s)', shrink=0.6, ax = ax1)

# dv h3 map
    ax2 = fig.add_subplot(gs[0,1], projection=dv_h3_map, label='b)')
    norm = colors.Normalize(vmin = -20, vmax = 20)
    dv_h3_map.plot_settings['norm'] = norm
    dv_h3_map.plot(axes=ax2, cmap = mpl.colormaps['coolwarm'], title='', aspect=aspect_ratio)
    ax2.tick_params('x', labelrotation=30)
    ax2.set_xlabel("Solar X (arcsec)")
    ax2.set_ylabel(" ")
    y = ax2.coords[1]
    y.set_ticklabel_visible(False)
    
    plt.colorbar(location='top', label=r'b) $\Delta v_{h3}$ (km/s)', shrink=0.6, ax = ax2)

# k2 sep map
    ax3 = fig.add_subplot(gs[0,2], projection=k2_sep_map, label='c)')
    norm = colors.Normalize(vmin = 0, vmax = 50)
    k2_sep_map.plot_settings['norm'] = norm
    k2_sep_map.plot(axes=ax3, cmap = mpl.colormaps['cubehelix_r'], title='', aspect=aspect_ratio)
    ax3.tick_params('x', labelrotation=30)
    ax3.set_xlabel("Solar X (arcsec)")
    ax3.set_ylabel(" ")
    y = ax3.coords[1]
    y.set_ticklabel_visible(False)
    
    plt.colorbar(location='top', label=r'c) $\Delta v_{k2}$ (km/s)', shrink=0.6, ax = ax3)

# h2 sep map
    ax4 = fig.add_subplot(gs[0,3], projection=h2_sep_map, label='d)')
    norm = colors.Normalize(vmin = 0, vmax = 50)
    h2_sep_map.plot_settings['norm'] = norm
    h2_sep_map.plot(axes=ax4, cmap = mpl.colormaps['cubehelix_r'], title='', aspect=aspect_ratio)
    ax4.tick_params('x', labelrotation=30)
    ax4.set_xlabel("Solar X (arcsec)")
    ax4.set_ylabel(" ")
    y = ax4.coords[1]
    y.set_ticklabel_visible(False)
    
    plt.colorbar(location='top', label=r'd) $\Delta v_{h2}$ (km/s)', shrink=0.6, ax = ax4)
    
    plt.suptitle(iris_window+'; '+plot_time)

    plt.savefig(output_loc+event+'/IRIS_analysis_'+iris_window.replace(' ', '_')+'_'+file_time+'.png', bbox_inches='tight')
    plt.close(fig)

# Plot the output of the Mg II fitting routine
def plot_mgii_sns_fits(dv_k3_map,dv_h3_map,k2_sep_map,h2_sep_map,iris_window,event,main_header):

    fig = plt.figure(constrained_layout=True, figsize=(10, 10))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=4, ncols=1, hspace=0.05, wspace=0.01)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.95)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(dv_k3_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(dv_k3_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

# Set the plotting parameters
    cadence = main_header['STEPT_AV']
    t_array = np.arange(0, dv_k3_map.data.shape[1])*cadence
    slit_pos = dv_k3_map.meta['crval2'] + dv_k3_map.meta['cdelt2'] * (np.arange(dv_k3_map.data.shape[0]) - dv_k3_map.meta['crpix2'])

# dv k3 map
    ax1 = fig.add_subplot(gs[0,0], label='a)')
    norm = colors.Normalize(vmin = -20, vmax = 20)
    plt.imshow(dv_k3_map.data.T, norm=norm, cmap = mpl.colormaps['coolwarm'], origin='lower', axes=ax1, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax1.set_ylabel("Solar Y (arcsec)")
    ax1.set_xlabel(" ")
    ax1.set_xticklabels([])
    plt.colorbar(location='right', label=r'a) $\Delta v_{k3}$ (km/s)', shrink=0.6, ax = ax1)

# dv h3 map
    ax2 = fig.add_subplot(gs[1,0], label='b)')
    norm = colors.Normalize(vmin = -20, vmax = 20)
    plt.imshow(dv_h3_map.data.T, norm=norm, cmap = mpl.colormaps['coolwarm'], origin='lower', axes=ax2, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax2.set_ylabel("Solar Y (arcsec)")
    ax2.set_xlabel(" ")
    ax2.set_xticklabels([])
    plt.colorbar(location='right', label=r'b) $\Delta v_{h3}$ (km/s)', shrink=0.6, ax = ax2)

# k2 sep map
    ax3 = fig.add_subplot(gs[2,0], label='c)')
    norm = colors.Normalize(vmin = 0, vmax = 50)
    plt.imshow(k2_sep_map.data.T, norm=norm, cmap = mpl.colormaps['cubehelix_r'], origin='lower', axes=ax3, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax3.set_ylabel("Solar Y (arcsec)")
    ax3.set_xlabel(" ")
    ax3.set_xticklabels([])
    plt.colorbar(location='right', label=r'c) $\Delta v_{k2}$ (km/s)', shrink=0.6, ax = ax3)

# h2 sep map
    ax4 = fig.add_subplot(gs[3,0], label='d)')
    norm = colors.Normalize(vmin = 0, vmax = 50)
    plt.imshow(h2_sep_map.data.T, norm=norm, cmap = mpl.colormaps['cubehelix_r'], origin='lower', axes=ax4, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax4.set_ylabel("Solar Y (arcsec)")
    ax4.set_xlabel("Time from raster start (s)")

    plt.colorbar(location='right', label=r'd) $\Delta v_{h2}$ (km/s)', shrink=0.6, ax = ax4)
    
    plt.suptitle('Mg II analysis; '+plot_time)

    plt.savefig(output_loc+event+'/IRIS_analysis_'+iris_window.replace(' ', '_')+'_'+file_time+'.png', bbox_inches='tight')
    plt.close(fig)

# Plot the output of the Mg II fitting routine
def plot_mgii_quartiles(mgii_k_integ_int,mgii_h_integ_int,mgii_k_vdopp,mgii_h_vdopp,mgii_k_width,mgii_h_width,
                        mgii_k_asym,mgii_h_asym,mg_integ_int,aspect_ratio,iris_window,event,smooth=False):

    fig = plt.figure(constrained_layout=True, figsize=(9, 9))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=3, ncols=3, hspace=0.1, wspace=0.01)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.96)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(mgii_k_vdopp.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(mgii_k_vdopp.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

    dopp_rng = 20
    max_wid = 1.5
    asym_rng = 0.5
    min_rat = 0.7
    max_rat = 1.7

# k integrated intensity map
    ax00 = fig.add_subplot(gs[0,0], projection=mgii_k_integ_int, label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(mgii_k_integ_int.data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    mgii_k_integ_int.plot_settings['norm'] = norm
    mgii_k_integ_int.plot(axes=ax00, cmap = mpl.colormaps['Reds_r'], title='', aspect=aspect_ratio)
    ax00.set_ylabel("Solar Y (arcsec)")
    ax00.set_xlabel(" ")
    x = ax00.coords[0]
    x.set_ticklabel_visible(False)

    plt.colorbar(location='top', label=r'a) Mg II k Intensity', shrink=0.6, ax = ax00)

# k doppler map
    ax01 = fig.add_subplot(gs[0,1], projection=mgii_k_vdopp, label='b)')
    norm = colors.Normalize(vmin = -dopp_rng, vmax = dopp_rng)
    mgii_k_vdopp.plot_settings['norm'] = norm
    mgii_k_vdopp.plot(axes=ax01, cmap = mpl.colormaps['seismic'], title='', aspect=aspect_ratio)
    ax01.set_ylabel(" ")
    ax01.set_xlabel(" ")
    x = ax01.coords[0]
    x.set_ticklabel_visible(False)

    plt.colorbar(location='top', label=r'b) Mg II k V$_{dopp}$ (km/s)', shrink=0.6, ax = ax01)

# k line width map
    ax02 = fig.add_subplot(gs[0,2], projection=mgii_k_width, label='c)')
    norm = colors.Normalize(vmin = 0, vmax = max_wid)
    mgii_k_width.plot_settings['norm'] = norm
    mgii_k_width.plot(axes=ax02, cmap = mpl.colormaps['cubehelix'], title='', aspect=aspect_ratio)
    ax02.set_xlabel(" ")
    ax02.set_ylabel(" ")
    x = ax02.coords[0]
    x.set_ticklabel_visible(False)
    y = ax02.coords[1]
    y.set_ticklabel_visible(False)
    
    plt.colorbar(location='top', label=r'c) Mg II k line width ($\AA$)', shrink=0.6, ax = ax02)

# k asym map
    ax10 = fig.add_subplot(gs[1,0], projection=mgii_k_asym, label='d)')
    norm = colors.Normalize(vmin = -asym_rng, vmax = asym_rng)
    mgii_k_asym.plot_settings['norm'] = norm
    mgii_k_asym.plot(axes=ax10, cmap = mpl.colormaps['bwr'], title='', aspect=aspect_ratio)
    ax10.set_xlabel(" ")
    ax10.set_ylabel("Solar Y (arcsec)")
    x = ax10.coords[0]
    x.set_ticklabel_visible(False)
    
    plt.colorbar(location='top', label=r'd) Mg II k asymmetry', shrink=0.6, ax = ax10)

# h integrated intensity map
    ax11 = fig.add_subplot(gs[1,1], projection=mgii_h_integ_int, label='e)')
    alpha = 1
    upr_bnd = np.nanpercentile(mgii_h_integ_int.data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    mgii_h_integ_int.plot_settings['norm'] = norm
    mgii_h_integ_int.plot(axes=ax11, cmap = mpl.colormaps['Reds_r'], title='', aspect=aspect_ratio)
    ax11.set_ylabel(" ")
    ax11.set_xlabel(" ")
    x = ax11.coords[0]
    x.set_ticklabel_visible(False)
    y = ax11.coords[1]
    y.set_ticklabel_visible(False)

    plt.colorbar(location='top', label=r'e) Mg II h Intensity', shrink=0.6, ax = ax11)

# h doppler map
    ax12 = fig.add_subplot(gs[1,2], projection=mgii_h_vdopp, label='f)')
    norm = colors.Normalize(vmin = -dopp_rng, vmax = dopp_rng)
    mgii_h_vdopp.plot_settings['norm'] = norm
    mgii_h_vdopp.plot(axes=ax12, cmap = mpl.colormaps['seismic'], title='', aspect=aspect_ratio)
    ax12.set_ylabel(" ")
    ax12.set_xlabel(" ")
    x = ax12.coords[0]
    x.set_ticklabel_visible(False)
    y = ax12.coords[1]
    y.set_ticklabel_visible(False)

    plt.colorbar(location='top', label=r'f) Mg II h V$_{dopp}$ (km/s)', shrink=0.6, ax = ax12)

# h line width map
    ax20 = fig.add_subplot(gs[2,0], projection=mgii_h_width, label='g)')
    norm = colors.Normalize(vmin = 0, vmax = max_wid)
    mgii_h_width.plot_settings['norm'] = norm
    mgii_h_width.plot(axes=ax20, cmap = mpl.colormaps['cubehelix'], title='', aspect=aspect_ratio)
    ax20.set_xlabel("Solar X (arcsec)")
    ax20.set_ylabel("Solar y (arcsec)")
    
    plt.colorbar(location='top', label=r'g) Mg II h line width ($\AA$)', shrink=0.6, ax = ax20)

# h asym map
    ax21 = fig.add_subplot(gs[2,1], projection=mgii_h_asym, label='h)')
    norm = colors.Normalize(vmin = -asym_rng, vmax = asym_rng)
    mgii_h_asym.plot_settings['norm'] = norm
    mgii_h_asym.plot(axes=ax21, cmap = mpl.colormaps['bwr'], title='', aspect=aspect_ratio)
    ax21.set_xlabel("Solar X (arcsec)")
    ax21.set_ylabel(" ")
    y = ax21.coords[1]
    y.set_ticklabel_visible(False)
    
    plt.colorbar(location='top', label=r'h) Mg II h asymmetry', shrink=0.6, ax = ax21)

# Intensity ratio map
    ax22 = fig.add_subplot(gs[2,2], projection=mg_integ_int, label='i)')
    norm = colors.Normalize(vmin = min_rat, vmax = max_rat)
    mg_integ_int.plot_settings['norm'] = norm
    mg_integ_int.plot(axes=ax22, cmap = mpl.colormaps['cubehelix_r'], title='', aspect=aspect_ratio)
    ax22.set_xlabel("Solar X (arcsec)")
    ax22.set_ylabel(" ")
    y = ax22.coords[1]
    y.set_ticklabel_visible(False)

    plt.colorbar(location='top', label=r'i) k/h ratio', shrink=0.6, ax = ax22)


    plt.suptitle(iris_window+'; '+plot_time)

    if smooth:
        filename = output_loc+event+'/IRIS_quartiles_smooth_'+iris_window.replace(' ', '_')+'_'+file_time+'.png'
    else:
        filename = output_loc+event+'/IRIS_quartiles_'+iris_window.replace(' ', '_')+'_'+file_time+'.png'

    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)

# Plot the output of the Mg II fitting routine
def plot_mgii_sns_quartiles(mgii_k_integ_int,mgii_h_integ_int,mgii_k_vdopp,mgii_h_vdopp,mgii_k_width,mgii_h_width,
                            mgii_k_asym,mgii_h_asym,mg_integ_int,iris_window,event,main_header,smooth=False):

    fig = plt.figure(constrained_layout=True, figsize=(18, 12))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=3, ncols=3, hspace=0.1, wspace=0.05)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.95)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(mgii_k_vdopp.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(mgii_k_vdopp.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

# Set the plotting parameters
    cadence = main_header['STEPT_AV']
    t_array = np.arange(0, mgii_k_vdopp.data.shape[1])*cadence
    slit_pos = mgii_k_vdopp.meta['crval2'] + mgii_k_vdopp.meta['cdelt2'] * (np.arange(mgii_k_vdopp.data.shape[0]) - mgii_k_vdopp.meta['crpix2'])

    dopp_rng = 30
    max_wid = 1.5
    asym_rng = 1
    min_rat = 0.7
    max_rat = 1.7

# k integrated intensity map
    ax00 = fig.add_subplot(gs[0,0], label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(mgii_k_integ_int.data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    plt.imshow(mgii_k_integ_int.data, norm=norm, cmap = mpl.colormaps['Reds_r'], axes=ax00, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax00.set_ylabel("Solar Y (arcsec)")
    ax00.set_xlabel(" ")
    ax00.set_xticklabels([])

    plt.colorbar(location='top', label=r'Mg II k Intensity', shrink=0.9, ax = ax00)

# k doppler map
    ax01 = fig.add_subplot(gs[0,1], label='b)')
    norm = colors.Normalize(vmin = -dopp_rng, vmax = dopp_rng)
    plt.imshow(mgii_k_vdopp.data, norm=norm, cmap = mpl.colormaps['seismic'], axes=ax01, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax01.set_ylabel(" ")
    ax01.set_xlabel(" ")
    ax01.set_xticklabels([])
    ax01.set_yticklabels([])

    plt.colorbar(location='top', label=r'Mg II k v$_{dopp}$ (km/s)', shrink=0.9, ax = ax01)

# k line width map
    ax02 = fig.add_subplot(gs[0,2], label='c)')
    norm = colors.Normalize(vmin = 0, vmax = max_wid)
    plt.imshow(mgii_k_width.data, norm=norm, cmap = mpl.colormaps['cubehelix'], axes=ax02, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax02.set_ylabel(" ")
    ax02.set_xlabel(" ")
    ax02.set_xticklabels([])
    ax02.set_yticklabels([])

    plt.colorbar(location='top', label=r'Mg II k Line width ($\AA$)', shrink=0.9, ax = ax02)

# k asymmetry map
    ax10 = fig.add_subplot(gs[1,0], label='d)')
    norm = colors.Normalize(vmin = -asym_rng, vmax = asym_rng)
    plt.imshow(mgii_k_asym.data, norm=norm, cmap = mpl.colormaps['coolwarm'], axes=ax10, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax10.set_ylabel("Solar Y (arcsec)")
    ax10.set_xlabel(" ")
    ax10.set_xticklabels([])

    plt.colorbar(location='top', label=r'Mg II k Asymmetry', shrink=0.9, ax = ax10)

# h integrated intensity map
    ax11 = fig.add_subplot(gs[1,1], label='e)')
    alpha = 1
    upr_bnd = np.nanpercentile(mgii_h_integ_int.data, 100-alpha)

    norm = colors.Normalize(vmin = 0, vmax = upr_bnd)
    plt.imshow(mgii_h_integ_int.data, norm=norm, cmap = mpl.colormaps['Reds_r'], axes=ax11, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax11.set_ylabel(" ")
    ax11.set_xlabel(" ")
    ax11.set_xticklabels([])
    ax11.set_yticklabels([])

    plt.colorbar(location='top', label=r'Mg II h Intensity', shrink=0.9, ax = ax11)

# h doppler map
    ax12 = fig.add_subplot(gs[1,2], label='f)')
    norm = colors.Normalize(vmin = -dopp_rng, vmax = dopp_rng)
    plt.imshow(mgii_h_vdopp.data, norm=norm, cmap = mpl.colormaps['seismic'], axes=ax12, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax12.set_ylabel(" ")
    ax12.set_xlabel(" ")
    ax12.set_xticklabels([])
    ax12.set_yticklabels([])

    plt.colorbar(location='top', label=r'Mg II h v$_{dopp}$ (km/s)', shrink=0.9, ax = ax12)

# h line width map
    ax20 = fig.add_subplot(gs[2,0], label='g)')
    norm = colors.Normalize(vmin = 0, vmax = max_wid)
    plt.imshow(mgii_h_width.data, norm=norm, cmap = mpl.colormaps['cubehelix'], axes=ax20, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax20.set_ylabel("Solar Y (arcsec)")
    ax20.set_xlabel("Time from raster start (s)")

    plt.colorbar(location='top', label=r'Mg II h Line width ($\AA$)', shrink=0.9, ax = ax20)

# h asymmetry map
    ax21 = fig.add_subplot(gs[2,1], label='h)')
    norm = colors.Normalize(vmin = -asym_rng, vmax = asym_rng)
    plt.imshow(mgii_h_asym.data, norm=norm, cmap = mpl.colormaps['coolwarm'], axes=ax21, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax21.set_ylabel(" ")
    ax21.set_xlabel("Time from raster start (s)")
    ax21.set_yticklabels([])

    plt.colorbar(location='top', label=r'Mg II h Asymmetry', shrink=0.9, ax = ax21)

# Intensity ratio map
    ax22 = fig.add_subplot(gs[2,2], label='i)')
    norm = colors.Normalize(vmin = min_rat, vmax = max_rat)
    plt.imshow(mg_integ_int.data, norm=norm, cmap = mpl.colormaps['cubehelix_r'], axes=ax22, 
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax22.set_xlabel("Time from raster start (s)")
    ax22.set_ylabel(" ")
    ax22.set_yticklabels([])

    plt.colorbar(location='top', label=r'Mg II k/h ratio', shrink=0.9, ax = ax22)

    plt.suptitle('Raster start time = '+plot_time)

    if smooth:
        filename = output_loc+event+'/IRIS_quartiles_smooth_'+iris_window.replace(' ', '_')+'_'+file_time+'.png'
    else:
        filename = output_loc+event+'/IRIS_quartiles_'+iris_window.replace(' ', '_')+'_'+file_time+'.png'

    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)

# Fit the Si IV spectral lines
def fit_iris(file, iris_window, event, do_fit=False):

# Check the wavelength window
    match iris_window:
        case "Si IV 1394 1403":
            open_window='Si IV 1403'
            save_window = 'Si IV 1403'
        case "C II 1334 1336":
            open_window='C II 1336'
            save_window = 'C II 1334'
        case "C II 1335 1336":
            open_window='C II 1336'
            save_window = 'C II 1335'
        case _:
            open_window = iris_window
            save_window = iris_window

    a = fit_raster(file, iris_window, fulldisk=False)

    if do_fit:
        results_array, int_map, dopp_map, width_map, vnt_map, asym_map = a.fit_iris_data(v_nontherm=True)
        plot_time = int_map.date.strftime('%Y%m%d_%H%M%S')
# Save the outputs as an asdf file for each run
        tree = {'results_array':results_array, 'int_map':int_map, 'dopp_map':dopp_map, 'width_map':width_map, 
                'vnt_map':vnt_map, 'asym_map':asym_map}    
        with asdf.AsdfFile(tree) as asdf_file:  
            asdf_file.write_to(output_loc+event+'/IRIS_fitting_'+save_window.replace(' ', '_')+'_'+plot_time+'.asdf',
                            all_array_compression='zlib')
    else:
        data, main_header, header, wavelength = a.open_iris_file(open_window)
        img_time = dt.datetime.strptime(main_header['date_obs'], '%Y-%m-%dT%H:%M:%S.%f').strftime('%Y%m%d_%H%M%S')
        with asdf.open(output_loc+event+'/IRIS_fitting_'+save_window.replace(' ', '_')+'_'+img_time+'.asdf') as af:
            int_map = af.tree['int_map']
            dopp_map = af.tree['dopp_map']
            width_map = af.tree['width_map']
            vnt_map = af.tree['vnt_map']
            asym_map = af.tree['asym_map']

# Plot the outputs
    if int_map.meta['cdelt1'] == 1e-9:
        sp = fits.open(file)
        main_header = sp[0].header
        plot_iris_sns_fits(int_map,dopp_map,width_map,vnt_map,asym_map,save_window,event,main_header)
    else:
        aspect_ratio = np.abs(int_map.meta['cdelt2']) / np.abs(int_map.meta['cdelt1'])
        plot_iris_fits(int_map,dopp_map,width_map,vnt_map,asym_map,aspect_ratio,save_window,event)

    return

# Get the Mg II line properties
def fit_iris_mgii(file, iris_window, event, plot_time, do_fit=False, smooth=False):

    if smooth:
        int_25 = 3
        int_50 = 4
        int_75 = 5
        filename = output_loc+event+'/IRIS_fitting_smooth_MgII_'+plot_time+'.asdf'
    else:
        int_25 = 0
        int_50 = 1
        int_75 = 2
        filename = output_loc+event+'/IRIS_fitting_MgII_'+plot_time+'.asdf'

    a = fit_raster(file, iris_window)
    sp = fits.open(file, memmap=False)
    main_header = sp[0].header
    lines = []
    for i in range(main_header['NWIN']):
        lines.append(main_header['TDESC{}'.format(i+1)])
    ind = lines.index(iris_window)+1
    header = sp[ind].header

    if do_fit:
        # First, get the Mg II line properties using the iris_get_mg_features_lev2 code
        lc, rp, bp = get_mg.iris_get_mg_features_lev2(file)

        lc_k_v = lc[0,:,:,0]
        lc_k_i = lc[0,:,:,1]
        lc_h_v = lc[1,:,:,0]
        lc_h_i = lc[1,:,:,1]

        rp_k_v = rp[0,:,:,0]
        rp_k_i = rp[0,:,:,1]
        rp_h_v = rp[1,:,:,0]
        rp_h_i = rp[1,:,:,1]

        bp_k_v = bp[0,:,:,0]
        bp_k_i = bp[0,:,:,1]
        bp_h_v = bp[1,:,:,0]
        bp_h_i = bp[1,:,:,1]
    
        dv_k3_map = a.mk_iris_map(lc_k_v, header, main_header) # Upper chromosphere velocity

        dv_h3_map = a.mk_iris_map(lc_h_v, header, main_header) # Upper chromosphere velocity

        k2_diff = rp_k_v - bp_k_v
        k2_sep_map = a.mk_iris_map(k2_diff, header, main_header) # Mid chromosphere velocity

        h2_diff = rp_h_v - bp_h_v
        h2_sep_map = a.mk_iris_map(h2_diff, header, main_header) # Mid chromosphere velocity

    # Now get the Mg II line properties using the quartile approach. First do Mg II k
        # Give the rest wavelength to get an estimate of where the line is
        quartiles_k = get_mgii_quartiles.get_mgii_quartiles(file, 'Mg II k')

        y_size = quartiles_k.shape[0]
        x_size = quartiles_k.shape[1]
        line_pos = quartiles_k[:,:,int_50]

        integ_int = quartiles_k[:,:,6]
        if header['CDELT3'] == 0:
            mgii_k_integ_int = a.mk_iris_map(integ_int.T, header, main_header)
        else:
            mgii_k_integ_int = a.mk_iris_map(integ_int, header, main_header)

        # Doppler velocity
        # We're going to take the average median of the datacube as the rest wavelength in the region of interest
        wvl_ref = np.mean(line_pos)
        ref_wave = np.full((y_size, x_size), wvl_ref)
        v_dopp = ((line_pos - ref_wave)/(ref_wave))*(speed_of_light/1e3)
        if header['CDELT3'] == 0:
            mgii_k_vdopp = a.mk_iris_map(v_dopp.T, header, main_header)
        else:
            mgii_k_vdopp = a.mk_iris_map(v_dopp, header, main_header)

        # Line width
        line_width = quartiles_k[:,:,int_75] - quartiles_k[:,:,int_25]
        if header['CDELT3'] == 0:
            mgii_k_width = a.mk_iris_map(line_width.T, header, main_header)
        else:
            mgii_k_width = a.mk_iris_map(line_width, header, main_header)

        # Asymmetry
        asym = (((quartiles_k[:,:,int_75]-quartiles_k[:,:,int_50])-
                 (quartiles_k[:,:,int_50]-quartiles_k[:,:,int_25]))/(quartiles_k[:,:,int_75]-quartiles_k[:,:,int_25]))
        if header['CDELT3'] == 0:
            mgii_k_asym = a.mk_iris_map(asym.T, header, main_header)
        else:
            mgii_k_asym = a.mk_iris_map(asym, header, main_header)

    # Next do the Mg II h line
        quartiles_h = get_mgii_quartiles.get_mgii_quartiles(file, 'Mg II h')

        y_size = quartiles_h.shape[0]
        x_size = quartiles_h.shape[1]
        line_pos = quartiles_h[:,:,int_50]

        integ_int = quartiles_h[:,:,6]
        if header['CDELT3'] == 0:
            mgii_h_integ_int = a.mk_iris_map(integ_int.T, header, main_header)
        else:
            mgii_h_integ_int = a.mk_iris_map(integ_int, header, main_header)

        # Doppler velocity
        # We're going to take the average median of the datacube as the rest wavelength in the region of interest
        wvl_ref = np.mean(line_pos)
        ref_wave = np.full((y_size, x_size), wvl_ref)
        v_dopp = ((line_pos - ref_wave)/(ref_wave))*(speed_of_light/1e3)
        if header['CDELT3'] == 0:
            mgii_h_vdopp = a.mk_iris_map(v_dopp.T, header, main_header)
        else:
            mgii_h_vdopp = a.mk_iris_map(v_dopp, header, main_header)

        # Line width
        line_width = quartiles_h[:,:,int_75] - quartiles_h[:,:,int_25]
        if header['CDELT3'] == 0:
            mgii_h_width = a.mk_iris_map(line_width.T, header, main_header)
        else:
            mgii_h_width = a.mk_iris_map(line_width, header, main_header)

        # Asymmetry
        asym = (((quartiles_h[:,:,int_75]-quartiles_h[:,:,int_50]) - 
                 (quartiles_h[:,:,int_50]-quartiles_h[:,:,int_25]))/(quartiles_h[:,:,int_75]-quartiles_h[:,:,int_25]))
        if header['CDELT3'] == 0:
            mgii_h_asym = a.mk_iris_map(asym.T, header, main_header)
        else:
            mgii_h_asym = a.mk_iris_map(asym, header, main_header)

        # Integrated intensity
        int_ratio = np.divide(mgii_k_integ_int.data, mgii_h_integ_int.data)
        mg_int_rat = a.mk_iris_map(int_ratio, header, main_header)

# Save the output
        tree = {'dv_k3_map':dv_k3_map, 'dv_h3_map':dv_h3_map, 'k2_sep_map':k2_sep_map, 'h2_sep_map':h2_sep_map, 
                'mgii_k_integ_int':mgii_k_integ_int, 'mgii_h_integ_int':mgii_h_integ_int, 'mgii_k_vdopp':mgii_k_vdopp, 
                'mgii_h_vdopp':mgii_h_vdopp, 'mgii_k_width':mgii_k_width, 'mgii_h_width':mgii_h_width, 
                'mgii_k_asym':mgii_k_asym, 'mgii_h_asym':mgii_h_asym, 'mg_int_rat':mg_int_rat}
        with asdf.AsdfFile(tree) as asdf_file:  
            asdf_file.write_to(filename, all_array_compression='zlib')
    
    else:
        # If the data has already been processed, open the asdf file
        with asdf.open(filename) as af:
            dv_k3_map = af.tree['dv_k3_map']
            dv_h3_map = af.tree['dv_h3_map']
            k2_sep_map = af.tree['k2_sep_map']
            h2_sep_map = af.tree['h2_sep_map']
            mgii_k_integ_int = af.tree['mgii_k_integ_int']
            mgii_h_integ_int = af.tree['mgii_h_integ_int']
            mgii_k_vdopp = af.tree['mgii_k_vdopp']
            mgii_h_vdopp = af.tree['mgii_h_vdopp']
            mgii_k_width = af.tree['mgii_k_width']
            mgii_h_width = af.tree['mgii_h_width']
            mgii_k_asym = af.tree['mgii_k_asym']
            mgii_h_asym = af.tree['mgii_h_asym']
            mg_int_rat = af.tree['mg_int_rat']

# Plot the outputs
    if dv_k3_map.meta['cdelt1'] == 1e-9:
        sp = fits.open(file)
        main_header = sp[0].header
        plot_mgii_sns_fits(dv_k3_map,dv_h3_map,k2_sep_map,h2_sep_map,iris_window,event,main_header)
        plot_mgii_sns_quartiles(mgii_k_integ_int,mgii_h_integ_int,mgii_k_vdopp,mgii_h_vdopp,mgii_k_width,
                                mgii_h_width,mgii_k_asym,mgii_h_asym,mg_int_rat,iris_window,event,main_header,
                                smooth=smooth)
    else:
        aspect_ratio = np.abs(dv_k3_map.meta['cdelt2']) / np.abs(dv_k3_map.meta['cdelt1'])
        plot_mgii_fits(dv_k3_map,dv_h3_map,k2_sep_map,h2_sep_map,aspect_ratio,iris_window,event)
        plot_mgii_quartiles(mgii_k_integ_int,mgii_h_integ_int,mgii_k_vdopp,mgii_h_vdopp,mgii_k_width,
                            mgii_h_width,mgii_k_asym,mgii_h_asym,mg_int_rat,aspect_ratio,iris_window,event,
                            smooth=smooth)

    return 

# Main fitting routine
def fitdata(event):

    # Make the output directory
    os.makedirs(output_loc+event+'/', exist_ok='True')

    # Find the files to be processed
    f_iris_raster = glob.glob(IRIS_data_loc+event+"/*.fits")
    f_iris_raster.sort()

    for indiv_file in f_iris_raster:

        # Get the file time
        main_header = extract_irisL2data.only_header(indiv_file,extension=0,verbose=False)
        file_time = dt.datetime.strftime(dt.datetime.strptime(main_header['DATE_OBS'],'%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')
        
        # Define the windows to be processed. Can process both Si IV lines and both C II lines. 
        # Note that C II lines need to be called as: C II 1334 1336 & C II 1335 1336
        iris_window_list = ['C II 1334 1336']#, 'C II 1335 1336', 'Si IV 1394', 'Si IV 1403']

        for iris_window in iris_window_list:

# Check the wavelength window
            match iris_window:
                case 'Si IV 1394 1403': 
                    save_window = 'Si IV 1403'
                case 'C II 1334 1336': 
                    save_window = 'C II 1334'
                case 'C II 1335 1336': 
                    save_window = 'C II 1335'
                case _:
                    save_window = iris_window

            print('Processing '+save_window)

            try:
                files = os.path.exists(output_loc+event+'/IRIS_fitting_'+save_window.replace(' ', '_')+'_'+file_time+'.asdf')
                if files == True:
                    print(save_window+' data already processed. File exists at: '+
                          output_loc+event+'/IRIS_fitting_'+save_window.replace(' ', '_')+'_'+file_time+'.asdf')
                    print("Updating plot...")
                    fit_iris(indiv_file, iris_window, event, do_fit=False)
                if files == False:
                    print('Fitting '+save_window)
                    fit_iris(indiv_file, iris_window, event, do_fit=True)
            except (UnboundLocalError, IndexError):
                print('Wavelength window not available')

#        # Mg II
        files = os.path.exists(output_loc+event+'/IRIS_fitting_MgII_'+file_time+'.asdf')
        if files == True:
            print('Mg II data already processed. File exists at: '+
                output_loc+event+'/IRIS_fitting_MgII_'+file_time+'.asdf')
            print("Updating plot...")
            fit_iris_mgii(indiv_file, 'Mg II k 2796', event, file_time, do_fit=False, smooth=True)
        if files == False:
            print('Fitting Mg II data')
            fit_iris_mgii(indiv_file, 'Mg II k 2796', event, file_time, do_fit=True, smooth=True)

# Download the IRIS data
def getdata(iris_event):
    evt_date = dt.datetime.strptime(iris_event, '%Y%m%d_%H%M%S')
    attrs_time = a.Time(evt_date-timedelta(minutes=5), evt_date+timedelta(minutes=5))
    passband = a.Wavelength(int(1332)*u.Angstrom, int(1332)*u.Angstrom)
    instr = a.Instrument('IRIS')
    source = a.Source('IRIS')

    data_dir = IRIS_data_loc+iris_event+'/'

    results = Fido.search(attrs_time,instr,source,passband)
    
    fls = Fido.fetch(results, path = data_dir, overwrite=False, progress=True)

    gzip_files = glob.glob(data_dir+'*.tar.gz')

    files = tarfile.open(gzip_files[0])
    files.extractall(data_dir)
    files.close()

    os.remove(gzip_files[0])

# Call the routine to enable fitting.
if __name__ == "__main__":
    __spec__ = None

    iris_evts = ['20230503_072923']

    for event in iris_evts:
        print('')
#        print('Downloading data for event '+event)
#        getdata(event)
        fitdata(event)
        print('Event '+event+' processed.')
        print('')


