import glob
import os

from fit_iris_lines import fit_raster
import iris_get_mg_features_lv2 as get_mg
import get_quartiles
import asdf
from astropy.io import fits
import datetime as dt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

event = '20230503_070923'
output_loc = r"C:\Users\molly\OneDrive\OneDrive - Dublin City University personal\PHA4\Final_Year_Project\outputs"

# Plot the output of the fitting routine if fitting sit-and-stare data
def plot_iris_sns_quartile_fits(int_map, dopp_map, width_map, asym_map, iris_window, event, main_header):
    fig = plt.figure(constrained_layout=True, figsize=(10, 10))
    plt.rcParams['font.size'] = '10'

    gs = GridSpec(nrows=4, ncols=1, hspace=0.05, wspace=0.05)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.95)

    plot_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'),
                                     '%Y/%m/%dT%H:%M:%S')
    file_time = dt.datetime.strftime(dt.datetime.strptime(int_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f'),
                                     '%Y%m%d_%H%M%S')

    dopp_rng = 10
    max_wid = 0.2
    asym_rng = 1

    # Set the plotting parameters
    cadence = main_header['STEPT_AV']
    t_array = np.arange(0, int_map.data.shape[1]) * cadence
    slit_pos = int_map.meta['crval2'] + int_map.meta['cdelt2'] * (
                np.arange(int_map.data.shape[0]) - int_map.meta['crpix2'])

    # Intensity map
    ax1 = fig.add_subplot(gs[0, 0], label='a)')
    alpha = 1
    upr_bnd = np.nanpercentile(int_map.data, 100 - alpha)

    norm = colors.Normalize(vmin=0, vmax=upr_bnd)

    plt.imshow(int_map.data, norm=norm, cmap=mpl.colormaps['Reds_r'], axes=ax1, origin='lower',
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax1.set_ylabel(" ")
    ax1.set_xlabel(" ")
    ax1.set_xticklabels([])

    plt.colorbar(location='right', label=r'a) Integrated Intensity', shrink=0.6, ax=ax1, ticks=[0, upr_bnd])

    # Asymmetry map
    ax2 = fig.add_subplot(gs[1, 0], label='b)')
    norm = colors.Normalize(vmin=-asym_rng, vmax=asym_rng)
    plt.imshow(asym_map.data, norm=norm, cmap=mpl.colormaps['seismic'], axes=ax2, origin='lower',
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax2.set_ylabel(" ")
    ax2.set_xlabel(" ")
    ax2.set_xticklabels([])

    plt.colorbar(location='right', label=r'b) RB Asym.', shrink=0.6, ax=ax2)

    # Doppler map
    ax3 = fig.add_subplot(gs[2, 0], label='c)')
    norm = colors.Normalize(vmin=-dopp_rng, vmax=dopp_rng)
    plt.imshow(dopp_map.data, norm=norm, cmap=mpl.colormaps['coolwarm'], axes=ax3, origin='lower',
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax3.set_ylabel("Solar Y (arcsec)")
    ax3.set_xlabel(" ")
    ax3.set_xticklabels([])

    plt.colorbar(location='right', label=r'c) v$_{dopp}$ ($km~s^{-1}$)', shrink=0.6, ax=ax3)

    # Line width
    ax4 = fig.add_subplot(gs[3, 0], label='d)')
    norm = colors.Normalize(vmin=0, vmax=max_wid)
    plt.imshow(width_map.data, norm=norm, cmap=mpl.colormaps['cubehelix'], axes=ax4, origin='lower',
               extent=[t_array.min(), t_array.max(), slit_pos.min(), slit_pos.max()], aspect='auto')
    ax4.set_ylabel(" ")
    ax4.set_xlabel("Time from raster start (s)")

    plt.colorbar(location='right', label=r'd) Width ($\AA$)', shrink=0.6, ax=ax4)

    plt.suptitle(iris_window + r'$\AA$; ' + plot_time)
    plt.savefig(
        output_loc + event + 'IRIS_analysis_quartile_' + iris_window.replace(' ', '_') + '_' + file_time + '.png',
        bbox_inches='tight')
    plt.close(fig)

asdf_files = glob.glob(os.path.join(output_loc, 'IRIS_fitting_C_II_1335_20230503_072923.asdf'))

if not asdf_files:
    raise FileNotFoundError(f"No ASDF files found in {output_loc}")



# Get map names
#----------------------
for asdf_file in asdf_files:
    iris_window_underscore = os.path.basename(asdf_file).replace('IRIS_fitting_', '').replace(event, '').replace('.asdf', '')
    iris_window = iris_window_underscore.replace('_', ' ').strip()


    print(iris_window)
    with asdf.open(asdf_file) as af:
        possible_int_keys = [key for key in af.tree.keys() if 'int' in key.lower()]
        if not possible_int_keys:
            raise KeyError(f"No intensity map found in {asdf_file}")
        int_key = possible_int_keys[0]
        int_map = af.tree[int_key]

        possible_dopp_keys = [key for key in af.tree.keys() if 'dopp' in key.lower()]
        if not possible_dopp_keys:
            raise KeyError(f"No Doppler map found in {asdf_file}")
        dopp_key = possible_dopp_keys[0]
        dopp_map = af.tree[dopp_key]

        possible_width_keys = [key for key in af.tree.keys() if 'width' in key.lower()]
        if not possible_width_keys:
            raise KeyError(f"No width map found in {asdf_file}")
        width_key = possible_width_keys[0]
        width_map = af.tree[width_key]

        possible_vnt_keys = [key for key in af.tree.keys() if 'vnt' in key.lower()]
        if possible_vnt_keys:
            vnt_key = possible_vnt_keys[0]
            vnt_map = af.tree[vnt_key]
        else:
            vnt_map = None
            print(f'Warning: No VNT map found in {asdf_file}. Skipping vnt plot.')


        possible_asym_keys = [key for key in af.tree.keys() if 'asym' in key.lower()]
        if not possible_asym_keys:
            raise KeyError(f"No asym map found in {asdf_file}")
        asym_key = possible_asym_keys[0]
        asym_map = af.tree[asym_key]

    # Find ranges

    # Define per-window plotting ranges
    plot_ranges = {
        "C II 1334": {"dopp_rng": 10, "max_wid": 0.1, "asym_rng": 1, "max_vnt": 30},
        "C II 1335": {"dopp_rng": 10, "max_wid": 0.1, "asym_rng": 1, "max_vnt": 30},
    }

    # Default if a window is not in the dictionary
    default_ranges = {"dopp_rng": 10, "max_wid": 0.1, "asym_rng": 1, "max_vnt": 30}

    # Get cadence from .fits file as asdf does not contain it
    iris_fits = glob.glob(os.path.join(output_loc, "*.fits"))
    main_header = fits.getheader(iris_fits[0], 0)

    # Ensure output directory exists
    ranges = plot_ranges.get(iris_window, default_ranges)
    dopp_rng = ranges["dopp_rng"]
    max_wid = ranges["max_wid"]
    asym_rng = ranges["asym_rng"]
    max_vnt = ranges["max_vnt"]

#def plot_iris_sns_quartile_fits(int_map, dopp_map, width_map, asym_map, iris_window, event, main_header):
plot_iris_sns_quartile_fits(int_map, dopp_map, width_map, asym_map, 'C II 1335', event, main_header)