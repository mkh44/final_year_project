import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import IRIS_fitting_C_II_1334_20230503_072923.asdf as input_file

swap_map = sunpy.map.Map(input_file)

top_right = SkyCoord(0 * u.arcsec, -200 * u.arcsec, frame=swap_map.coordinate_frame)
bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame=swap_map.coordinate_frame)
swap_submap = swap_map.submap(bottom_left, top_right=top_right)
