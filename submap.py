import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import asdf

#find file and read intensity

top_right = SkyCoord(#coords , frame=swap_map.coordinate_frame)
bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame=swap_map.coordinate_frame)
swap_submap = swap_map.submap(bottom_left, top_right=top_right)
