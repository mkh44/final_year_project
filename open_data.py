##!/usr/bin/env python3
# -*- coding: utf-8 -*-

import asdf
import os

if os.uname().sysname == 'Linux':
    IRIS_data_loc = '/mnt/nas/ug/hurlem24/iris_data/iris_input_data/'
    output_loc = '/mnt/nas/ug/hurlem24/iris_data/iris_output/iris_output/'
else:
    IRIS_data_loc = '/mnt/nas/ug/hurlem24/iris_data/iris_input_data/'
    output_loc = '/mnt/nas/ug/hurlem24/iris_data/iris_output/iris_output/'

def open_iris(event, file):

    file_path = os.path.join(output_loc, event, file)
    
    with asdf.open(file) as af:
        int_map = af.tree['int_map']
        dopp_map = af.tree['dopp_map']
        width_map = af.tree['width_map']
        vnt_map = af.tree['vnt_map']
        asym_map = af.tree['asym_map']

    return int_map, dopp_map, width_map, vnt_map, asym_map


iris_evts = ['20230503_072923']
file = 'IRIS_fitting_C_II_1334_20230503_072923.asdf'

for event in iris_evts:
    open_iris(event, file)

