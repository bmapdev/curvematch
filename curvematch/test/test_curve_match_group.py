#!/usr/local/epd/bin/python
"""Test Curve Group Matching function"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch import match
from curvematch import plotting
from curvematch import geodesics
import numpy as np
import os

srclist = 'curvematch/test/data/group/curvelist.txt'
target = 'curvematch/test/data/group/curve1_cc_JOINED.ucf'
odirlist = 'curvematch/test/data/group/odirlist.txt'

def test_curvematch_group():

    src_curve_paths_file = open(srclist, 'r')
    src_curve_paths = src_curve_paths_file.read().split('\n')

    output_dirs_file = open(odirlist, 'r')
    output_dirs = output_dirs_file.read().split('\n')

    geodesic_array, src_curve_matched_to_target_array, target_curve_array = \
        match.match_curve_group(src_curve_paths, target, openflag=True, linearflag=False)

    for idx, geodesic in enumerate(geodesic_array):
        curve_path = geodesics.to_curve_path(geodesic.path)

        if not os.path.exists(output_dirs[idx]):
            os.mkdir(output_dirs[idx])

        plotting.plot_path(curve_path, filename=os.path.join(output_dirs[idx], 'geodesic_path.pdf'))

        _, input_src_curve_name = os.path.split(src_curve_paths[idx])
        _, input_target_curve_name = os.path.split(target)

        src_curve_matched_to_target_array[idx].writecurve(os.path.join(output_dirs[idx], os.path.splitext(input_src_curve_name)[0] + '_matched.ucf'))
        target_curve_array[idx].writecurve(os.path.join(output_dirs[idx], input_target_curve_name))
        plotting.plot_matching(os.path.join(output_dirs[idx], "geodesic_alignment"), target_curve_array[idx], src_curve_matched_to_target_array[idx])
        plotting.scalar_function_plot(os.path.join(output_dirs[idx], "gamma"), geodesic.gamma)
        np.savetxt(os.path.join(output_dirs[idx], 'geodesic_shape_distance.txt'), [geodesic.geodesic_distance], fmt='%.8f')
        np.savetxt(os.path.join(output_dirs[idx], 'reparameterization.txt'), geodesic.gamma, fmt='%.8f')
        np.savetxt(os.path.join(output_dirs[idx], 'tangent_vector.txt'), np.transpose(geodesic.tangent_vect[1].coords), fmt='%.8f')

test_curvematch_group()
