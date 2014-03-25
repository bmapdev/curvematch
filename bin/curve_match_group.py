#!/usr/local/epd/bin/python
"""Curve Match"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import argparse
from curvematch import match
import numpy as np
from curvematch.settings import Settings
import os
from curvematch import plotting
from curvematch import geodesics


def main():
    parser = argparse.ArgumentParser(description='Elastic Shape matching of a pair of curves.\n')
    parser.add_argument('-list', help='source curve list')
    parser.add_argument('-target', help='target curve filename')
    parser.add_argument('-open', dest='open', help='match open curves', action='store_true', default='False', required=False)
    parser.add_argument('-linear', dest='linear', help='use linear matching', action='store_true', default=False, required=False)
    parser.add_argument('-norotate', dest='norotate', help='do not align rotations', action='store_false', default=True, required=False)
    parser.add_argument('-noplot', dest='noplot', help='do not save graphical plots', action='store_true', default=False, required=False)
    parser.add_argument('-odirlist', dest='odirlist', help='output directory list', required=True)
    parser.add_argument('-resize', dest='resize', help='resize to the specified number of vertices', required=False, default=200)
    args = parser.parse_args()
    curve_match_group(args.list, args.target, args.odirlist, args.open, args.linear, args.norotate, args.resize, args.noplot)


def curve_match_group(srclist, target, odirlist, openflag=False, linearflag=False, norotateflag=False, resize=200, noplotflag=False):
    settings = Settings()

    if openflag:
        settings.closed = False

    src_curve_paths_file = open(srclist, 'r')
    src_curve_paths = src_curve_paths_file.read().split('\n')

    output_dirs_file = open(odirlist, 'r')
    output_dirs = output_dirs_file.read().split('\n')

    if len(output_dirs) != len(src_curve_paths):
        os.sys.stdout.write('The number of output directories should be same as number of input curves. Now exiting.\n')
        os.sys.exit()

    geodesic_array, src_curve_matched_to_target_array, target_curve_array = \
        match.match_curve_group(src_curve_paths, target, openflag=openflag, linearflag=linearflag, norotateflag=norotateflag, resize=resize)

    for idx, geodesic in enumerate(geodesic_array):
        curve_path = geodesics.to_curve_path(geodesic.path)

        if not os.path.exists(output_dirs[idx]):
            os.mkdir(output_dirs[idx])

        if not noplotflag:
            plotting.plot_path(curve_path, filename=os.path.join(output_dirs[idx], 'geodesic_path.pdf'))

        _, input_src_curve_name = os.path.split(src_curve_paths[idx])
        _, input_target_curve_name = os.path.split(target)

        src_curve_matched_to_target_array[idx].writecurve(os.path.join(output_dirs[idx], os.path.splitext(input_src_curve_name)[0] + '_matched.ucf'))
        target_curve_array[idx].writecurve(os.path.join(output_dirs[idx], input_target_curve_name))
        plotting.plot_matching(os.path.join(output_dirs[idx], "geodesic_alignment"), target_curve_array[idx], src_curve_matched_to_target_array[idx])
        np.savetxt(os.path.join(output_dirs[idx], 'geodesic_shape_distance.txt'), [geodesic.geodesic_distance], fmt='%.8f')
        np.savetxt(os.path.join(output_dirs[idx], 'reparameterization.txt'), geodesic.gamma, fmt='%.8f')
        np.savetxt(os.path.join(output_dirs[idx], 'tangent_vector.txt'), np.transpose(geodesic.tangent_vect[1].coords), fmt='%.8f')

if __name__ == '__main__':
    main()
