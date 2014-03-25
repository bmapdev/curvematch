#!/usr/local/epd/bin/python
"""Curve Match"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import argparse
import os
from curvematch.settings import Settings
from curvematch import match
from curvematch import plotting
from curvematch import geodesics
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='Elastic Shape matching of a pair of curves.\n')
    parser.add_argument('src', help='source curve')
    parser.add_argument('target', help='target curve')
    parser.add_argument('-open', dest='open', help='match open curves', action='store_true', default='False', required=False)
    parser.add_argument('-linear', dest='linear', help='use linear matching', action='store_true', default=False, required=False)
    parser.add_argument('-norotate', dest='norotate', help='do not align rotations', action='store_false', default=True, required=False)
    parser.add_argument('-noplot', dest='noplot', help='do not save graphical plots', action='store_true', default=False, required=False)
    parser.add_argument('-odir', dest='odir', help='output directory', required=True)
    parser.add_argument('-resize', dest='resize', help='resize to the specified number of vertices', required=False, default=100)
    args = parser.parse_args()
    curve_match(args.src, args.target, args.open, args.linear, args.norotate, args.resize, args.noplot, args.odir)


def curve_match(src, target, openflag, linearflag, norotateflag, resize, noplotflag, odir):
    settings = Settings()

    if openflag:
        settings.closed = False

    geodesic, target_curve, src_curve_matched_to_target = match.match_curve_pair(target, src, settings, rotation=norotateflag, siz=resize, return_curves=True, linear=linearflag)
    curve_path = geodesics.to_curve_path(geodesic.path)

    if not os.path.exists(odir):
        os.mkdir(odir)

    if not noplotflag:
        plotting.plot_path(curve_path, filename=os.path.join(odir, 'geodesic_path.pdf'))

    _, input_src_curve_name = os.path.split(src)
    _, input_target_curve_name = os.path.split(target)

    src_curve_matched_to_target.writecurve(os.path.join(odir, os.path.splitext(input_src_curve_name)[0] + '_matched.ucf'))
    target_curve.writecurve(os.path.join(odir, input_target_curve_name))
    np.savetxt(os.path.join(odir, 'geodesic_shape_distance.txt'), [geodesic.geodesic_distance], fmt='%.8f')
    np.savetxt(os.path.join(odir, 'reparameterization.txt'), geodesic.gamma, fmt='%.8f')
    np.savetxt(os.path.join(odir, 'tangent_vector.txt'), np.transpose(geodesic.tangent_vect[1].coords), fmt='%.8f')

if __name__ == '__main__':
    main()
