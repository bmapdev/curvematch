"""Curve Mean"""

__author__ = "Suchit Panjiyar"
__copyright__ = "Copyright 2016, Shantanu H. Joshi, Suchit Panjiyar, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "suchit@ucla.edu"


import argparse
import numpy as np
from curvematch.settings import Settings
import os
import scipy            #used for saving .mat files
from scipy import io    #not currently in use
from curvematch.curve import Curve
from curvematch.qshape import QShape
from curvematch import mean


def main():
    parser = argparse.ArgumentParser(description='Elastic Shape matching of a pair of curves.\n')
    parser.add_argument('-list', help='source curve list')
    parser.add_argument('-odir', dest='odir', help='output directory', required=True)
    parser.add_argument('-open', dest='open', help='match open curves', action='store_true', default=False, required=False)
    parser.add_argument('-linear', dest='linear', help='use linear matching', action='store_true', default=False, required=False)
    parser.add_argument('-norotate', dest='norotate', help='do not align rotations', action='store_false', default=True, required=False)
    parser.add_argument('-noplot', dest='noplot', help='do not save graphical plots', action='store_true', default=False, required=False)
    parser.add_argument('-resize', dest='resize', help='resize to the specified number of vertices', required=False, default=200)
    args = parser.parse_args()
    curve_mean(args.list, args.odir, args.open, args.linear, args.norotate, args.resize, args.noplot)


def curve_mean(srclist, odir, openflag=False, linearflag=False, norotateflag=False, resize=200, noplotflag=False):
    settings = Settings()

    settings.closed = False

    src_curve_paths_file = open(srclist, 'r')
    src_curve_paths = src_curve_paths_file.read().split('\n')

    if not os.path.exists(odir):
        os.mkdir(odir)

    qarray = []

    for i in src_curve_paths:
        c_i = Curve(file=i)
        c_i.resample_curve_uniform(resize)
        q_i = QShape()
        q_i.from_curve(c_i)
        qarray.append(q_i)

    qmean,tangent_vectors = mean.karcher_mean(qarray, settings)

    curve = qmean.to_curve()

    curve.writecurve(os.path.join(odir, "mean.ucf"))
    np.savetxt(os.path.join(odir, 'mean.txt'), qmean.coords)
    for i in range(0,len(tangent_vectors)):
        temp_str = 'tangent_vector_'+ str(i) +'.txt'
        np.savetxt(os.path.join(odir, temp_str), tangent_vectors[i].coords)

if __name__ == '__main__':
    main()





