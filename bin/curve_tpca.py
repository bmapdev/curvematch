"""Curve TPCA"""

__author__ = "Suchit Panjiyar"
__copyright__ = "Copyright 2016, Shantanu H. Joshi, Suchit Panjiyar, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "suchit@ucla.edu"


import argparse
import numpy as np
from curvematch.settings import Settings
import os
from curvematch.curve import Curve
from curvematch.qshape import QShape
from curve_mean import curve_mean
from curvematch import curvestats
from scipy import io

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
    qmean, tangent_vectors = curve_mean(args.list, args.odir, args.open, args.linear, args.norotate, args.resize, args.noplot)
    a = curvestats.pca(qmean, tangent_vectors, args.odir)
    io.savemat(os.path.join(args.odir, 'tpca'), a)


if __name__ == '__main__':
    main()


