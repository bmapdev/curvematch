#!/usr/local/epd/bin/python
"""Test Dynamic Programming for matching qshapes"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi,  \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from shapeio import curveio
from curvematch.qshape import QShape
import numpy as np

import time
from curvematch.DPmatch import DPmatchcy
from curvematch.settings import Settings
from sys import stdout
import pkg_resources

Settings.output_dir = pkg_resources.resource_filename('curvematch', 'test/data')
curve1 = pkg_resources.resource_filename('curvematch', 'test/data/curve1.ucf')
curve2 = pkg_resources.resource_filename('curvematch','test/data/curve2.ucf')


def test_DP_matchcy():

    X1, attributes, isMultilevelUCF = curveio.readcurve(curve1)
    X2, attributes, isMultilevelUCF = curveio.readcurve(curve2)

    q1 = QShape(X1)
    q2 = QShape(X2)
    t = time.time()
    gamma = DPmatchcy.match(q1.coords, q2.coords)
    gamma = 2 * np.pi * gamma/np.max(gamma)
    print gamma

    elapsed = time.time() - t
    stdout.write('\nElapsed time ' + str(elapsed))
    stdout.flush()
    return gamma

