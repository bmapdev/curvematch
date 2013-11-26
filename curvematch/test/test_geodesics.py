#!/usr/local/epd/bin/python
"""Test Geodesics"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi,  \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch import geodesics
from curvematch.qshape import QShape
from shapeio import curveio
import time
from curvematch.settings import Settings
from sys import stdout
import numpy as np

Settings.output_dir = 'curvematch/test/test/data'
Settings.steps = 7
curve1 = 'curvematch/test/data/q1_partial.ucf'
curve2 = 'curvematch/test/data/q2_partial.ucf'


def test_geodesic_on_sphere():

    X1,attributes,isMultilevelUCF = curveio.readcurve(curve1)
    X2,attributes,isMultilevelUCF = curveio.readcurve(curve2)

    q1 = QShape(X1)
    q2 = QShape(X2)
    t = time.time()

    geodesic = geodesics.compute_on_sphere(q1, q2, Settings.steps)
    for i in geodesic:
        print np.linalg.norm(i.coords)
    elapsed = time.time() - t
    stdout.write('\nElapsed time ' + str(elapsed))
    stdout.flush()


test_geodesic_on_sphere()