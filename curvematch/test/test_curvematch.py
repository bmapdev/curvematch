#!/usr/local/epd/bin/python
"""Test Curve Matching functions"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

# import curvematch
from curvematch.settings import Settings
# from shapeio import curveio
# from curvematch.qshape import QShape
from curvematch import match
from curvematch import plotting
from curvematch import geodesics
import time
from sys import stdout

Settings.output_dir = 'curvematch/test/data'
curve1 = 'curvematch/test/data/04_Curve16.ucf'
curve2 = 'curvematch/test/data/05_Curve16.ucf'
# curve1 = 'curvematch/test/data/curve1.ucf'
# curve2 = 'curvematch/test/data/curve2.ucf'

settings = Settings()
settings.closed = False


def test_curvematch():
    t = time.time()
    geodesic = match.match_curve_pair(curve1, curve2, settings, rotation=False)
    print geodesic.geodesic_distance
    curve_path = geodesics.to_curve_path(geodesic.path)
    plotting.plot_path(curve_path, 'elastic_path')
    elapsed = time.time() - t
    stdout.write('\nElapsed time ' + str(elapsed))
    stdout.flush()

test_curvematch()