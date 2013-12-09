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


Settings.output_dir = 'curvematch/test/data'
curve1 = 'curvematch/test/data/curve1.ucf'
curve2 = 'curvematch/test/data/curve2.ucf'
settings = Settings()
settings.closed = False

def test_curvematch():
    geodesic = match.match_curve_pair(curve1, curve2, settings)
    curve_path = geodesics.to_curve_path(geodesic.path)
    plotting.plot_path(curve_path)

test_curvematch()