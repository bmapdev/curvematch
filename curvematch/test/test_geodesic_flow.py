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
import pkg_resources
from scipy.io import savemat

Settings.output_dir = pkg_resources.resource_filename('curvematch', 'test/test/data')
Settings.steps = 7
curve1 = pkg_resources.resource_filename('curvematch', 'test/data/q1_partial.ucf')
curve2 = pkg_resources.resource_filename('curvematch', 'test/data/q2_partial.ucf')

q1_array = np.loadtxt('~/gdrive/research/projects/shapes/curves/elastic-path-straightening/open/matlab/ver4.0DP+gradientMEX/q1.txt')
w_array = np.loadtxt('~/gdrive/research/projects/shapes/curves/elastic-path-straightening/open/matlab/ver4.0DP+gradientMEX/w.txt')

def test_geodesic_flow():
    q1 = QShape(coords=q1_array)
    w, qt, geodesic_path = geodesics.compute_flow_open(q1, w_array, Settings.steps)
    savemat('~/Desktop/data.mat', {'qt': qt, 'geodesic_path': geodesic_path, 'w': w})


test_geodesic_flow()