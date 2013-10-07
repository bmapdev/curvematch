#!/usr/local/epd/bin/python
"""Nose testing for Plotting"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch import plotting
from curvematch.settings import Settings
from shapeio import curveio
from curvematch.qshape import QShape

Settings.output_dir = 'test/data'
curve1 = 'test/data/curve1.ucf'


def test_plotcurve_with_name():
    coords,atts,ismultiUCF = curveio.readcurve(curve1)
    q1 = QShape(coords)
    plotting.plotcurve(q1, filename='plot_test')
    pass


def test_plotcurve_without_name():
    coords,atts,ismultiUCF = curveio.readcurve(curve1)
    q1 = QShape(coords)
    plotting.plotcurve(q1)
    pass
