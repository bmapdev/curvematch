#!/usr/local/epd/bin/python
"""Test the Curve class"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch.settings import Settings
from curvematch.qshape import QShape
from shapeio import curveio

Settings.output_dir = 'curvematch/test/data'
curve1 = 'curvematch/test/data/curve1.ucf'


