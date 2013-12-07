#!/usr/local/epd/bin/python
"""Test Curve Matching functions"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch.settings import Settings
from shapeio import curveio
from curvematch.qshape import QShape
import geodesics


def match_curve_pair(curvefilename1, curvefilename2, settings):

    X1, attributes, isMultilevelUCF = curveio.readcurve(curvefilename1)
    X2, attributes, isMultilevelUCF = curveio.readcurve(curvefilename2)

    q1 = QShape(X1)
    q2 = QShape(X2)

    geodesic = geodesics.Geodesic()

    if settings.closed:
        geodesics.compute_for_closed_curves(q1, q2, settings)
    else:
        geodesic = geodesics.compute_for_open_curves_elastic(q1, q2, settings)

    return geodesics
