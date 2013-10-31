#!/usr/local/epd/bin/python
"""Test Dynamic Programming for matching qshapes"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi,  \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch.qshape import QShape
from shapeio import curveio
import time
from curvematch.DPmatch import DPmatchpy
import numpy as np
from curvematch.settings import Settings
from os import path
from sys import stdout

Settings.output_dir = 'curvematch/test/test/data'
curve1 = 'curvematch/test/data/q1_partial.ucf'
curve2 = 'curvematch/test/data/q2_partial.ucf'


def test_DP_match():

    X1,attributes,isMultilevelUCF = curveio.readcurve(curve1)
    X2,attributes,isMultilevelUCF = curveio.readcurve(curve2)

    q1 = QShape(X1)
    q2 = QShape(X2)
    t = time.time()
    gamma, Energy = DPmatchpy.matchq(q1,q2)
    np.savetxt(path.join(Settings.output_dir, 'gamma.txt'), gamma)
    np.savetxt(path.join(Settings.output_dir, 'Energy.txt'), Energy)

    elapsed = time.time() - t
    stdout.write('\nElapsed time ' + str(elapsed))
    stdout.flush()
