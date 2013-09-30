#!/usr/local/epd/bin/python

"""Test Dynamic Programming for matching qshapes"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from qshape import qshape
from shapeio import curveio
import time
import DPmatch
import numpy as np

curve1 = 'q1_partial.ucf'
curve2 = 'q2_partial.ucf'

X1,attributes,isMultilevelUCF = curveio.readcurve(curve1)
X2,attributes,isMultilevelUCF = curveio.readcurve(curve2)

q1 = qshape(X1)
q2 = qshape(X2)
t = time.time()
gamma, Energy = DPmatch.matchq(q1,q2)
np.savetxt('gamma_new.txt',gamma)
np.savetxt('Energy_new.txt',Energy)

elapsed = time.time() - t
print elapsed
