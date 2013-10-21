#!/usr/local/epd/bin/python

"""Declaration of the qshape class"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from numpy import transpose


class QShape():

    def __init__(self, coords=[], dim=0, siz=0):
        self.coords = transpose(coords)
        self.dim = dim
        self.siz = siz

        if coords is not None:
            self.dim = self.coords.shape[0]
            self.siz = self.coords.shape[1]

    def fromcurve(self, curve):
        pass

