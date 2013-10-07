#!/usr/local/epd/bin/python

"""Declaration of the qshape class"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                   Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from numpy import transpose


class QShape():

    def __init__(self, coords=[], n=0, T=0):
        self.coords = transpose(coords)
        self.n = n
        self.T = T

        if coords is not None:
            self.n = self.coords.shape[0]
            self.T = self.coords.shape[1]
