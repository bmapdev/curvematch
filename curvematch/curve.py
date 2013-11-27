#!/usr/local/epd/bin/python

"""Declaration of the Curve class"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


import numpy as np

from shapeio import curveio
from numpy import transpose, array, empty, diff, cumsum, interp, linspace, logical_or, gradient, sqrt
from numpy.linalg import norm


class Curve():

    def __init__(self, coords=np.array([]), attributes=[], dim=0, siz=0):
        self.coords = transpose(coords)
        self.attributes = attributes
        self.dim = dim
        self.siz = siz
        self.shape = (self.dim, self.siz)

        if coords is not None:
            self.dim = self.coords.shape[0]
            self.siz = self.coords.shape[1]
            self.shape = (self.dim, self.siz)
    def winding_number(self):
        pass

    def estimate_pose(self):
        pass

    def repose_curve(self):
        pass

    def length(self):
        arc_length = empty((self.dim, self.siz))
        curve_gradient = diff(self.coords)
        for i in range(0, self.siz):
            arc_length[:, i] = norm(curve_gradient[:, i])
        return sum(arc_length)

    def compute_curvature(self):
        xprime = gradient(self.coords[0,:])
        yprime = gradient(self.coords[1,:])
        xdoubleprime = gradient(xprime)
        ydoubleprime = gradient(yprime)

        if self.dim == 2:
            return (xprime*ydoubleprime - yprime*xdoubleprime)/ \
                        (xprime**2 + yprime**2)**1.5
        elif self.dim == 3:
            zprime = gradient(self.coords[2,:])
            zdoubleprime = gradient(zprime)
            return sqrt((zdoubleprime*yprime - ydoubleprime*zprime)**2 +
                            (xdoubleprime*zprime - zdoubleprime*xprime)**2 +
                            (ydoubleprime*xprime - xdoubleprime*yprime)**2)/ \
                            (xprime**2 + yprime**2 + zprime**2)**1.5

    def smooth_curve(self):
        pass

    def remove_duplicate_vertices(self):
        epsilon = 1/20*1/self.siz
        curve_gradient = diff(self.coords)
        duplicate_status = curve_gradient > epsilon
        bitwise_or_flag = array(self.siz*[True], bool)
        # Logical OR all the components. Where False, the index is duplicate
        if self.dim == 2:
            bitwise_or_flag = logical_or(duplicate_status[0, :], duplicate_status[1, :])
        elif self.dim == 3:
            bitwise_or_flag = logical_or(duplicate_status[0, :], duplicate_status[1, :])
            bitwise_or_flag = logical_or(bitwise_or_flag, duplicate_status[2, :])
        self.coords = self.coords[:, bitwise_or_flag]

    def remove_traceover_defects(self):
        pass

    def resample_curve_uniform(self):
        self.remove_duplicate_vertices()
        #curve_gradient = empty((self.dim, self.siz))
        arc_length = empty((self.dim, self.siz))

        curve_gradient = diff(self.coords)
        #for i in range(0, self.dim):
        #    curve_gradient[i, :] = diff(self.coords[i, :])
        for i in range(0, self.siz):
            arc_length[:, i] = norm(curve_gradient[:, i])

        cumulative_arc_length = cumsum(arc_length)

        for i in range(0, self.dim):
            self.coords[i, :] = interp(cumulative_arc_length[i, :], curve_gradient,
                                       linspace(cumulative_arc_length[0], cumulative_arc_length[-1], self.siz))

    def readcurve(self,filename):
        self.coords, self.attributes, ismultiUCF = curveio.readcurve(filename)
        self.dim, self.siz = self.coords.shape

