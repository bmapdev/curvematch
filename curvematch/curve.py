#!/usr/local/epd/bin/python

"""Declaration of the Curve class"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


import numpy as np
from numpy import linalg as LA
from shapeio import curveio


class Curve():

    def __init__(self, coords=np.array([]), attributes=[], dim=0, siz=0, file=None):
        if file is not None:
            self.readcurve(file)
            coords = self.coords

        if coords.size != 0:
            n, T = coords.shape
            if n > T:
                self.coords = np.transpose(coords)
            else:
                self.coords = coords
            self.dim, self.siz = self.coords.shape
        else:
            self.dim = 0
            self.siz = 0
        self.shape = (self.dim, self.siz)

        # self.coords = np.transpose(coords)
        # self.attributes = attributes
        # self.dim = dim
        # self.siz = siz
        # self.shape = (self.dim, self.siz)
        #
        # if coords.size != 0:
        #     self.dim = self.coords.shape[0]
        #     self.siz = self.coords.shape[1]
        #     self.shape = (self.dim, self.siz)

    def winding_number(self):
        pass

    def estimate_pose(self):
        pass

    def repose_curve(self):
        pass

    def length(self):
        arc_length = np.empty((self.dim, self.siz))
        curve_gradient = np.diff(self.coords)
        for i in range(0, self.siz):
            arc_length[:, i] = LA.norm(curve_gradient[:, i])
        return sum(arc_length)

    def compute_curvature(self):
        xprime = np.gradient(self.coords[0,:])
        yprime = np.gradient(self.coords[1,:])
        xdoubleprime = np.gradient(xprime)
        ydoubleprime = np.gradient(yprime)

        if self.dim == 2:
            return (xprime*ydoubleprime - yprime*xdoubleprime)/ \
                        (xprime**2 + yprime**2)**1.5
        elif self.dim == 3:
            zprime = np.gradient(self.coords[2,:])
            zdoubleprime = np.gradient(zprime)
            return np.sqrt((zdoubleprime*yprime - ydoubleprime*zprime)**2 +
                            (xdoubleprime*zprime - zdoubleprime*xprime)**2 +
                            (ydoubleprime*xprime - xdoubleprime*yprime)**2)/ \
                            (xprime**2 + yprime**2 + zprime**2)**1.5

    def smooth_curve(self):
        pass

    def remove_duplicate_vertices(self):
        epsilon = 1/20.0*1.0/self.siz
        curve_gradient = np.zeros((self.dim, self.siz))
        for i in xrange(self.dim):
            curve_gradient[i, :] = np.gradient(self.coords[i, :])

        arc_length = np.zeros(self.siz)
        for i in range(0, self.siz):
            arc_length[i] = LA.norm(curve_gradient[:, i])

        duplicate_status = arc_length > epsilon  # Where False, the index is duplicate
        self.coords = self.coords[:, duplicate_status]
        self.siz = self.coords.shape[1]
        self.shape = self.coords.shape

    def remove_traceover_defects(self):
        pass

    def resample_curve_uniform(self, newsiz=None):
        if newsiz is None:
            newsiz = self.siz

        self.remove_duplicate_vertices()
        curve_gradient = np.zeros((self.dim, self.siz))
        arc_length = np.zeros(self.siz)

        for i in xrange(self.dim):
            curve_gradient[i, :] = np.gradient(self.coords[i, :])

        for i in range(0, self.siz):
            arc_length[i] = LA.norm(curve_gradient[:, i])

        cumulative_arc_length = np.cumsum(arc_length)
        newcoords = np.zeros((self.dim, newsiz))
        for i in range(0, self.dim):
            newcoords[i, :] = np.interp(np.transpose(np.linspace(cumulative_arc_length[0],
                                                      cumulative_arc_length[-1], newsiz)),
                                                      cumulative_arc_length, self.coords[i, :])
        self.coords = newcoords
        self.siz = newsiz

    def readcurve(self,filename):
        self.coords, self.attributes, ismultiUCF = curveio.readcurve(filename)
        self.dim, self.siz = self.coords.shape

    def append_curve(self, other_curve):
        if self.dim != other_curve.dim:
            raise ValueError("Cannot connect curves with mismatched dimensions!")
        for col in other_curve.coords.T:
            if col not in self.coords:
                self.coords = np.append(self.coords, col)
                self.siz += 1
                





