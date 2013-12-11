#!/usr/local/epd/bin/python

"""Declaration of the Utility module"""

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import numpy as np
from math import pi
from numpy import linalg as LA


def inner_prod(coords1, coords2):
    # Assume coords1 and coords2 are numpy arrays
    value = np.trapz(sum(coords1*coords2), np.linspace(0, 2*pi, coords1.shape[1]))
    return value


def shift_vertices(p, tau):
    shape_p = np.shape(p)
    if tau == 0:
        pn = p
    elif tau > 0:
        t = abs(tau)
        pn = np.zeros(shape_p, np.int32)
        pn[:, :-t] = p[:, t:]
        pn[:, -t:] = p[:, :t]
    elif tau < 0:
        t = abs(tau)
        pn = np.zeros(shape_p,np.int32)
        pn[:, :t] = p[:, -t:]
        pn[:, t:] = p[:, :-t]
    return pn


def find_best_rotation(q1, q2):
    #assumes starting points are fixed
    A = np.dot(q1.coords, q2.coords.transpose())
    # A = q1.coords*q2.coords.transpose()
    U, S, V = LA.svd(A)
    if np.absolute(LA.det(U)*LA.det(V) -1) < 10*np.spacing(1):
        S = np.eye(q1.dim)
    else:
        S = np.eye(q1.dim)
        S[:, -1] = -S[:, -1]
    R = np.dot(np.dot(U, S), np.transpose(V))  # R=U*S*V' (matrix multiplication)
    q2new = q2.copy()
    q2new.coords = np.dot(R, q2.coords)
    return q2new, R


