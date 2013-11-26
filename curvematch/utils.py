#!/usr/local/epd/bin/python

"""Declaration of the Utility module"""

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import numpy as np
from math import pi
from numpy import linalg as LA


def inner_prod(q1, q2):
    value = np.trapz(sum(q1.coords*q2.coords), np.linspace(0, 2*pi, q1.siz))
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


def Find_Best_Rotation(q1, q2):
    #assumes starting points are fixed
    A = q1*q2.transpose()
    U, S, V = LA.svd(A)
    if np.absolute(np.det(U*np.det(V) -1)) < 10*np.spacing(1) :
        S = np.eye(q1.dim)
    else:
        S = np.eye(q1.dim)
        S[:, -1] = -S[:, -1]
    R = U*S*V
    q2new = R*q2
    return q2new


