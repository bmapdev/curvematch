#!/usr/local/epd/bin/python
"""Geodesics for matching qshapes"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


import numpy as np

import utils
from curvematch.settings import Settings
from curvematch.DPmatch import DPmatchcy
from qshape import QShape

class GeodesicsClosed():
    pass


class GeodesicsOpen():
    pass


def compute_flow(q1, tangent_vect, step):
    # TODO uncomment project to tangent
    # TODO uncomment ProjectB


    norm = np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))
    if norm <= 0.001:
        return q1, tangent_vect

    geodesic_path = np.zeros(q1.dim, q1.siz, step)
    qt = q1.copy()
    for i in range(1, step+1):
        qt = qt + norm/step
        #qt = ProjectB(qt)
        geodesic_path.append(qt)

        #tangent_vect = project_tangent(tangent_vect, qt)
        tangent_vect = tangent_vect*norm/np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))

    return qt, geodesic_path


def compute_Palais_inner_prod(tangent_vect1, tangent_vect2):
    steps = len(tangent_vect1)
    inner_prod_val = np.zeros((0, steps))
    for i in range(0, steps):
        inner_prod_val[i] = utils.inner_prod(tangent_vect1[i].coords, tangent_vect2[i].coords)

    return np.trapz(inner_prod_val, np.linspace(0, 1, steps))

def parallel_transport_tangent(tangent_vect, q1, q2):
    # TODO uncomment project to tangent
    norm = np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))

    if norm <= 0.0001:
        parallel_vect = tangent_vect
    else:
        #parallel_vect = project_tangent(tangent_vect, q2)
        #parallel_vect = parallel_vect*norm/np.sqrt(utils.inner_prod(parallel_vect, parallel_vect))
        raise ValueError("Not yet implemented for cases where norm > 0.0001")

    return parallel_vect


def compute_path_derivative(geodesic_path):
# TODO uncomment the project to tangent
    q = geodesic_path[:, :, 0]
    step = geodesic_path.shape[2]
    path_derivative = np.zeros((q.dim, q.siz, step))

    for tau in range(2, step+1):
        path_derivative[:, :, tau] = step*(geodesic_path[:, :, tau] - geodesic_path[:, :, tau-1])
        #path_derivative[:, :, tau] = project_on_tangent_space(path_derivative[:, :, tau])
    return path_derivative


def compute_path_length(path_dt):
    steps = len(path_dt)
    sqrt_inner_prod_val = np.zeros((0, steps))
    for tau in range(0, steps):
        sqrt_inner_prod_val[tau] = utils.inner_prod(path_dt[tau].coords, path_dt[tau].coords)

    return np.trapz(sqrt_inner_prod_val, np.linspace(0, 1, steps))

def compute_ambient(q1, q2, stp):
    geodesic_path = np.zeros((q1.dim, q1.siz, stp))
    for tau in range(0, stp):
        geodesic_path[:, :, tau] = (1-tau*1.0/stp)*q1.coords + tau*1.0/stp*q2.coords

    return geodesic_path


def compute_on_sphere(q1, q2, steps):

    theta = np.arccos(utils.inner_prod(q1, q2))
    f = q2 - utils.inner_prod(q1, q2)*q1
    f = theta*f/np.sqrt(utils.inner_prod(f, f))
    qt = QShape(q1.coords)
    geodesic_sphere = []
    for tau in range(0, steps):
        dt = 1.0*tau/steps
        vnorm = np.sqrt(utils.inner_prod(f, f))
        qt = np.cos(dt*vnorm)*q1 + np.sin(dt*vnorm)*f/vnorm
        qt.project_B()
        geodesic_sphere.append(qt)

    return geodesic_sphere


def compute_for_closed_curves(q1, q2, settings):


    gamma = DPmatchcy.match(q1.coords, q2.coords)
    #q2n = GroupactionGamma(q2,gamma)

    return


def compute_for_open_curves(q1, q2, settings):

    gamma = DPmatchcy.match(q1.coords, q2.coords)
    #q2n = GroupactionGamma(q2,gamma)

    return


def project_tangent(f,q):
    return f - q * utils.inner_prod(f, q)


def compute_cov_derivative_path(alpha):
    n, T, k = np.shape(alpha)
    stp = k - 1
    alpha_dt = np.zeros(np.shape(alpha))
    for tau in xrange(1,k):
        alpha_dt[:, :, tau] = stp*(alpha[:, :, tau] - alpha[:, :, tau-1])
        alpha_dt[:, :, tau] = project_tangent(alpha_dt[:,:, tau],alpha[:, :, tau])
    return alpha_dt
