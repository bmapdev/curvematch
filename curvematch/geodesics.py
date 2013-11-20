#!/usr/local/epd/bin/python
"""Geodesics for matching qshapes"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from curvematch.settings import Settings
from curvematch.DPmatch import DPmatchcy
import numpy as np
import utils

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
    step = tangent_vect1.shape[2]
    inner_prod_val = np.zeros((0, step))
    for i in range(0, step):
        inner_prod_val[i] = utils.inner_prod(tangent_vect1[:, :, i], tangent_vect2[:, :, i])

    return np.trapz(inner_prod_val, np.linspace(0, 1, step))

def parallel_transport_tangent(tangent_vect, q1, q2):
    # TODO uncomment project to tangent
    norm = np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))

    if norm <= 0.0001:
        parallel_vect = tangent_vect
    else:
        #parallel_vect = project_tangent(tangent_vect, q2)
        #parallel_vect = parallel_vect*norm/np.sqrt(utils.inner_prod(parallel_vect, parallel_vect))
        pass

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


def compute_ambient(q1, q2, stp):
    geodesic_path = np.zeros((q1.dim, q1.siz, stp))
    for tau in range(0, stp):
        geodesic_path[:, :, tau] = (1-tau*1.0/stp)*q1.coords + tau*1.0/stp*q2.coords

    return geodesic_path


def compute_on_sphere():
    pass


def compute_for_closed_curves(q1, q2, settings):


    gamma = DPmatchcy.match(q1.coords, q2.coords)
    #q2n = GroupactionGamma(q2,gamma)

    return


def compute_for_open_curves(q1, q2, settings):

    gamma = DPmatchcy.match(q1.coords, q2.coords)
    #q2n = GroupactionGamma(q2,gamma)

    return
