#!/usr/local/epd/bin/python
"""Geodesics for matching qshapes"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2012, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


import numpy as np
import numpy.linalg as LA
from math import pi
import utils
from curvematch.settings import Settings
from curvematch.DPmatch import DPmatchcy
from qshape import QShape


class Geodesic(object):
    path = []
    gamma = []
    tangent_vect = []
    geodesic_distance = 0


class GeodesicsClosed(Geodesic):
    pass


class GeodesicsOpen(Geodesic):
    pass


def compute_flow(q1, tangent_vect, step):
    # TODO uncomment project to tangent
    # TODO uncomment ProjectB


    norm = np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))
    if norm <= 0.001:
        return q1, tangent_vect

    geodesic_path = np.zeros(q1.dim(), q1.siz(), step)
    qt = q1.copy()
    for i in range(1, step+1):
        qt = qt + norm/step
        #qt = ProjectB(qt)
        geodesic_path.append(qt)

        #tangent_vect = project_tangent(tangent_vect, qt)
        tangent_vect = tangent_vect*norm/np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))

    return qt, geodesic_path


def compute_flow_open(q1, tangent_vect, step):
    # TODO uncomment project to tangent
    # TODO uncomment ProjectB

    norm = np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))
    if norm <= 0.001:
        return q1, tangent_vect

    geodesic_path = []
    qt = q1.copy()
    geodesic_path.append(q1)
    for i in range(1, step+1):
        qt.coords += tangent_vect/step
        qt.project_b()
        geodesic_path.append(qt.copy())

        tangent_vect = utils.project_tangent_q(tangent_vect, qt)
        tangent_vect = tangent_vect*norm/np.sqrt(utils.inner_prod(tangent_vect, tangent_vect))

    return tangent_vect, qt, geodesic_path


def compute_Palais_inner_prod(tangent_vect1, tangent_vect2):
    steps = len(tangent_vect1)
    inner_prod_val = np.zeros(steps)
    for i in range(0, steps):
        inner_prod_val[i] = utils.inner_prod(tangent_vect1[i].coords, tangent_vect2[i].coords)

    return np.trapz(inner_prod_val, np.linspace(0, 2*pi, steps))

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
    path_derivative = np.zeros((q.dim(), q.siz(), step))

    for tau in range(2, step+1):
        path_derivative[:, :, tau] = step*(geodesic_path[:, :, tau] - geodesic_path[:, :, tau-1])
        # path_derivative[:, :, tau] = project_tangent(path_derivative[:, :, tau])
    return path_derivative


def compute_path_length(path_dt):
    steps = len(path_dt)
    sqrt_inner_prod_val = np.zeros(steps)
    for tau in range(0, steps):
        sqrt_inner_prod_val[tau] = np.sqrt(utils.inner_prod(path_dt[tau].coords, path_dt[tau].coords))

    return np.trapz(sqrt_inner_prod_val, np.linspace(0, 2*pi, steps))


def compute_ambient(q1, q2, stp):
    geodesic_path = np.zeros((q1.dim(), q1.siz(), stp))
    for tau in range(0, stp):
        geodesic_path[:, :, tau] = (1-tau*1.0/stp)*q1.coords + tau*1.0/stp*q2.coords

    return geodesic_path


def compute_on_sphere(q1, q2, steps):

    theta = np.arccos(utils.inner_prod(q1.coords, q2.coords))
    f = q2 - utils.inner_prod(q1.coords, q2.coords)*q1
    f = theta*f/np.sqrt(utils.inner_prod(f.coords, f.coords))
    qt = QShape(q1.coords)
    geodesic_sphere = []
    for tau in xrange(0, steps+1):
        dt = 1.0*tau/steps
        vnorm = np.sqrt(utils.inner_prod(f.coords, f.coords))
        qt = np.cos(dt*vnorm)*q1 + np.sin(dt*vnorm)*f/vnorm
        qt.project_b()
        geodesic_sphere.append(qt)
    #print geodesic_sphere
    return geodesic_sphere


def compute_for_closed_curves(q1, q2, settings):

    alpha = compute_on_sphere(q1, q2, settings.steps)
    alpha = project_space_closed_curves(alpha)
    alpha_t = compute_cov_derivative_path(alpha)
    alpha_pip = compute_Palais_inner_prod(alpha_t, alpha_t)
    alpha_path_len = compute_path_length(alpha_t)
    return alpha, alpha_t, alpha_pip, alpha_path_len


def project_tangent(f,q):
    return f - q * utils.inner_prod(f.coords, q.coords)


def compute_cov_derivative_path(alpha):
    n, T = np.shape(alpha[0].coords)
    stp = len(alpha) - 1
    alpha_dt = []
    alpha_dt.append(QShape(np.zeros((n, T))))
    for tau in xrange(1,len(alpha)):
        tmp = stp*(alpha[tau] - alpha[tau-1])
        alpha_dt.append(project_tangent(tmp, alpha[tau]))
    return alpha_dt


def compute_for_closed_curves_approx(q1, q2, settings):
    alpha = compute_on_sphere(q1, q2, settings.steps)
    for shape in alpha:
        shape.project_to_space_closed_curves()
    alpha_t = compute_cov_derivative_path(alpha)
    alpha_pip = compute_Palais_inner_prod(alpha_t, alpha_t)
    alpha_path_len = compute_path_length(alpha_t)
    return alpha, alpha_t, alpha_pip, alpha_path_len


def compute_for_open_curves(q1, q2, steps):

    alpha = compute_on_sphere(q1, q2, steps)
    alpha_t = compute_cov_derivative_path(alpha)
    alpha_pip = compute_Palais_inner_prod(alpha_t, alpha_t)
    alpha_path_len = compute_path_length(alpha_t)
    return alpha, alpha_t, alpha_pip, alpha_path_len


def compute_for_open_curves_elastic(q1, q2, settings, rotation=True, linear=False, gradient_decent=True):

    if rotation:
        q2, R = utils.find_best_rotation(q1, q2)

    if not linear:
        gamma = DPmatchcy.match(q1.coords, q2.coords)
    else:
        gamma = np.linspace(0, 2*pi, q1.siz())
    geodesic = Geodesic()
    geodesic.gamma = gamma*2*pi
    q2n = q2.group_action_by_gamma(geodesic.gamma)
    #alpha, alpha_t, alpha_pip, alpha_path_len = compute_for_open_curves(q1, q2n, settings.steps)
    geodesic.path, geodesic.tangent_vect, alpha_pip, geodesic.geodesic_distance \
        = compute_for_open_curves(q1, q2n, settings.steps)

    if not gradient_decent:
        #geodesic = Geodesic()
        #geodesic.path = alpha
        #geodesic.tangent_vect = alpha_t
        #geodesic.geodesic_distance = alpha_path_len
        #geodesic.gamma = gamma
        return geodesic
    else:
        V = q1.form_basis_d()
        epsilon = 0.1
        if utils.inner_prod(q1.coords - q2n.coords, q1.coords - q2n.coords) < epsilon * (1.0/15.0):
            # Default to linear matching if curves are very close together
            geodesic.gamma = np.linspace(0, 2*pi, q1.size())
            geodesic.geodesic_distance = 0
            for i in xrange(settings.steps):
                geodesic.path[:, :, i] = q1
            geodesic.tangent_vect = np.zeros(q1.dim(), q1.siz(), settings.steps+1)
            return geodesic
        i = 1
        dist_iter = []
        egeo_iter = []
        ednorm_sq_iter = []
        Anormiter = []
        while i < 45 and geodesic.geodesic_distance > epsilon * (1.0/100.0):
            q2n.project_b()
            geodesic.path, geodesic.tangent_vect, alpha_pip, geodesic.geodesic_distance\
                = compute_for_open_curves(q1, q2n, settings.steps)
            dist_iter.append(geodesic.geodesic_distance)
            dist_iter.append(alpha_pip)
            u = geodesic.path[settings.steps]
            u = u.coords
            D_q = q2n.form_basis_d_shape()
            uproj, a = utils.project_tangent_d_q(u, D_q)
            ednorm_sq = utils.inner_prod(uproj, uproj)
            ednorm_sq_iter.append(ednorm_sq)
            g = np.dot(a, V).transpose()
            s = np.linspace(0, 2*pi, q2n.siz())
            s = np.reshape(s, [s.shape[0], 1])
            Anormiter.append(np.trapz(g**2, s[0, :]))

            # Form gamma
            gamma_n = s - epsilon*g
            gamma_n -= gamma_n[0]

            #if sum(np.greater(0, gamma_n)) >= 1 or sum(np.greater(0, np.diff(gamma_n))):
            #    raise ValueError("Gamma is INVALID")
            q2n.group_action_by_gamma(gamma_n.flatten())
            q2n.project_b()
            i += 1
            if i > 1:
                geodesic.geodesic_distance = abs(dist_iter[i-1] - dist_iter[i-2])
        geodesic.gamma = q2n.estimate_gamma()
        return geodesic









def compute_for_closed_curves_elastic(q1, q2, settings, rotation=True):

    if rotation:
        q2, R = utils.find_best_rotation(q1, q2)

    gamma = DPmatchcy.match(q1.coords, q2.coords)
    gamma = gamma*2*pi
    q2n = q2.group_action_by_gamma(gamma)
    alpha, alpha_t, alpha_pip, alpha_path_len = compute_for_closed_curves_approx(q1, q2n, settings)
    geodesic = Geodesic()
    geodesic.path = alpha
    geodesic.tangent_vect = alpha_t
    geodesic.geodesic_distance = alpha_path_len
    geodesic.gamma = gamma
    return geodesic


def to_curve_path(shape_path):
    curve_path = []
    for i in xrange(0, len(shape_path)):
        curve_path.append(shape_path[i].to_curve())
    return curve_path