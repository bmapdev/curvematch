"""Mean for qshapes"""

__author__ = "Suchit Panjiyar"
__copyright__ = "Copyright 2016, Shantanu H. Joshi, Suchit Panjiyar, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "suchit@ucla.edu"


import numpy as np
import geodesics


def karcher_mean(qarray, setting):
    n_curves = len(qarray)
    q_mean = qarray[0]
    q_mean -= qarray[0]
    for i in range(0,n_curves):
        q_mean += qarray[i]
    q_mean /= n_curves
    if setting.closed:
        q_mean.project_to_space_closed_curves()
    else:
        q_mean.project_b()
    for i in range(0,6):
        alpha_t_mean = np.zeros(qarray[0].shape())
        for j in range(0,n_curves):
            if setting.closed:
                [alpha, alpha_t, alpha_pip, alpha_path_len] = geodesics.compute_for_closed_curves(q_mean, qarray[j], setting)
            else:
                [alpha, alpha_t, alpha_pip, alpha_path_len] = geodesics.compute_for_open_curves(q_mean, qarray[j], setting.steps)
            alpha_t_mean += alpha_t[1].coords
        alpha_t_mean /= n_curves
        _, q_mean, _ = geodesics.compute_flow_open(q_mean, alpha_t_mean, setting.steps)

    tangent_vectors = []
    for j in range(0,n_curves):
        if setting.closed:
            [alpha, alpha_t, alpha_pip, alpha_path_len] = geodesics.compute_for_closed_curves(q_mean, qarray[j], setting)
        else:
            [alpha, alpha_t, alpha_pip, alpha_path_len] = geodesics.compute_for_open_curves(q_mean, qarray[j], setting.steps)
        tangent_vectors.append(alpha_t[1])

    return q_mean,tangent_vectors







