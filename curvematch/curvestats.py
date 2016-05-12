"""Curve Statistics"""

__author__ = "Suchit Panjiyar"
__copyright__ = "Copyright 2016, Shantanu H. Joshi, Suchit Panjiyar, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "suchit@ucla.edu"

import numpy as np
from scipy import linalg
from scipy import stats
import utils


def pca(qmean, tangent_vectors, odir):
    epsilon = 0.0001
    Y = utils.Gram_Schmidt(tangent_vectors)

    Xproj = utils.Project_To_Basis(tangent_vectors,Y)
    Xproj = np.array(Xproj)
    C = np.cov(Xproj.T)
    [U, S, V] = linalg.svd(C)

    sDiag = np.diag(S)

    tmp = np.identity(len(S))
    tmp = epsilon*tmp

    Cn = U*(tmp+sDiag)*U.T

    ret_dict = {}
    ret_dict['Cn'] = Cn
    ret_dict['qmean'] = qmean
    ret_dict['alpha_t_array'] = tangent_vectors
    ret_dict['Y'] = Y
    ret_dict['Xproj'] = Xproj
    ret_dict['U'] = U
    ret_dict['S'] = S
    ret_dict['V'] = V
    ret_dict['C'] = C

    return ret_dict



# covdata.qmean = qmean;
# covdata.alpha_t_array = alpha_t_array;
# covdata.Y = Y;
# covdata.Xproj = Xproj;
# covdata.U = U;
# covdata.S = S;
# covdata.V = V;
# covdata.C = C;
# covdata.Cn = Cn;

