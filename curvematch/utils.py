#!/usr/local/epd/bin/python

"""Declaration of the Utility module"""

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import numpy as np
import math
from numpy import linalg as LA
import qshape

def inner_prod(coords1, coords2):
    value = np.trapz(sum(coords1*coords2), np.linspace(0, 2*math.pi, coords1.shape[1]))
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
    U, S, V_transpose = LA.svd(A) # linalg.svd actually returns U, S, and transpose(V)
    if np.absolute(LA.det(U)*LA.det(V_transpose.transpose()) - 1) < 10*np.spacing(1):
        S = np.eye(q1.dim())
    else:
        S = np.eye(q1.dim())
        S[:, -1] = -S[:, -1]
    R = np.dot(np.dot(U, S), V_transpose)  # R=U*S*V' (matrix multiplication)
    q2new = q2.copy()
    q2new.coords = np.dot(R, q2.coords)
    return q2new, R


def reparameterize_by_gamma(coords, gamma):
    dim = np.shape(coords)[0]
    if len(np.shape(coords)) > 1:
        siz = np.shape(coords)[1]
    else:
        siz = 0
    if siz > 0:
        for i in xrange(0, dim):
                coords[i, :] = np.interp(gamma, np.linspace(0, 2*pi, siz), coords[i, :])
    else:
        coords = np.interp(gamma, np.linspace(0, 2*pi, dim), coords[:])
    return coords

def project_tangent_d_q(u, d_q):
    n, T, d = d_q.shape
    uproj = 0
    a = np.zeros(d_q.shape[2])
    for i in xrange(d):
        a[i] = inner_prod(u, d_q[:, :, i])
        uproj += a[i] * d_q[:, :, i]
    a = np.reshape(a, (1, a.shape[0]))
    return uproj, a


def project_tangent_q(f, q):
    fnew = f - inner_prod(f, q.coords)*q.coords
    return fnew

def Gram_Schmidt(X):
    epsilon = 0.00005
    columns = len(X)
    i = 0
    r = 0
    Y = []
    Y.append(X[0])
    while i < columns:
        tempvect = 0
        for j in range(0,i):
            tempvect = (Y[j]*inner_prod(Y[j].coords, X[r].coords)) + tempvect

        Atemp = X[r] - tempvect
        if i == 0:
            Y[0] = Atemp
        else:
            Y.append(Atemp)
        tmp = inner_prod(Y[i].coords,Y[i].coords)
        if tmp>epsilon:
            Y[i] = Y[i]/math.sqrt(tmp)
            i = i+1
            r = r+1
        else:
            if r < i:
                r = r+1
            else:
                break
    return Y

def Project_To_Basis(alpha_t_array,Y):
    n,T = Y[0].shape()
    Xproj = []
    for i in range(0,len(alpha_t_array)):
        arr = []
        for j in range(0,len(Y)):
            temp = inner_prod(alpha_t_array[i].coords, Y[j].coords)
            arr.append(temp)
        Xproj.append(arr)
    return Xproj



'''
function [Xproj,X] = Project_To_Basis(alpha_t_array,Y)

[n,T] = size(Y{1});
for i = 1:length(alpha_t_array);
    X{i} = zeros(n,T);
    for j = 1:length(Y)
        Xproj(i,j) = InnerProd_Q(alpha_t_array{i},Y{j});
        X{i} = X{i} + Xproj(i,j) * Y{j};
    end
end

return;
'''