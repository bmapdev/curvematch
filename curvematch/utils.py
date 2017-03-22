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
                coords[i, :] = np.interp(gamma, np.linspace(0, 2*math.pi, siz), coords[i, :])
    else:
        coords = np.interp(gamma, np.linspace(0, 2*math.pi, dim), coords[:])
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

def Form_Basis_D(d,T):
    x = np.linspace(0,1,T)
    xdarray = np.zeros((1,d))
    for i in range(d):
        xdarray[0,i] = i+1
    xdarray = np.dot(np.transpose(xdarray), x)
    V = [np.ones((1,T))*np.sqrt(2),np.sqrt(2)*np.cos(2*math.pi*xdarray), np.sqrt(2)*np.sin(2*math.pi*xdarray) ]
    return V

def Form_Basis_L2_R(d, T):
    x = np.linspace(0,1,T)
    k = 0
    B = []
    for j in range(1,d+1):
        B.append([np.sqrt(2)*np.cos(2*math.pi*j*x), np.zeros((1, T)), np.zeros((1, T))])
        B.append([np.zeros((1, T)), np.sqrt(2)*np.cos(2*math.pi*j*x), np.zeros((1, T))])
        B.append([np.zeros((1, T)), np.zeros((1, T)), np.sqrt(2)*np.cos(2*math.pi*j*x)])
        B.append([np.sqrt(2)*np.sin(2*math.pi*j*x), np.zeros((1, T)), np.zeros((1, T))])
        B.append([np.zeros((1, T)), np.sqrt(2)*np.sin(2*math.pi*j*x), np.zeros((1, T))])
        B.append([np.zeros((1, T)), np.zeros((1, T)), np.sqrt(2)*np.sin(2*math.pi*j*x)])
        k=k+1

    return B

def Form_Basis_L2_R3(d, T):
    B = []
    constB = []
    constB.append([np.sqrt(2)*np.ones((1, T)), np.zeros((1, T)), np.zeros((1, T))])
    constB.append([np.zeros((1, T)), np.sqrt(2)*np.ones((1, T)), np.zeros((1, T))])
    constB.append([np.zeros((1, T)), np.zeros((1, T)), np.sqrt(2)*np.ones((1, T))])
    B.append(constB)
    B.append(Form_Basis_L2_R(d,T))
    return B

def Form_Basis_O_q(B,q):
    d = len(B)
    #Needs to be completed

def Form_Basis_of_Tangent_Space_of_S_at_q(Bnew, G_O_q):
    G = []
    for j in range(len(Bnew)):
        tmp = 0
        for k in range(len(G_O_q)):
            tmp += np.dot(inner_prod(Bnew[j], G_O_q[k]), G_O_q[k])
        G.append(Bnew[j] - tmp)
    return G