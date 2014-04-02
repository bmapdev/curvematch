#!/usr/local/epd/bin/python

"""Declaration of the qshape class"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from math import pi
import numpy as np
from numpy import linalg as LA
from scipy import integrate
import curve
import utils
import copy

class QShape():

    def __init__(self, coords=np.array([])):

        if coords.size != 0:
            if coords.shape[0] > coords.shape[1]:
                self.coords = np.transpose(coords)
            else:
                self.coords = coords

    def dim(self):
        return self.coords.shape[0]

    def siz(self):
        if len(self.coords.shape) < 2:
            return 0
        else:
            return self.coords.shape[1]

    def shape(self):
        if len(self.coords.shape) > 1:
            return self.coords.shape
        else:
            return 0, 0


    def __add__(self, q):
        return QShape(self.coords + q.coords)

    def __sub__(self, q):
        return QShape(self.coords - q.coords)

    def __mul__(self, x):
        return QShape(self.coords * x)

    def __rmul__(self, x):
        return QShape(self.coords * x)

    def __div__(self, x):
        return QShape(self.coords / x)

    def copy(self):
        return QShape(copy.copy(self.coords))

    def from_curve_file(self, filename):
        curve1 = curve.Curve()
        curve1.readcurve(filename)

        self.from_curve(curve1)

    def from_curve(self, curve):

        dim, siz = curve.dim(), curve.siz()
        if dim > siz:
            coords = np.transpose(curve.coords)
            dim, siz = coords.shape
        else:
            coords = curve.coords

        pdiff = np.zeros(coords.shape)
        for i in xrange(dim):
            pdiff[i, :] = np.gradient(coords[i, :], 2*pi / siz)
        v = np.zeros(coords.shape)
        for i in xrange(dim):
            v[i, :] = (2 * pi / siz) * pdiff[i, :]
        coords = np.zeros(coords.shape)
        for i in xrange(siz):
            coords[:, i] = v[:, i] / np.sqrt(LA.norm(v[:, i], ord=2))

        self.coords = coords
        self.project_b()
        #self.ProjectC(q)

    def to_curve(self):  # Return a curve object? Should we import curve?
        s = np.linspace(0, 2*pi, self.siz())
        qnorm = np.zeros(self.siz())
        for i in xrange(self.siz()):
            qnorm[i] = LA.norm(self.coords[:,i ], ord=2)
        p = curve.Curve(np.zeros(self.shape()), [], self.dim(), self.siz())
        for i in xrange(self.dim()):
            temp = self.coords[i, :] * qnorm
            p.coords[i, :] = integrate.cumtrapz(temp, s, initial=0)
        return p

    def project_b(self):
        self.coords = self.coords/(np.sqrt(utils.inner_prod(self.coords, self.coords)) + np.spacing(1))
        #qnew = q/np.sqrt(utils.inner_prod(q,q))
        #return qnew

    def project_to_space_closed_curves(self):
        dt = 0.3
        epsilon = (1.0/60.0) * (2.0*pi / self.siz())

        itr = 0
        res = np.ones((1, self.dim()))
        J = np.zeros((self.dim(), self.dim()))

        s = np.linspace(0, 2*pi, self.siz())

        qnew = QShape()
        qnew.coords = self.coords/(np.sqrt(utils.inner_prod(self.coords, self.coords)) + np.spacing(1))
        C = []
        while LA.norm(res, ord=2) > epsilon:
            if itr > 300:
                print "Warning: Shape failed to project.  Geodesics will be incorrect."
                self = qnew
                return  #(qnew)

            # Compute Jacobian
            for i in xrange(self.dim()):
                    for j in xrange(self.dim()):
                        J[i, j] = 3 * np.trapz(qnew.coords[i, :] * qnew.coords[j, :], s)
            J += np.eye(self.dim())

            qnorm = np.zeros(self.siz())
            for i in xrange(self.siz()):
                qnorm[i] = LA.norm(qnew.coords[:, i], ord=2)
            ##################################################
            # Compute the residue
            G = np.zeros(self.dim())
            for i in xrange(self.dim()):
                G[i] = np.trapz((qnew.coords[i, :] * qnorm), s)
            res = -G

            if LA.norm(res, ord=2) < epsilon:
                self = qnew
                return  # qnew

            #C.append( LA.norm(res, ord=2))
            cond_J = LA.cond(J)

            if np.isnan(cond_J) or np.isinf(cond_J) or (cond_J < 0.1):
                print '\nProjection may not be accurate\n'
                #qnew = q
                #qnew = q / np.sqrt(InnerProd_Q(q, q))
                self.coords = self.coords/(np.sqrt(utils.inner_prod(self.coords, self.coords)) + np.spacing(1))
                return #qnew
            else:
                x = LA.solve(J, res.T)
                delG = self.form_basis_normal_a()
                temp = 0
                for i in xrange(self.dim()):
                    temp = temp + x[i] * delG[i] * dt
                qnew.coords += temp
                itr += 1  #Iterator for while loop

    def form_basis_normal_a(self):
        e = np.eye(self.dim())
        ev = []
        for i in xrange(self.dim()):
            ev.append( np.tile(e[:, i], (self.siz(), 1)).transpose())
            qnorm = np.zeros(self.siz())
        for i in xrange(self.siz()):
            qnorm[i] = LA.norm(self.coords[:, i], ord=2)
        del_g = {}
        for i in xrange(self.dim()):
            tmp1 = np.tile((self.coords[i, :] / qnorm), (self.dim(), 1))
            tmp2 = np.tile(qnorm, (self.dim(), 1))
            del_g[i] = tmp1*self.coords + tmp2*ev[i]
        return del_g

    def group_action_by_gamma(self, gamma):
        #gamma_orig = Estimate_Gamma(q)
        gamma_t = np.gradient(gamma, 2*pi/self.siz())
        q_composed_gamma = np.zeros(self.shape())
        for i in xrange(self.dim()):
            q_composed_gamma[i, :] = np.interp(gamma, np.linspace(0, 2*pi, self.siz()), self.coords[i, :])  # Can not specify "nearest" method
        sqrt_gamma_t = np.tile(np.sqrt(gamma_t), (self.dim(), 1))  # Possible error?
        return QShape(q_composed_gamma * sqrt_gamma_t)

    def estimate_gamma(self):
        p = self.to_curve()
        #Evaluate the arc-length function
        s = np.linspace(0,2*pi,self.siz)
        pgrad = np.gradient(p.coords, 2*pi/np.siz)
        ds = np.sqrt(np.sum(pgrad**2, 0)) * np.siz
        gamma = integrate.cumtrapz(ds, s) * 2*pi/np.max(integrate.cumtrapz(ds,s))
        return gamma



