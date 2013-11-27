#!/usr/local/epd/bin/python

"""Declaration of the qshape class"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from math import pi

import numpy as np
import scipy as sp
from scipy import linalg as LA
from scipy import integrate

import utils

class QShape():

    def __init__(self, coords=[], dim=0, siz=0):
        n, T = coords.shape
        if n > T:
            self.coords = np.transpose(coords)
        else:
            self.coords = coords
        self.dim = dim
        self.siz = siz
        self.shape = (self.dim,self.siz)
        if coords is not None:
            self.dim = self.coords.shape[0]
            self.siz = self.coords.shape[1]
            self.shape = (self.dim,self.siz)

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

    def from_curve(self, curve):
        self.dim,self.siz = curve.dim, curve.siz
        #n = shape_p[0]
        #Tcoord = shape_p[1]
        # Compute the gradient of all row vectors in p
        pdiff = np.zeros(self.shape)
        for i in xrange(self.dim):
            pdiff[i,:] = np.gradient( curve.coords[i,:], 2*pi /self.siz )
        v = np.zeros(self.shape)
        for i in xrange(self.dim):
            v[i,:] = ( 2 * pi / self.siz ) * pdiff[i,:]
        q = np.zeros(self.shape)
        for i in xrange(self.siz):
            q[:,i] = v[:,i]/ np.sqrt(LA.norm(v[:,i]))

        self.coords = self.ProjectC(q)
        return(q)


    def to_curve(self,q):  #Return a curve object????
        s = np.linspace(0,2*pi,self.siz)
        qnorm = np.zeros(self.siz)
        for i in xrange(self.siz):
            qnorm[i] = LA.norm(q[:,i],2)
        p = np.zeros(np.shape(q))
        for i in xrange(self.dim):
            temp = q[i,:] * qnorm
            p[i,:] = integrate.cumtrapz(temp,s,initial=0)
        return(p)


    def project_B(self):
        self.coords = self.coords/np.sqrt(utils.inner_prod(self, self))
        #qnew = q/np.sqrt(utils.inner_prod(q,q))
        #return qnew


    def projectC(self,q):
        shape_q = np.shape(q)
        n = shape_q[0]
        T = shape_q[1]
        dt = 0.3
        epsilon = (1.0/60.0) * 2.0 * (pi / T )

        itr = 0
        res = np.ones((1,n))
        J = np.zeros( (n,n) )
        s = np.linspace(0,2*pi,T)

        qnew = q / np.sqrt(utils.inner_prod(q,q))
        C = {}
        while LA.norm(res) > epsilon:
            if itr > 300:
                print "Warning: Shape failed to project.  Geodesics will be incorrect."
                break

            # Compute Jacobian
            for i in xrange(n):
                for j in xrange(n):
                    J[i,j] = 3 * np.trapz( qnew[i,:] * qnew[j,:] , s)
            J += np.eye(n)

            qnorm = np.zeros(T)
            for i in xrange(T):
                qnorm[i] = LA.norm(qnew[:,i])
                ##################################################
            # Compute the residue
            G = np.zeros(n)
            for i in xrange(n):
                G[i] = np.trapz( (qnew[i,:] * qnorm) , s)
            res = -G
            C[itr] = LA.norm(res)
            cond_J = LA.cond(J)
            if LA.norm(res) < epsilon:
                qnew = qnew/np.sqrt(utils.inner_prod(qnew,qnew))
                return(qnew)

            if  np.isnan(cond_J) or np.isinf(cond_J) or (cond_J < 0.1):
                print '\nProjection may not be accurate\n'
                qnew = q
                qnew = qnew/np.sqrt(utils.inner_prod(qnew,qnew))
                return(qnew)
            else:
                x = LA.solve(J ,res.T)
                delG = self.form_basis_normal_a(qnew)
                temp = 0
                for i in xrange(n):
                    temp = temp + x[i] * delG[i] * dt
                qnew += temp

                itr += 1  #Iterator for while loop
        qnew = qnew/np.sqrt(utils.inner_prod(qnew,qnew))
        print("What?")

        return(qnew)


    def form_basis_normal_a(self,q):
        shape_q = np.shape(q)
        n = shape_q[0]
        T = shape_q[1]
        e = np.eye(n)

        Ev = np.zeros((n,T,n))
        for i in xrange(n):
            Ev[:,:,i] = np.tile(e[:,i],(T,1) ).transpose()
            qnorm = np.zeros(T)
        for i in xrange(T):
            qnorm[i] = LA.norm(q[:,i])
        delG = {}
        for i in xrange(n):
            tmp1 = np.tile( (q[i,:] / qnorm) ,(n,1))
            tmp2 = np.tile( qnorm ,(n,1) )
            delG[i] = tmp1*q + tmp2*Ev[:,:,i]
        return delG


    def Group_Action_by_Gamma(self,q,gamma):
        shape_q = np.shape(q)
        n,T = shape_q
        #gamma_orig = Estimate_Gamma(q)
        gamma_t = np.gradient(gamma,2*pi/(T-1))
        q_composed_gamma = np.zeros(shape_q)
        for i in xrange(n):
            q_composed_gamma[i,:] = np.interp(np.linspace(0,2*pi,T),q[i,:],gamma)  #Can not specify "nearest" method
        sqrt_gamma_t = np.tile(np.sqrt(gamma_t),(5,1))
        qn = q_composed_gamma * sqrt_gamma_t;
        return(qn)


    def Estimate_Gamma(self,q):
        p = self.to_curve(q)
        shape_q = np.shape(p)
        n,T = shape_q

        #Evaluate the arc-length function
        pdiff = np.diff(p,1,1)
        ds = np.sqrt(np.sum(pdiff**2,1)) * T
        gamma = integrate.cumsum(ds)*2*pi/max(integrate.cumsum(ds))
        return(gamma)



