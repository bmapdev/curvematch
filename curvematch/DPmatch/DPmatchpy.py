#!/usr/local/epd/bin/python
"""Dynamic Programming for matching qshapes"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from qshape import QShape
import numpy as np
from sys import stdout


def matchq( q1,  q2):
    gamma = []
    Energy = []

    n1 = q1.n
    T1 = q1.T
    n2 = q2.n
    T2 = q2.T
    if n1 != n2 or T1 != T2:
        stdout.write('\nDimension mismatch. Shape coordinates should have same length.\n')
        return

    Nbrs = np.array([[1, 1], [1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1], [1, 4], [3, 4], [4, 3], [4, 1], [1, 5], [2, 5],
                    [3, 5], [4, 5], [5, 4], [5, 3], [5, 2], [5, 1], [1, 6], [5, 6], [6, 5], [6, 1], [1, 7], [2, 7], [3, 7], [4, 7], [5, 7], [6, 7],
                    [4, 10], [4, 30], [4, 40], [4, 50], [5, 50],[6, 50],[8, 50],[10, 50],[15, 50],[20, 50],[25, 50], [30, 50], [35, 50], [40, 50]])

    NBR_SIZ = Nbrs.shape[0]
    T = q1.T

    Energy = np.zeros(shape=(T, T), dtype=float, order='C')
    Energy[0,:] = 5
    Energy[:,0] = 5
    Energy[0,0] = 0

    Path_x = np.zeros(shape=(T, T), dtype=float, order='C')
    Path_y = np.zeros(shape=(T, T), dtype=float, order='C')
    CandE = np.zeros(shape=(NBR_SIZ,),dtype=float)
    xx1 = np.arange(0,T)
    x = np.zeros(shape=(T,),dtype=float)
    y = np.zeros(shape=(T,),dtype=float)
    xnew = np.zeros(shape=(T,),dtype=float)
    ynew = np.zeros(shape=(T,),dtype=float)

    for i in np.arange(1,T):
        stdout.write("Index: %d   \r" % i)
        stdout.flush()
        for j in np.arange(1,T):
            minCandE = 10000
            for Num in np.arange(0,NBR_SIZ):
                k = i - Nbrs[Num][0]
                l = j - Nbrs[Num][1]
                if k >= 0 and l >= 0:
                    CandE[Num] = Energy[k,l] + match_costq(q1, q2,k,l,i,j)
                else:
                    CandE[Num] = 10000

                if CandE[Num] < minCandE:
                    minCandE = CandE[Num]
                    minCandE_idx = Num
                Energy[i][j] = minCandE

                Path_x[i][j] = i - Nbrs[minCandE_idx][0]
                Path_y[i][j] = j - Nbrs[minCandE_idx][1]

    # Displaying the energies and the minimum path
    x[0] = T-1; y[0] = T-1
    cnt = 0
    while x[cnt] > 0:
        y[cnt+1] = Path_x[y[cnt],x[cnt]]
        x[cnt+1] = Path_y[y[cnt],x[cnt]]
        cnt += 1

    x = x[0:cnt+1]
    y = y[0:cnt+1]
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]

    gamma = np.interp(xx1,x,y)
    gamma = (gamma - gamma[0])/(gamma[-1] - gamma[0])
    return gamma, Energy


def match_costq(q1, q2,k,l,i,j):

    n1 = q1.n
    T1 = q1.T

    slope = (float)(i - k)/(j - l)
    if slope == 0:
        stdout.write("\nslope zero\n")

    tmp = 0
    Energy = 0
    for x in np.arange(l,j+1):
        y = k + (x-l)*slope
        y1 = np.floor(y)
        y2 = np.ceil(y)
        f = y - y1

        for kk in np.arange(0,n1):
            tmp = (f*q2.coords[kk,y2] + (1 - f)*q2.coords[kk,y1])*np.sqrt(slope)
            Energy += (q1.coords[kk,x] - tmp)*(q1.coords[kk,x] - tmp)

    Energy = Energy/T1

    return Energy
