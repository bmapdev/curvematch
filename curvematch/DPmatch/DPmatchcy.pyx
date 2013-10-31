# distutils: language = c++
# distutils: sources = DPmatch.cpp

import numpy as np
cimport numpy as np
cimport cython
np.import_array()

# Wrapper class for DPmatch
# Declare constructors and matching functions
cdef extern from "DPmatch.h":
    cdef cppclass DPmatch:
        DPmatch() except +
        DPmatch(int dim, int siz) except +
        void Initialize(int dim, int siz, float *v11, float *v12, float *v13, float *v21, float *v22, float *v23)
        void MatchQ(float *)

# Package level functions
@cython.boundscheck(False)
@cython.wraparound(False)
def match(coordsarray1, coordsarray2):

    coordsarray1 = coordsarray1.astype(np.dtype('f4'))
    cdef np.ndarray[ndim=2,dtype=np.float32_t] coords1 = np.ascontiguousarray(coordsarray1)

    coordsarray2 = coordsarray2.astype(np.dtype('f4'))
    cdef np.ndarray[ndim=2,dtype=np.float32_t] coords2 = np.ascontiguousarray(coordsarray2)

    cdef int dim = coords1.shape[0]
    cdef int siz = coords1.shape[1]

    cdef DPmatch* dpmatchobj = new DPmatch(dim,siz)

    if dim == 2:
        dpmatchobj.Initialize(dim,siz,&coords1[0,0], &coords1[1,0], NULL, &coords2[0,0], &coords2[1,0], NULL)
    elif dim == 3:
        dpmatchobj.Initialize(dim,siz,&coords1[0,0], &coords1[1,0], &coords1[2,0], &coords2[0,0], &coords2[1,0], &coords2[2,0])

    cdef np.ndarray[ndim=1,dtype=np.float32_t] gamma = np.zeros(siz,np.dtype('f4'))

    dpmatchobj.MatchQ(&gamma[0])
    del dpmatchobj

    return gamma

