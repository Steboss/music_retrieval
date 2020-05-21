import scipy.io.wavfile
import matplotlib.pyplot  as plt
import sys
from libc.stdlib cimport malloc,free
from cpython cimport PyObject, Py_INCREF, array


# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()
# Import the Python-level symbols of numpy
import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np


cdef extern from "../c_code/nmf.c":
    double* nmf(int n_components, int max_iter, int random_state, double alpha,
               double *X, int X_rows, int X_columns)




cpdef play(X):
    row_size,column_size = np.shape(X)
    print(row_size, column_size)
    cdef double[:,:] a = X

    #cdef np.ndarray[double, ndim=2, mode="c"] X_cython = np.asarray(X, dtype = float, order="C")
    #Create our helper array
    #cdef double** point_to_X = <double **>malloc(row_size * sizeof(double*))
    #for i in range(row_size):
        #point_to_a[i] = &a_cython[i, 0]
        #print(X_cython[i])
        #print(X_cython[i,1])
#        point_to_X[i] = X_cython[i,0]


    nmf(2,10,10, 10.0, &a[0,0], row_size, column_size)
