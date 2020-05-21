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


cdef extern from "../c_code/nmf_2.c":
    double* nmf(int n_components,
                int max_iter,
                int random_state,
                double *X,
                int X_rows,
                int X_columns,
                double *W_cython,
                double *H_cython,
                int verbose)




cpdef play(X, n_components, max_iter, random_state, verbose):
    r"""Cython transformation for nmf
    Parameters
    ----------
    X:  numpy array
        input signal, like magnitude from stft

    n_components: int
        number of components nmf has to compute

    max_iter: int
        maximum number of iterations

    random_state: int
        random state for reproducibility

    verbose: int
        if verbose or not

    Returns
    -------
    """
    #get signal size
    row_size,column_size = np.shape(X)
    #convert the input array in a memory view
    print("Conversion of X into memoryview")
    cdef double[:,:] a = X
    #TODO: What about if we have a high dimensional signal?????
    print("Creating W array with %d elements"% (n_components*row_size))
    #initialize a vectorize version of W
    W_elems = n_components*row_size
    cdef double[:] W_cython = np.zeros(W_elems)
    #intiialize a vectorize version of H
    print("Creating H array with %d elements"%(n_components*column_size))
    H_elems = n_components*column_size
    cdef double[:] H_cython = np.zeros(H_elems)
    print("Calling C")
    #call C-nmf
    nmf(n_components,
        max_iter,
        random_state,
        &a[0,0],
        row_size,
        column_size,
        &W_cython[0],
        &H_cython[0],
        verbose)
    #translate the pointer to a numpy array
    W_array =np.asarray(W_cython)
    W_array = np.reshape(W_array, [row_size, n_components])
    #print(W_array)
    H_array =np.asarray(H_cython)
    H_array = np.reshape(H_array, [n_components, column_size])
    #print(H_array)
    return W_array, H_array
