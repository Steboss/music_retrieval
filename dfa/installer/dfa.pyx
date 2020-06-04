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

#REMEMBER TO DEFINE ALL THE C CODES, AS MULTIPLE DEFINITION ERROR WILL APPEAR OTHERWISE
cdef extern from "../c_code/dfa.c":
    void dfa(float* X,
              int n_elems,
              int* scales,
              int scales_len,
              float* fluct,
              float* coeffs)

cdef extern from "../c_code/polyfit.c":
    int polyfit(float* dependentValues,
                float* independentValues,
                int    countOfElements,
                int    order,
                float* coefficients)



cpdef play(X, scale_low, scale_high, scale_dense):
    r"""Cython transformation for dfa
    Parameters
    ----------
    For the moment just take X from python, then we'll pass scales and dens

    Returns
    -------
    """
    #get signal size
    n_elems = len(X)
    #convert to float
    X_cython = np.zeros(n_elems, dtype=np.float32)
    for i in range(n_elems):
      X_cython[i] = X[i]
    cdef float[:] a = X_cython
    #print(X_cython)
    #create the scales  intc otherwise np.int it's long which is not what we want and we'll have a type mismatch
    #scales = (2**np.arange(scale_low, scale_high, scale_dense)).astype(np.intc)
    scales = (2**np.arange(scale_low, scale_high, scale_dense)).astype(np.intc)
    print(scales)
    scale_length = len(scales)
    #convert to memory vview
    cdef int[:] memory_scales = scales#scales_cython

    #define flucuations
    cdef float[:] fluct = np.zeros(scale_length, dtype=np.float32)
    cdef float[:] coeffs = np.zeros(2, dtype=np.float32)
    print("Calling C")


    #call C-nmf
    dfa(&a[0],
        n_elems,
        &memory_scales[0],
        scale_length,
        &fluct[0],
        &coeffs[0])

    coeff_to_array = np.asarray(coeffs)
    fluct_to_array = np.asarray(fluct)
    return scales, fluct_to_array, coeff_to_array
