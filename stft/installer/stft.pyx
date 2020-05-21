import scipy.io.wavfile
import matplotlib.pyplot  as plt
import sys
from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF


# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()


cdef extern from "c_code/stft.c":
    double* stft(double *wav_data, int samples, int windowSize, int hop_size,\
                 double *magnitude, int sample_freq, int length)

cpdef play(audData, rate,windowSize,hopSize):

    #create a memory view, pointer, that can be processed in C
    length = len(audData)
    cdef double[:] audData_view = audData
    #the total number of elements to retrieve from the C code is:
    n_elements = int( (length/(windowSize/2))*((windowSize/2)+1) )
    print("Total number of samples %d"% n_elements)

    #Create a view for the magnitude
    cdef double[:] magnitude = np.zeros(n_elements)

    #test with a complex element fftw_complex *stft_data,
    #cdef fftw_complex[:] test = np.zeroes(n_elements)


    #view for the frequencies
    n_freqs = int(windowSize/2)+1
    #cdef double[:] frequencies = np.zeros(n_freqs)
    print("Total number of freq in Cython %d\n" %(n_freqs))
    print("Analysis with these parameters:")
    print("Total number of samples %d" % n_elements)
    print("Windows length %d" % windowSize)
    print("Hopping size %d" % hopSize)

    #compute the stft
    stft(&audData_view[0], n_elements, windowSize, hopSize,&magnitude[0],rate, length)

    #translate the pointer to a numpy array
    magnitude_array = np.asarray(magnitude)
    #frequencies_array = np.asarray(frequencies)

    #return magnitude_array,frequencies_array
    #try to implement the python part here

    cols = int( (length/(windowSize/2)) -1)
    rows = int(windowSize/2)+1
    print(cols, rows)
    new_array = np.zeros([cols,rows])
    counter = 0
    for i in range(0,cols):
        for j in range(0,rows):
            new_array[i][j] = magnitude_array[counter]
            counter+=1
            #print(counter)

    return  new_array
