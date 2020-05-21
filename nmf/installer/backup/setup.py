import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import scipy
import numpy

libraries = [ "m"]
extra_compile_args = ['-O3', '-std=c99']
extra_link_args = ['-O3', '-std=c99']

cython_directives = {
    'embedsignature': True,
    'boundscheck': False,
    'wraparound': False,
}

ext = Extension(
        "nmf.nmf",
        ["nmf.pyx",],
        libraries=libraries,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        )

setup(
    name='nmf',
    packages=['nmf'],
    version='0.0.1',
    ext_modules=cythonize([ext, ],
        compiler_directives=cython_directives,
    ),
)
