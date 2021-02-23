import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import scipy
import numpy

os.environ['ARCHFLAGS'] = '-arch x86_64'
libraries = [ "m"]
extra_compile_args = []
extra_link_args = []

cython_directives = {
    'embedsignature': True,
    'boundscheck': False,
    'wraparound': False,
}

ext = Extension(
        "dfa.dfa",
        ["dfa.pyx"],
        libraries=libraries,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        include_dirs=[numpy.get_include()]
        )

setup(
    name='dfa',
    packages=['dfa'],
    version='0.0.1',
    ext_modules=cythonize([ext, ],
        compiler_directives=cython_directives,
    ),
)
