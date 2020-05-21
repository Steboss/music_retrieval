import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import scipy
import numpy

#mdgxhome = os.environ.get("M")

#if mdgxhome is not None:
#    library_dirs = [os.path.join(mdgxhome, 'lib'),]
#    include_dirs = [os.path.join(mdgxhome, 'include'),]
#else:
#    amberhome = os.environ.get('AMBERHOME')
#    library_dirs = [os.path.join(amberhome, 'lib'),]
#    include_dirs = [os.path.join(amberhome, 'include'),]
libraries = [ "fftw3","m"]
extra_compile_args = ['-O3', '-std=c99']
extra_link_args = ['-O3', '-std=c99']

cython_directives = {
    'embedsignature': True,
    'boundscheck': False,
    'wraparound': False,
}

ext = Extension(
        "stft.stft",
        ["stft.pyx",],
        libraries=libraries,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        )

setup(
    name='stft',
    packages=['stft'],
    version='0.0.1',
    ext_modules=cythonize([ext, ],
        compiler_directives=cython_directives,
    ),
)
