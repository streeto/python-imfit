# Licensed under GPLv3 - see LICENSE.rst

import os

from distutils.core import Extension

IMFITROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [os.path.join(IMFITROOT, 'test_extlib.pyx')]
    libraries = ['imfit', 'fftw3', 'fftw3_threads']

    test_ext = Extension(
        name='imfit.test_extlib',
        sources=sources,
        libraries=libraries,
        language='c++')

    return [test_ext]
