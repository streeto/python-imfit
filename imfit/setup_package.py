# Licensed under GPLv3 - see LICENSE.rst

import os
from glob import glob
from distutils.core import Extension
from astropy import setup_helpers

IMFITROOT = os.path.relpath(os.path.dirname(__file__))

def get_extensions():
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    cfg['sources'] = ['imfit/fit_wrapper.pyx']
    cfg['libraries'].extend(['imfit', 'fftw3', 'fftw3_threads', 'gsl'])
    cfg['language'] = 'c++'

    return [Extension('imfit.fit_wrapper', **cfg)]
