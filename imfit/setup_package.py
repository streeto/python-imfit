# Licensed under GPLv3 - see LICENSE.rst

import os
from distutils.core import Extension
from astropy import setup_helpers
import subprocess

IMFITROOT = os.path.relpath(os.path.dirname(__file__))

def get_extensions():
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    cfg['sources'] = ['imfit/fit_wrapper.pyx']
    cfg['language'] = 'c++'

    # FIXME: more portable way to check the libraries.    
    libs_str = subprocess.check_output(['pkg-config', '--libs-only-l', 'imfit'])
    libs =  [l.strip() for l in libs_str.split('-l') if l != '']
    cfg['libraries'].extend(libs)

    return [Extension('imfit.fit_wrapper', **cfg)]
