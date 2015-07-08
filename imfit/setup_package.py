# Licensed under GPLv3 - see LICENSE.rst

from distutils.core import Extension
import subprocess

def get_extensions():
    cfg = {}
    cfg['include_dirs'] = ['numpy']
    cfg['sources'] = ['imfit/lib/lib_wrapper.pyx']
    cfg['language'] = 'c++'

    # FIXME: more portable way to check the libraries.    
    libs_str = subprocess.check_output(['pkg-config', '--libs-only-l', 'imfit'])
    libs =  [l.strip() for l in libs_str.split('-l') if l != '']
    cfg['libraries'] = libs

    return [Extension('imfit.lib.lib_wrapper', **cfg)]
