'''
Created on 13/03/2014

@author: andre
'''

from imfit.fitting import Imfit
from imfit import SimpleModelDescription, function_description
from imfit.psf import gaussian_psf
import numpy as np
from numpy.testing import assert_allclose
import time
from multiprocessing import cpu_count


def create_model():
    model = SimpleModelDescription()
    model.x0.setValue(50, limits=(40, 60))
    model.y0.setValue(50, limits=(40, 60))
    
    bulge = function_description('Sersic', name='bulge')
    bulge.I_e.setValue(1.0, limits=(0.5, 1.5))
    bulge.r_e.setValue(10, limits=(5, 15))
    bulge.n.setValue(4, [3,5])
    bulge.PA.setValue(45, [30, 60])
    bulge.ell.setValue(0.5, [0, 1])
    
    disk = function_description('Exponential', name='disk')
    disk.I_0.setValue(0.7, limits=(0.4, 0.9))
    disk.h.setValue(15, limits=(10, 20))
    disk.PA.setValue(60, [45, 90])
    disk.ell.setValue(0.2, [0, 0.5])
    
    model.addFunction(bulge)
    model.addFunction(disk)
    return model


def test_model_image(imsize, nproc, count, mode='new', chunk=8):
    #psf = gaussian_psf(2.5, size=9)
    psf = None
    model_orig = create_model()
    imfit = Imfit(model_orig, psf=psf, quiet=True, nproc=nproc, chunk_size=chunk)
    shape = (imsize, imsize)
    imfit.getModelImage(shape)
    t1 = time.time()
    for _ in xrange(count):
        imfit._modelObject._testCreateModelImage(mode=mode, count=1)
    return time.time() - t1
    

if __name__ == '__main__':
    for imsize in [50, 64, 100, 128, 200, 256, 400, 512, 800, 1024, 1600, 2048]:
        for chunk in [8, 10, 16, 20, 32]:
            print 'imsize = %d, chunk = %d' % (imsize, chunk)
            with open('perf_%d_%d.dat' % (imsize, chunk), 'w') as f:
                for nproc in xrange(1, cpu_count() + 1):
                    t_new = test_model_image(imsize, nproc, 1000, mode='new', chunk=chunk)
                    t_old = test_model_image(imsize, nproc, 1000, mode='orig', chunk=chunk)
                    print nproc, t_new, t_old
                    f.write('%d %f %f\n' % (nproc, t_new, t_old))
    