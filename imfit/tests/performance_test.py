'''
Created on 13/03/2014

@author: andre
'''

from imfit import Imfit, SimpleModelDescription, function_description
import time
from multiprocessing import cpu_count


def create_model():
    model = SimpleModelDescription()
    model.x0.setValue(50, vmin=40, vmax=60)
    model.y0.setValue(50, vmin=40, vmax=60)
    
    bulge = function_description('Sersic', name='bulge')
    bulge.I_e.setValue(1.0, vmin=0.5, vmax=1.5)
    bulge.r_e.setValue(10, vmin=5, vmax=15)
    bulge.n.setValue(4, vmin=3, vmax=5)
    bulge.PA.setValue(45, vmin=30, vmax=60)
    bulge.ell.setValue(0.5, vmin=0, vmax=1)
    
    disk = function_description('Exponential', name='disk')
    disk.I_0.setValue(0.7, vmin=0.4, vmax=0.9)
    disk.h.setValue(15, vmin=10, vmax=20)
    disk.PA.setValue(60, vmin=45, vmax=90)
    disk.ell.setValue(0.2, vmin=0, vmax=0.5)
    
    model.addFunction(bulge)
    model.addFunction(disk)
    return model


def test_model_image(imsize, nproc, count, chunk=8):
    psf = None
    model_orig = create_model()
    imfit = Imfit(model_orig, psf=psf, quiet=True, nproc=nproc, chunk_size=chunk)
    shape = (imsize, imsize)
    imfit.getModelImage(shape)
    t1 = time.time()
    imfit._modelObject._testCreateModelImage(count)
    return time.time() - t1
    

if __name__ == '__main__':
    for imsize in [50, 64, 100, 128, 200, 256, 400, 512, 800, 1024, 1600, 2048]:
        for chunk in [8, 10, 16, 20, 32]:
            print 'imsize = %d, chunk = %d' % (imsize, chunk)
            with open('perf_%d_%d.dat' % (imsize, chunk), 'w') as f:
                for nproc in xrange(1, cpu_count() + 1):
                    t = test_model_image(imsize, nproc, 1000, chunk=chunk)
                    print nproc, t
                    f.write('%d %f\n' % (nproc, t))
    