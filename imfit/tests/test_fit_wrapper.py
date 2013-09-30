'''
Created on Sep 17, 2013

@author: andre
'''


from imfit import ModelDescription
from imfit import Imfit
from imfit.lib_wrapper import ModelObjectWrapper 
from imfit.psf import moffat_psf

def read_image(fname):
    import pyfits
    arr = pyfits.getdata(fname)
    return arr.astype('float64')


def fit1():
    image = read_image('imfit/tests/data/K0846_0.3.6_qSignal.fits')
    noise = read_image('imfit/tests/data/K0846_0.3.6_qNoise.fits')
    mask = read_image('imfit/tests/data/K0846_0.3.6_mask.fits')
    psf = read_image('imfit/tests/data/psf_moffat36.fits')
    model_desc = ModelDescription.load('imfit/tests/data/config_sersic_K0846.dat')
    model = ModelObjectWrapper(model_desc)
    model.setPSF(psf)
    model.setData(image, noise, mask,
                  n_combined=1,
                  exp_time=1.0,
                  gain=1.0,
                  read_noise=0.0,
                  original_sky=0.0)
    for i in xrange(1):
        print '#'*10, i
        model.fit(verbose=1)
        print 'chi2 = %f' % model.getFitStatistic(mode='chi2')
        print 'Reduced chi2 = %f' % model.getFitStatistic(mode='reduced_chi2')
        print 'AIC = %f' % model.getFitStatistic(mode='AIC')
        print 'BIC = %f' % model.getFitStatistic(mode='BIC')
        print 'FEV = %d' % model.nFev
        print 'Iterations = %d' % model.nIter
        print 'Pegged params = %d' % model.nPegged
        print 'Error? %s' % model.error
        print 'Converged? %s' % model.converged


def fit3():
    image = read_image('imfit/tests/data/K0846_0.3.6_qSignal.fits')
    noise = read_image('imfit/tests/data/K0846_0.3.6_qNoise.fits')
    mask = read_image('imfit/tests/data/K0846_0.3.6_mask.fits')
    psf = moffat_psf(3.6, size=51)
    model_desc = ModelDescription.load('imfit/tests/data/config_sersic_K0846.dat')
    my_imfit = Imfit(model_desc, psf)
    for i in xrange(1):
        print '#'*10, i
        my_imfit.fit(image, noise, mask)
        print my_imfit.getRawParameters()
        print my_imfit.getModelDescription()
    

for i in xrange(1):
    print '#'*50
    print 'Model wrapper:'
    fit1()
    print '#'*50
    print 'fitter:'
    fit3()
