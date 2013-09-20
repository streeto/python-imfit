'''
Created on Sep 17, 2013

@author: andre
'''


from imfit.lib_wrapper import ModelObjectWrapper, fit_wrapper
from imfit.model import ModelDescription

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
    


def fit2():
    image = read_image('imfit/tests/data/K0846_0.3.6_qSignal.fits')
    noise = read_image('imfit/tests/data/K0846_0.3.6_qNoise.fits')
    mask = read_image('imfit/tests/data/K0846_0.3.6_mask.fits')
    psf = read_image('imfit/tests/data/psf_moffat36.fits')
    model_desc = ModelDescription.load('imfit/tests/data/config_sersic_K0846.dat')
    fit_wrapper(image=image,
                mask=mask,
                noise=noise,
                psf=psf,
                model=model_desc,
                nCombined=1,
                expTime=1.0,
                gain=1.0,
                readNoise=0.0,
                originalSky=0.0,
                ftol=1.0e-8,
                verbose=1)



for i in xrange(1):
    print '#'*50
    print 'Model wrapper:'
    fit1()
    print '#'*50
    print 'fit_wrapper:'
    fit2()
