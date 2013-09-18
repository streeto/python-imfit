'''
Created on Sep 17, 2013

@author: andre
'''


from imfit.fit_wrapper import fit, fit_config_file
from imfit.model import ModelDescription

def read_image(fname):
    import pyfits
    arr = pyfits.getdata(fname)
    return arr.astype('float64')

image = read_image('lala/K0846_0.3.6_qSignal.fits')
noise = read_image('lala/K0846_0.3.6_qNoise.fits')
mask = read_image('lala/K0846_0.3.6_mask.fits')
psf = read_image('lala/psf_moffat36.fits')

for i in xrange(1):
    fit(image=image,
        mask=mask,
        noise=noise,
        psf=psf,
        model=ModelDescription(),
        nCombined=1,
        expTime=1.0,
        gain=1.0,
        readNoise=0.0,
        originalSky=0.0,
        ftol=1.0e-8)

    fit_config_file(image=image,
                    mask=mask,
                    noise=noise,
                    psf=psf,
                    configFileName='lala/config_sersic_K0846.dat',
                    nCombined=1,
                    expTime=1.0,
                    gain=1.0,
                    readNoise=0.0,
                    originalSky=0.0,
                    ftol=1.0e-8)
