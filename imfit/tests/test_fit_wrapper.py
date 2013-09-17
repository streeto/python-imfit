'''
Created on Sep 17, 2013

@author: andre
'''


from imfit.fit_wrapper import fit
import pyfits
image = pyfits.getdata('lala/K0846_0.3.6_qSignal.fits').astype('float64')
noise = pyfits.getdata('lala/K0846_0.3.6_qNoise.fits').astype('float64')
mask = pyfits.getdata('lala/K0846_0.3.6_mask.fits').astype('float64')
psf = pyfits.getdata('lala/psf_moffat36.fits').astype('float64')

for i in xrange(1):
    fit(image, mask, noise, psf, 'lala/config_sersic_K0846.dat', 1, 1.0, 1.0, 0.0, 0.0, 1.0e-8)

