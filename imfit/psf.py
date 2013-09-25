'''
Created on Sep 25, 2013

@author: andre
'''
from .fitting import Imfit
from .model import SimpleModelDescription
from .lib_wrapper import function_description
import numpy as np

__all__ = ['gaussian_psf', 'moffat_psf']


FWHM_to_sigma_factor = 2.0 * np.sqrt(2.0 * np.log(2.0))


def gaussian_psf(width, type='fwhm', PA=0.0, ell=0.0, size=31):
    if size % 2 != 1:
        raise ValueError('Size must be an odd number.')
    if type == 'fwhm':
        sigma = width / FWHM_to_sigma_factor
    elif type == 'sigma':
        sigma = width
    else:
        ValueError('type must be either fwhm or sigma.')
    center = (size+1) / 2
    model = SimpleModelDescription()
    model.x0.setValue(center)
    model.y0.setValue(center)
    gaussian = function_description('Gaussian')
    gaussian.I_0.setValue(1.0)
    gaussian.sigma.setValue(sigma)
    gaussian.PA.setValue(PA)
    gaussian.ell.setValue(ell)
    model.addFunction(gaussian)
    imfit = Imfit(model)
    psf_image = imfit.getModelImage(shape=(size,size))
    return psf_image


def moffat_psf(fwhm, beta=3.1, PA=0.0, ell=0.0, size=31):
    if size % 2 != 1:
        raise ValueError('Size must be an odd number.')
    center = (size+1) / 2
    model = SimpleModelDescription()
    model.x0.setValue(center)
    model.y0.setValue(center)
    moffat = function_description('Moffat')
    moffat.I_0.setValue(1.0)
    moffat.fwhm.setValue(fwhm)
    moffat.beta.setValue(beta)
    moffat.PA.setValue(PA)
    moffat.ell.setValue(ell)
    model.addFunction(moffat)
    imfit = Imfit(model)
    psf_image = imfit.getModelImage(shape=(size,size))
    return psf_image
