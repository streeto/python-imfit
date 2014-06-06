'''
Created on 13/03/2014

@author: andre
'''

from imfit import Imfit, SimpleModelDescription, function_description, gaussian_psf
import numpy as np
from numpy.testing import assert_allclose


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


def get_model_param_array(model):
    params = []
    for p in model.parameterList():
        params.append(p.value)
    return np.array(params)


def test_fitting():
    psf = gaussian_psf(2.5, size=9)
    model_orig = create_model()
    imfit = Imfit(model_orig, psf=psf, quiet=False)

    noise_level = 0.1
    shape = (100, 100)
    image = imfit.getModelImage(shape)
    noise = image * noise_level
    image += (np.random.random(shape) * noise)
    mask = np.zeros_like(image, dtype='bool')
    
    imfit.fit(image, noise, mask)
    model_fitted = imfit.getModelDescription()

    orig_params = get_model_param_array(model_orig)
    fitted_params = get_model_param_array(model_fitted)

    assert_allclose(orig_params, fitted_params, rtol=noise_level)
    

if __name__ == '__main__':
    test_fitting()
    