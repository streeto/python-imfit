'''
Created on Sep 20, 2013

@author: andre
'''
from .model import ModelDescription
from .lib_wrapper import ModelObjectWrapper
import numpy as np

__all__ = ['fitter']

class fitter(object):
    
    def __init__(self, model_descr, psf=None):
        if not isinstance(model_descr, ModelDescription):
            raise ValueError('model_descr must be a ModelDescription object.')
        self._modelDescr = model_descr
#         self._paramDtype = self._getParamDtype()
        self._psf = psf
    
    
    def _getParamDtype(self):
        dt = [(p.name, 'float64') for p in self._modelDescr.parameterList()]
        return np.dtype(dt)
    
    
    @property
    def paramDtype(self):
        return self._paramDtype


    @property
    def parameters(self):
#         return np.array(self._modelObject.values(), dtype=self._paramDtype)
        return np.array(self._modelObject.values())


    def _setupModel(self):
        self._modelObject = ModelObjectWrapper(self._modelDescr)
        if self._psf is not None:
            self._modelObject.setPSF(np.asarray(self._psf))
            
    
    def __call__(self, image, noise, mask=None):
        self._setupModel()
        if isinstance(image, np.ma.MaskedArray):
            image = image.filled()
            if mask is None:
                mask = image.mask.astype('float64')
        if mask is None:
            mask = np.ones_like(image)
        else:
            mask = mask.astype('float64')
        if isinstance(noise, np.ma.MaskedArray):
            noise = noise.filled()

        self._modelObject.setData(image, noise, mask,
                                  n_combined=1, exp_time=1.0, gain=1.0, read_noise=0.0, original_sky=0.0)
        self._modelObject.fit()
        