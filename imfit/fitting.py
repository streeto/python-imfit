'''
Created on Sep 20, 2013

@author: andre
'''
from .model import ModelDescription
from .lib_wrapper import ModelObjectWrapper
import numpy as np
from copy import deepcopy

__all__ = ['Imfit']

class Imfit(object):
    
    def __init__(self, model_descr, psf=None):
        if not isinstance(model_descr, ModelDescription):
            raise ValueError('model_descr must be a ModelDescription object.')
        self._modelDescr = model_descr
        self._psf = psf
        self._mask = None
        self._modelObject = None    
    

    def getModelDescription(self):
        return deepcopy(self._modelDescr)


    def getRawParameters(self):
        return np.array(self._modelObject.getRawParameters())


    def _setupModel(self):
        if self._modelObject is not None:
            # FIXME: Find a better way to free cython resources.
            self._modelObject.close()
        self._modelObject = ModelObjectWrapper(self._modelDescr)
        if self._psf is not None:
            self._modelObject.setPSF(np.asarray(self._psf))
            
    
    def fit(self, image, noise, mask=None, quiet=True):
        self._setupModel()
        if isinstance(image, np.ma.MaskedArray):
            if mask is None:
                mask = image.mask
            image = image.filled()
        if mask is None:
            mask = np.ones_like(image)
        else:
            self._mask = mask.astype('bool')

        if isinstance(noise, np.ma.MaskedArray):
            noise = noise.filled()

        image = image.astype('float64')
        noise = noise.astype('float64')
        mask = mask.astype('float64')
        
        self._modelObject.setData(image, noise, mask,
                                  n_combined=1, exp_time=1.0, gain=1.0, read_noise=0.0, original_sky=0.0)
        verbose = -1 if quiet else 1
        self._modelObject.fit(verbose=verbose)
        
        
    def getModelImage(self, shape=None):
        if self._modelObject is None:
            self._setupModel()
        if shape is not None:
            self._modelObject.setupModelImage(shape)

        image = self._modelObject.getModelImage()
        if self._mask is None:
            return image
        else:
            return np.ma.array(image, mask=self._mask)
        
        
    def __del__(self):
        if self._modelObject is not None:
            # FIXME: Find a better way to free cython resources.
            self._modelObject.close()
         
