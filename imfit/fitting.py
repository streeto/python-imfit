'''
Created on Sep 20, 2013

@author: andre
'''
from .model import ModelDescription
from .lib import ModelObjectWrapper
import numpy as np
from copy import deepcopy

__all__ = ['Imfit']

class Imfit(object):
    '''
    A class for fitting models to images using Imfit by Peter Erwin.
    Can also be used to create images based on models.
    
    Due to some library limitations, this object can only fit the model
    specified in the constructor. If you want to fit several models,
    you need to create one instance of :class:`Imfit` for each model.
    On the other hand, one instance can be used to fit the same model
    to any number of images, or to fit and then create the model image.
    
    Parameters
    ----------
    model_descr : :class:`ModelDescription`
        Template model to be fitted, an instance of :class:`ModelDescription`.
        It will be the template model to every subsequent fitting in this instance.
        
    psf : 2-D array
        Point Spread Function image to be convolved to the images.
        Default: ``None`` (no convolution).
        
    quiet : bool, optional
        Suppress output, only error messages will be printed.
        Default: ``True``.
        
    nproc : int
        Number of processors to use when fitting.
        Default: ``-1`` (use all processors).
        
    See also
    --------
    parse_config_file, fit
    '''
    
    def __init__(self, model_descr, psf=None, quiet=True, nproc=-1, chunk_size=8):
        if not isinstance(model_descr, ModelDescription):
            raise ValueError('model_descr must be a ModelDescription object.')
        self._modelDescr = model_descr
        self._psf = psf
        self._mask = None
        self._modelObject = None
        self._nproc = nproc
        self._chunkSize = chunk_size
        self._debugLevel = -1 if quiet else 1
    

    def getModelDescription(self):
        '''
        Returns
        -------
        model : :class:`ModelDescription`
            A copy of the currently fitted model, or a copy of
            the template model if no fitting has taken place yet.
        '''
        if self._modelObject is not None:
            return self._modelObject.getModelDescription()
        else:
            # FIXME: get rid of deepcopy
            return deepcopy(self._modelDescr)


    def getRawParameters(self):
        '''
        Model parameters for debugging purposes.
        
        Returns
        -------
        raw_params : array of floats
            A flat array containing all the model parameter values.  
        '''
        return np.array(self._modelObject.getRawParameters())


    def _setupModel(self):
        if self._modelObject is not None:
            # FIXME: Find a better way to free cython resources.
            self._modelObject.close()
        self._modelObject = ModelObjectWrapper(self._modelDescr, self._debugLevel)
        if self._psf is not None:
            self._modelObject.setPSF(np.asarray(self._psf))
        if self._nproc > 0:
            self._modelObject.setMaxThreads(self._nproc)
        if self._chunkSize > 0:
            self._modelObject.setChunkSize(self._chunkSize)
            
    
    def fit(self, image, noise, mask=None, mode='LM'):
        '''
        Fit the model to ``image``, using the inverse of ``noise`` as weight,
        optionally masking some pixels.
        
        Parameters
        ----------
        image : 2-D array
            Image to be fitted. Can be a masked array.
        
        noise : 2-D array
            Noise image, same shape as ``image``.
        
        mask : 2-D array, optional
            Array containing the masked pixels (``True`` is bad).
            If not set and ``image`` is a masked array, ``image.mask`` is used.
            
        mode : string
            One of the following algorithms:
                * ``'LM'`` : Levenberg-Marquardt least squares.
                * ``'DE'`` : Differential Evolution.
                * ``'NM'`` : Nelder-Mead Simplex.
            
        Examples
        --------
        TODO: Examples of fit().
        
        '''
        if mode not in ['LM', 'DE', 'NM']:
            raise Exception('Invalid mode: %s' % mode)
        self._setupModel()
        if isinstance(image, np.ma.MaskedArray):
            if mask is None:
                mask = image.mask
            image = image.filled(fill_value=0.0)
        if mask is None:
            mask = np.ones_like(image)
        else:
            self._mask = mask.astype('bool')

        if (image.shape != noise.shape) or (image.shape != mask.shape):
            raise Exception('Image, noise and mask shapes do not match.')

        if isinstance(noise, np.ma.MaskedArray):
            noise = noise.filled(fill_value=noise.max())

        image = image.astype('float64')
        noise = noise.astype('float64')
        mask = mask.astype('float64')
        
        self._modelObject.setData(image, noise, mask,
                                  n_combined=1, exp_time=1.0, gain=1.0, read_noise=0.0, original_sky=0.0)
        self._modelObject.fit(verbose=self._debugLevel, mode=mode)
        
    
    @property
    def fitConverged(self):
        return self._modelObject.fitConverged

    
    @property
    def fitError(self):
        return self._modelObject.fitError
    
    
    @property
    def fitTerminated(self):
        return self._modelObject.fitTerminated
    
    
    @property
    def nIter(self):
        return self._modelObject.nIter

    @property
    def nPegged(self):
        return self._modelObject.nPegged

    
    @property
    def nValidPixels(self):
        return self._modelObject.nValidPixels
    
    
    @property
    def chi2(self):
        '''
        The :math:`\\chi^2` statistic of the fit.
        '''
        return self._modelObject.getFitStatistic(mode='chi2')
    
    
    @property
    def AIC(self):
        '''
        Bias-corrected Akaike Information Criterion for the fit.
        '''
        return self._modelObject.getFitStatistic(mode='AIC')
    
    
    @property
    def BIC(self):
        '''
        Bayesian Information Criterion for the fit.
        '''
        return self._modelObject.getFitStatistic(mode='BIC')
    
    
    def getModelImage(self, shape=None):
        '''
        Computes an image from the currently fitted model.
        If not fitted, use the template model.
        
        Parameters
        ----------
        shape : tuple
            Shape of the image in (Y, X) format.
            
        Returns
        -------
        image : 2-D array
            Image computed from the current model.
        '''
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
         
