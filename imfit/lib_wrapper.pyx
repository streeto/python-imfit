'''
Created on Sep 15, 2013

@author: andre
'''

from .imfit_lib cimport ModelObject, mp_par, mp_result
from .imfit_lib cimport AddFunctions, LevMarFit
from .imfit_lib cimport GetFunctionParameters, GetFunctionNames as GetFunctionNames_lib 
from .imfit_lib cimport AIC_corrected, BIC
from .imfit_lib cimport MASK_ZERO_IS_GOOD, WEIGHTS_ARE_SIGMAS
from .model import ModelDescription, FunctionDescription, ParameterDescription

cimport numpy as np
import numpy as np
from os import path
from copy import deepcopy

import cython
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdlib cimport calloc, free
from libc.string cimport memcpy

################################################################################

def function_names():
    cdef vector[string] func_names
    GetFunctionNames_lib(func_names)
    return [f for f in func_names]

################################################################################

def function_description(func_type, name=None):
    cdef int status
    cdef vector[string] parameters
    status = GetFunctionParameters(func_type, parameters)
    if status < 0:
        raise ValueError('Function %s not found.' % func_type)
    func_desc = FunctionDescription(func_type, name)
    for p in parameters:
        param_desc = ParameterDescription(p, value=0.0)
        func_desc.addParameter(param_desc)
    return func_desc

################################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double *alloc_copy_from_ndarray(np.ndarray[np.double_t, ndim=2, mode='c'] orig):
    cdef double *dest
    cdef int n_rows, n_cols, n_pix, imsize
    n_rows = orig.shape[0]
    n_cols = orig.shape[1]
    n_pix = n_rows * n_cols
    imsize = n_pix * sizeof(double)
    dest = <double *> calloc(n_pix, sizeof(double))
    memcpy(dest, &orig[0,0], imsize)
    return dest

################################################################################

cdef class ModelObjectWrapper(object):

    cdef ModelObject *_model 
    cdef mp_par *_paramInfo
    cdef double *_paramVect
    cdef bool _paramLimitsExist
    cdef int _nParams
    cdef int _nFreeParams
    cdef object _modelDescr
    cdef object _parameterList
    cdef int _nPixels, _nRows, _nCols
    cdef mp_result _fitResult
    cdef int _fitStatus
    
    cdef double *_imageData, *_errorData, *_maskData, *_psfData
    cdef bool _inputDataLoaded
    cdef bool _fitted
    cdef bool _freed
    

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, model_descr):
        self._paramLimitsExist = False
        self._paramInfo = NULL
        self._paramVect = NULL
        self._model = NULL

        self._imageData = NULL
        self._errorData = NULL
        self._maskData = NULL
        self._psfData = NULL
        self._inputDataLoaded = False
        self._fitted = False
        self._freed = False
        self._fitStatus = 0
        
        if not isinstance(model_descr, ModelDescription):
            raise ValueError('model_descr must be a ModelDescription object.')
        self._modelDescr = model_descr

        self._model = new ModelObject()
        if self._model == NULL:
            raise MemoryError('Could not allocate ModelObject.')
        # TODO: subsampling as an option.
        self._addFunctions(self._modelDescr, subsampling=True)
        self._paramSetup(self._modelDescr)
        
        
    def setMaxThreads(self, int nproc):
        self._model.SetMaxThreads(nproc)
        
        
    def _paramSetup(self, model_descr):
        self._parameterList = model_descr.parameterList()
        self._nParams = self._nFreeParams = self._model.GetNParams()
        if self._nParams != len(self._parameterList):
            raise Exception('Number of input parameters (%d) does not equal required number of parameters for specified functions (%d).' % (len(self._parameterList), self._nParams))

        self._paramInfo = <mp_par *> calloc(self._nParams, sizeof(mp_par))
        if self._paramInfo == NULL:
            raise MemoryError('Could not allocate parameter info.')
        self._paramVect = <double *> calloc(self._nParams, sizeof(double))
        if self._paramVect == NULL:
            raise MemoryError('Could not allocate parameter initial values.')
    
        # Fill parameter info and initial value.
        for i, param in enumerate(self._parameterList):
            if param.fixed:
                self._paramInfo[i].fixed = True
                self._nFreeParams -= 1
            elif param.limits is not None:
                self._paramInfo[i].limited[0] = True
                self._paramInfo[i].limited[1] = True
                self._paramInfo[i].limits[0] = param.limits[0]
                self._paramInfo[i].limits[1] = param.limits[1]
                self._paramLimitsExist = True
            self._paramVect[i] = param.value


    cdef _addFunctions(self, model_descr, bool subsampling):
        cdef int status = 0
        status = AddFunctions(self._model, model_descr.functionList(),
                              model_descr.functionSetIndices(), subsampling)
        if status < 0:
            raise RuntimeError('Failed to add the functions.')

    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def setPSF(self, np.ndarray[np.double_t, ndim=2, mode='c'] psf):
        cdef int n_rows_psf, n_cols_psf

        # Maybe this was called before.
        if self._psfData != NULL:
            free(self._psfData)
        self._psfData = alloc_copy_from_ndarray(psf)
        n_rows_psf = psf.shape[0]
        n_cols_psf = psf.shape[1]
        self._model.AddPSFVector(n_cols_psf * n_rows_psf, n_cols_psf, n_rows_psf, self._psfData)
        

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def setData(self,
                np.ndarray[np.double_t, ndim=2, mode='c'] image,
                np.ndarray[np.double_t, ndim=2, mode='c'] noise,
                np.ndarray[np.double_t, ndim=2, mode='c'] mask,
                int n_combined,
                double exp_time,
                double gain,
                double read_noise,
                double original_sky,
                int error_type=WEIGHTS_ARE_SIGMAS,
                int mask_format=MASK_ZERO_IS_GOOD):
        cdef int n_rows, n_cols, n_rows_err, n_cols_err
        cdef int n_pixels

        # Maybe this was called before.
        if self._inputDataLoaded:
            raise RuntimeError('Data already set.')
        if self._freed:
            raise RuntimeError('Objects already freed.')
            
        self._imageData = alloc_copy_from_ndarray(image)
        self._errorData = alloc_copy_from_ndarray(noise)
        self._maskData = alloc_copy_from_ndarray(mask)
        
        self._nRows = image.shape[0]
        self._nCols = image.shape[1]
        self._nPixels = self._nRows * self._nCols
            
        self._model.AddImageDataVector(self._imageData, self._nCols, self._nRows)
        self._model.AddImageCharacteristics(gain, read_noise, exp_time, n_combined, original_sky)
        self._model.AddErrorVector(self._nPixels, self._nCols, self._nRows, self._errorData, error_type)
        self._model.AddMaskVector(self._nPixels, self._nCols, self._nRows, self._maskData, mask_format)
        self._model.ApplyMask()
        self._inputDataLoaded = True


    def setupModelImage(self, shape):
        if self._inputDataLoaded:
            raise Exception('Input data already loaded.')
        self._nRows = shape[0]
        self._nCols = shape[1]
        self._nPixels = self._nRows * self._nCols
        self._model.SetupModelImage(self._nCols, self._nRows)
        self._model.CreateModelImage(self._paramVect)
        self._inputDataLoaded = True
        
        
    def fit(self, ftol=1e-8, verbose=-1):
        self._fitStatus = LevMarFit(self._nParams, self._nFreeParams, self._nPixels,
                                    self._paramVect, self._paramInfo,
                                    self._model, ftol, self._paramLimitsExist,
                                    self._fitResult, verbose)
        self._fitted = True
    
    
    def getModelDescription(self):
        model_descr = deepcopy(self._modelDescr)
        for i, p in enumerate(model_descr.parameterList()):
            p.value = self._paramVect[i]
        return model_descr
    
        
    def getRawParameters(self):
        vals = []
        for i in xrange(self._nParams):
            vals.append(self._paramVect[i])
        return vals
            
            
    def getModelImage(self):
        cdef double *model_image
        cdef np.ndarray[np.double_t, ndim=2, mode='c'] output_array
        cdef int imsize = self._nPixels * sizeof(double)

        model_image = self._model.GetModelImageVector()
        output_array = np.empty((self._nRows, self._nCols), dtype='float64')
        memcpy(&output_array[0,0], model_image, imsize)

        return output_array
        
        
    def getFitStatistic(self, mode='chi2'):
        cdef double chi2
        if self._fitted:
            chi2 = self._fitResult.bestnorm
        else:
            chi2 = self._model.GetFitStatistic(self._paramVect)
        if mode == 'chi2':
            return chi2

        cdef int n_valid_pix = self._model.GetNValidPixels()
        cdef int deg_free = n_valid_pix - self._nFreeParams
        if mode == 'reduced_chi2':
            return chi2 / deg_free
        if mode == 'AIC':
            return AIC_corrected(chi2, self._nFreeParams, n_valid_pix, 1)
        if mode == 'BIC':
            return BIC(chi2, self._nFreeParams, n_valid_pix, 1);


    @property
    def nPegged(self):
        if not self._fitted:
            raise Exception('Not fitted yet.')
        return self._fitResult.npegged
    
    
    @property
    def nIter(self):
        if not self._fitted:
            raise Exception('Not fitted yet.')
        return self._fitResult.niter
    
    
    @property
    def nFev(self):
        if not self._fitted:
            raise Exception('Not fitted yet.')
        return self._fitResult.nfev
    
    
    @property
    def nValidPixels(self):
        return self._model.GetNValidPixels()
    

    @property
    def fitConverged(self):
        # See Imfit/src/mpfit.cpp for magic numbers.
        return self._fitted and (self._fitStatus > 0) and (self._fitStatus < 5)
    
    
    @property
    def fitError(self):
        # See Imfit/src/mpfit.cpp for magic numbers.
        return self._fitted and self._fitStatus <= 0
    
    
    @property
    def fitTerminated(self):
        # See Imfit/src/mpfit.cpp for magic numbers.
        return self._fitted and self._fitStatus >= 5
    
    
    def close(self):
        if self._model != NULL:
            del self._model
        if self._paramInfo != NULL:
            free(self._paramInfo)
        if self._paramVect != NULL:
            free(self._paramVect)
            
        if self._imageData != NULL:
            free(self._imageData)
        if self._errorData != NULL:
            free(self._errorData)
        if self._maskData != NULL:
            free(self._maskData)
        if self._psfData != NULL:
            free(self._psfData)
        self._freed = True
        
        