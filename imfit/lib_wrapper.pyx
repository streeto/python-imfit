'''
Created on Sep 15, 2013

@author: andre
'''

from .imfit_lib cimport ModelObject, mp_par, configOptions
from .imfit_lib cimport ReadConfigFile, AddFunctions, LevMarFit, PrintResults, mp_result
from .imfit_lib cimport GetFunctionParameters, GetFunctionNames as GetFunctionNames_lib 
from .imfit_lib cimport MASK_ZERO_IS_GOOD, WEIGHTS_ARE_SIGMAS
from .model import ModelDescription, FunctionDescription, ParameterDescription

cimport numpy as np
import numpy as np
from os import path

import cython
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdlib cimport calloc, free
from libc.string cimport memcpy


def getFunctionNames():
    cdef vector[string] func_names
    GetFunctionNames_lib(func_names)
    return [f for f in func_names]


def getFunctionDescription(func_name):
    cdef int status
    cdef string func_name_str = func_name
    cdef vector[string] parameters
    status = GetFunctionParameters(func_name_str, parameters)
    if status < 0:
        raise ValueError('Function %s not found.' % func_name)
    func_desc = FunctionDescription(func_name)
    for p in parameters:
        param_desc = ParameterDescription(p, value=0.0)
        func_desc.addParameter(param_desc)
    return func_desc


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



cdef class ModelObjectWrapper(object):

    cdef ModelObject *_model 
    cdef mp_par *_paramInfo
    cdef double *_paramVect
    cdef bool _paramLimitsExist
    cdef int _nParams
    cdef int _nFreeParams
    cdef int _nPixels
    
    cdef double *_imageData, *_errorData, *_maskData, *_psfData

    
    cdef object _modelDescr
    cdef object _parameterList

    def __init__(self, model_description):
        self._paramLimitsExist = False
        self._paramInfo = NULL
        self._paramVect = NULL
        self._model = NULL

        self._imageData = NULL
        self._errorData = NULL
        self._maskData = NULL
        self._psfData = NULL
        
        if not isinstance(model_description, ModelDescription):
            raise ValueError('model_description must be a ModelDescription object.')
        self._modelDescr = model_description
        self._parameterList = self._modelDescr.parameterList()

        cdef int status = 0
        
        self._model = new ModelObject()
        if self._model == NULL:
            raise MemoryError('Could not allocate ModelObject.')
        
        # TODO: subsampling as an option. 
        status = AddFunctions(self._model, self._modelDescr.functionList(),
                              self._modelDescr.functionSetIndices(), True)
        if status < 0:
            raise Exception('Failed to add the functions.')
        
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
        if self._imageData != NULL:
            raise ValueError('Data already set.')
            
        self._imageData = alloc_copy_from_ndarray(image)
        self._errorData = alloc_copy_from_ndarray(noise)
        self._maskData = alloc_copy_from_ndarray(mask)
        
        n_rows = image.shape[0]
        n_cols = image.shape[1]
        self._nPixels = n_rows * n_cols
        print 'nRows = %d, nColumns = %d' % (n_rows, n_cols)
            
        self._model.AddImageDataVector(self._imageData, n_cols, n_rows)
        self._model.AddImageCharacteristics(gain, read_noise, exp_time, n_combined, original_sky)
        self._model.AddErrorVector(self._nPixels, n_cols, n_rows, self._errorData, error_type)
        self._model.AddMaskVector(self._nPixels, n_cols, n_rows, self._maskData, mask_format)
        self._model.ApplyMask()


    @cython.boundscheck(False)
    @cython.wraparound(False)
    def fit(self,
            ftol=1e-8,
            verbose=-1):
                  
        print 'Calling Levenberg-Marquardt solver ...'
        print self._nParams, self._nFreeParams, self._nPixels, ftol, self._paramLimitsExist
        status = LevMarFit(self._nParams, self._nFreeParams, self._nPixels, self._paramVect, self._paramInfo,
                           self._model, ftol, self._paramLimitsExist, verbose)    
        self._updateParamValues() 
    
    
    def _updateParamValues(self):
        for i, p in enumerate(self._parameterList):
            p.value = self._paramVect[i]
    
        
    def __del__(self):
        if self._model != NULL:
            del self.theModel
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
    
    
    def __str__(self):
        return str(self._modelDescr)


@cython.boundscheck(False)
@cython.wraparound(False)
def fit_wrapper(np.ndarray[np.double_t, ndim=2, mode='c'] image,
                np.ndarray[np.double_t, ndim=2, mode='c'] mask,
                np.ndarray[np.double_t, ndim=2, mode='c'] noise,
                np.ndarray[np.double_t, ndim=2, mode='c'] psf,
                model,
                int nCombined,
                double expTime,
                double gain,
                double readNoise,
                double originalSky,
                double ftol,
                int verbose=-1,
                int errorType=WEIGHTS_ARE_SIGMAS,
                int maskFormat=MASK_ZERO_IS_GOOD,
                copy=True):
    '''
    Wrapper for the whole fitting, including input as a config file.
    '''
        
    cdef int nPixels_tot, nColumns, nRows
    cdef int nPixels_psf, nRows_psf, nColumns_psf
    cdef int nErrColumns, nErrRows, nMaskColumns, nMaskRows
    cdef int nDegFreedom
    cdef int nParamsTot, nFreeParams
    cdef double *allPixels
    cdef double *psfPixels
    cdef double *allErrorPixels
    cdef double *allMaskPixels
    cdef double *paramsVect
    cdef ModelObject *theModel
    
    cdef bool paramLimitsExist = False
    cdef mp_par *parameterInfo
    cdef int status
    cdef int i
    
    theModel = new ModelObject()
    theModel.SetMaxThreads(2)     
    
    status = AddFunctions(theModel, model.functionList(), model.functionSetIndices(), True)
    if status < 0:
        print '*** WARNING: Failure in AddFunctions!\n\n'
        return

    parameterList = model.parameterList()
    
    # Set up parameter vector(s), now that we know total # parameters
    nParamsTot = nFreeParams = theModel.GetNParams();
    if verbose >= 0:
        print '%d total parameters' % nParamsTot
    if nParamsTot != len(parameterList):
        print '*** WARNING: number of input parameters (%d) does not equal required number of parameters for specified functions (%d)!' % (len(parameterList), nParamsTot)
        return

    # The routines change the values of the pixels, not very nice.
    # By default we make a copy of the arrays.
    if copy:
        image = image.copy()
        mask = mask.copy()
        noise = noise.copy()
        psf = psf.copy()
    
    nRows = image.shape[0]
    nColumns = image.shape[1]
    nPixels_tot = nColumns * nRows
    allPixels = &image[0,0]

    if verbose >= 0:
        print 'nRows = %d, nColumns = %d' % (nRows, nColumns)
    
    nMaskRows = mask.shape[0]
    nMaskColumns = mask.shape[1]
    allMaskPixels = &mask[0,0]

    nErrRows = noise.shape[0]
    nErrColumns = noise.shape[1]
    allErrorPixels = &noise[0,0]

    nRows_psf = psf.shape[0]
    nColumns_psf = psf.shape[1]
    nPixels_psf = nColumns_psf * nRows_psf
    psfPixels = &psf[0,0]
    
    theModel.AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels)
    theModel.AddImageDataVector(allPixels, nColumns, nRows)
    theModel.AddImageCharacteristics(gain, readNoise, expTime, nCombined,
                                     originalSky)
    theModel.AddErrorVector(nPixels_tot, nColumns, nRows, allErrorPixels,
                            errorType)
    theModel.AddMaskVector(nPixels_tot, nColumns, nRows, allMaskPixels,
                           maskFormat)
    theModel.ApplyMask()

    if verbose >= 0:
        theModel.PrintDescription()

    '''
    START OF MINIMIZATION-ROUTINE-RELATED CODE
    Parameter limits and other info:
    First we create a C-style array of mp_par structures, containing parameter constraints
    (if any) *and* any other useful info (like X0,Y0 offset values).  This will be used
    by DiffEvolnFit (if called) and by PrintResults.  We also decrement nFreeParams for
    each *fixed* parameter.
    Then we point the mp_par-array variable mpfitParameterConstraints to this array *if*
    there are actual parameter constraints; if not, mpfitParameterConstraints is set = NULL,
    since that's what mpfit() expects when there are no constraints.
    '''
    if verbose >= 0:
        print 'Setting up parameter information array ...'
    parameterInfo = <mp_par *> calloc(nParamsTot, sizeof(mp_par))
    paramsVect = <double *> calloc(nParamsTot, sizeof(double))

    for i, param in enumerate(parameterList):
        if param.fixed:
            parameterInfo[i].fixed = True
            nFreeParams -= 1
        elif param.limits is not None:
            parameterInfo[i].limited[0] = True
            parameterInfo[i].limited[1] = True
            parameterInfo[i].limits[0] = param.limits[0]
            parameterInfo[i].limits[1] = param.limits[1]
            paramLimitsExist = True
        # Copy initial parameter values into C array
        paramsVect[i] = param.value
    
    if verbose >= 0:
        nDegFreedom = theModel.GetNValidPixels() - nFreeParams
        print '%d free parameters (%d degrees of freedom)' % (nFreeParams, nDegFreedom)
    
    # Check the model.
#     PrintResults(paramsVect, <double*>0, <mp_result*>0, theModel, nFreeParams, parameterInfo, 1)
    
    # DO THE FIT!
    print 'Calling Levenberg-Marquardt solver ...'
    status = LevMarFit(nParamsTot, nFreeParams, nPixels_tot, paramsVect, parameterInfo,
                       theModel, ftol, paramLimitsExist, verbose)
 
    # TODO: return a ModelDescription with the fitted result.

    for i, p in enumerate(parameterList):
        p.value = paramsVect[i]
        
#     values = []
#     for i in xrange(len(parameterList)):
#         values.append(paramsVect[i])
# 
#     model.updateParamValues(values)
    
    free(paramsVect)
    free(parameterInfo)
    del theModel
