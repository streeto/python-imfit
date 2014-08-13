'''
Created on Sep 15, 2013

@author: andre
'''

from .imfit_lib cimport ModelObject, mp_par, mp_result
from .imfit_lib cimport AddFunctions, LevMarFit, DiffEvolnFit, NMSimplexFit
from .imfit_lib cimport GetFunctionParameters, GetFunctionNames as GetFunctionNames_lib 
from .imfit_lib cimport AIC_corrected, BIC
from .imfit_lib cimport MASK_ZERO_IS_GOOD, MASK_ZERO_IS_BAD
from .imfit_lib cimport WEIGHTS_ARE_SIGMAS, WEIGHTS_ARE_VARIANCES, WEIGHTS_ARE_WEIGHTS
from .imfit_lib cimport Convolver

from ..model import ModelDescription, FunctionDescription, ParameterDescription

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


__all__ = ['function_types', 'function_description', 'convolve_image', 'ModelObjectWrapper']

################################################################################

def function_types():
    '''
    List the available model function types.
    
    Returns
    -------
    func_types : list
        A list containing the function types (string).
    '''
    cdef vector[string] func_names
    GetFunctionNames_lib(func_names)
    return [f for f in func_names]

################################################################################

def function_description(func_type, name=None):
    '''
    Get the description of a function type.
    
    Parameters
    ----------
    func_type : string
        The function type.
        
    name : string, optional
        Custom name to the function. Example: 'disk', 'bulge'.
        Default: None.
        
    Returns
    -------
    func_desc : :class:`FunctionDescription`
        Instance of :class:`FunctionDescription`.
        
    '''
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

def convolve_image(np.ndarray[np.double_t, ndim=2, mode='c'] image not None,
                   np.ndarray[np.double_t, ndim=2, mode='c'] psf not None,
                   int nproc=0,
                   verbose=False,
                   do_fftw_measure=False):
    '''
    Convolve an image with a given PSF.
    
    Parameters
    ----------
    image : array
        Image to be convolved.
        
    psf : array
        PSF to apply.
        
    nproc : int, optional
        Number of threads to use. If ```nproc <= 0``, use all available cores.
        Default: ``0``, use all cores.
    
    verbose : bool, optional
        Print diagnostic messages.
        Default: ``False`` ,be quiet.
        
    do_fftw_measure : bool, optional
        Tells FFTW to find an optimized plan by actually computing several
        FFTs and measuring their execution time. This can be slow.
        Default: ``False``, use a faster estimation of plan, which
        can actually be slower to execute.
        
    Returns
    -------
    convolved_image : array
        An array of same shape as ``image`` containing
        the convolved image.
        
    '''
    cdef double *psf_data = alloc_copy_from_ndarray(psf)
    cdef Convolver *convolver = new Convolver()
    convolver.SetupPSF(psf_data, psf.shape[1], psf.shape[0])
    if nproc >= 0:
        convolver.SetMaxThreads(nproc)
        
    convolver.SetupImage(image.shape[1], image.shape[0])
    cdef int debug_level = 1 if verbose else -1
    convolver.DoFullSetup(debug_level, do_fftw_measure)
    cdef double *image_data = alloc_copy_from_ndarray(image)
    convolver.ConvolveImage(image_data)
    
    shape = (image.shape[1], image.shape[0])
    cdef np.ndarray[np.double_t, ndim=2, mode='c'] convolved_image
    convolved_image = np.empty(shape, dtype='float64')
    memcpy(&convolved_image[0,0], image_data, image.nbytes)
    free(image_data)
    del convolver
    free(psf_data)
    return convolved_image

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
    
    cdef double *_imageData
    cdef double *_errorData
    cdef double *_maskData
    cdef double *_psfData
    cdef bool _inputDataLoaded
    cdef bool _fitted
    cdef object _fitMode
    cdef bool _freed
    

    def __init__(self, object model_descr, int debug_level=0, int verbose_level=-1, bool subsampling=True):
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
        self._fitMode = None
        self._freed = False
        self._fitStatus = 0
        
        if not isinstance(model_descr, ModelDescription):
            raise ValueError('model_descr must be a ModelDescription object.')
        self._modelDescr = model_descr

        self._model = new ModelObject()
        self._model.SetDebugLevel(debug_level)
        self._model.SetVerboseLevel(verbose_level)
        if self._model == NULL:
            raise MemoryError('Could not allocate ModelObject.')
        self._addFunctions(self._modelDescr, subsampling=subsampling, verbose=debug_level>0)
        self._paramSetup(self._modelDescr)
        
        
    def setMaxThreads(self, int nproc):
        self._model.SetMaxThreads(nproc)
        
        
    def setChunkSize(self, int chunk_size):
        self._model.SetOMPChunkSize(chunk_size)
        
        
    def _paramSetup(self, object model_descr):
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


    cdef _addFunctions(self, object model_descr, bool subsampling, bool verbose=False):
        cdef int status = 0
        status = AddFunctions(self._model, model_descr.functionList(),
                              model_descr.functionSetIndices(), subsampling, verbose)
        if status < 0:
            raise RuntimeError('Failed to add the functions.')

    
    def setPSF(self, np.ndarray[np.double_t, ndim=2, mode='c'] psf):
        cdef int n_rows_psf, n_cols_psf

        # Maybe this was called before.
        if self._psfData != NULL:
            free(self._psfData)
        self._psfData = alloc_copy_from_ndarray(psf)
        n_rows_psf = psf.shape[0]
        n_cols_psf = psf.shape[1]
        self._model.AddPSFVector(n_cols_psf * n_rows_psf, n_cols_psf, n_rows_psf, self._psfData)
        

    def loadData(self,
                 np.ndarray[np.double_t, ndim=2, mode='c'] image not None,
                 np.ndarray[np.double_t, ndim=2, mode='c'] error,
                 np.ndarray[np.double_t, ndim=2, mode='c'] mask,
                 **kwargs):
        
        cdef int n_rows, n_cols, n_rows_err, n_cols_err
        cdef int n_pixels

        # Maybe this was called before.
        if self._inputDataLoaded:
            raise RuntimeError('Data already loaded.')
        if self._freed:
            raise RuntimeError('Objects already freed.')
            
        # kwargs
        cdef int n_combined
        cdef double exp_time
        cdef double gain
        cdef double read_noise
        cdef double original_sky
        cdef int error_type
        cdef int mask_format
        cdef bool use_cash_statistics
        cdef bool use_model_for_errors
        
        if 'n_combined' in kwargs:
            n_combined = kwargs['n_combined']
        else:
            n_combined = 1
            
        if 'exp_time' in kwargs:
            exp_time = kwargs['exp_time']
        else:
            exp_time = 1.0

        if 'gain' in kwargs:
            gain = kwargs['gain']
        else:
            gain = 1.0

        if 'read_noise' in kwargs:
            read_noise = kwargs['read_noise']
        else:
            read_noise = 0.0

        if 'original_sky' in kwargs:
            original_sky = kwargs['original_sky']
        else:
            original_sky = 0.0

        if 'error_type' in kwargs:
            if kwargs['error_type'] == 'sigma':
                error_type = WEIGHTS_ARE_SIGMAS
            elif kwargs['error_type'] == 'variance':
                error_type = WEIGHTS_ARE_VARIANCES
            elif kwargs['error_type'] == 'weight':
                error_type = WEIGHTS_ARE_WEIGHTS
            else:
                raise Exception('Unknown error type: %s' % kwargs['error_type'])
        else:
            error_type = WEIGHTS_ARE_SIGMAS

        if 'mask_format' in kwargs:
            if kwargs['mask_format'] == 'zero_is_good':
                mask_format = MASK_ZERO_IS_GOOD
            elif kwargs['mask_format'] == 'zero_is_bad':
                mask_format = MASK_ZERO_IS_BAD
            else:
                raise Exception('Unknown mask format: %s' % kwargs['mask_format'])
        else:
            mask_format = MASK_ZERO_IS_GOOD
            
        if 'use_cash_statistics' in kwargs:
            use_cash_statistics = kwargs['use_cash_statistics']
        else:
            use_cash_statistics = False            
            
        if 'use_model_for_errors' in kwargs:
            use_model_for_errors = kwargs['use_model_for_errors']
        else:
            use_model_for_errors = False            
            
        self._imageData = alloc_copy_from_ndarray(image)
        self._nRows = image.shape[0]
        self._nCols = image.shape[1]
        self._nPixels = self._nRows * self._nCols
            
        self._model.AddImageDataVector(self._imageData, self._nCols, self._nRows)
        self._model.AddImageCharacteristics(gain, read_noise, exp_time, n_combined, original_sky)
        
        if use_cash_statistics:
            self._model.UseCashStatistic()
        else:
            if error is not None:
                self._errorData = alloc_copy_from_ndarray(error)
                self._model.AddErrorVector(self._nPixels, self._nCols, self._nRows, self._errorData, error_type)
            elif use_model_for_errors:
                self._model.UseModelErrors()
        
        if mask is not None:
            self._maskData = alloc_copy_from_ndarray(mask)
            success = self._model.AddMaskVector(self._nPixels, self._nCols, self._nRows, self._maskData, mask_format)
            if success != 0:
                raise Exception('Error adding mask vector, unknown mask format.')

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
        
        
    def _testCreateModelImage(self, int count=1):
        for _ from 0 <= _ < count:
            self._model.CreateModelImage(self._paramVect)
        
        
    def fit(self, double ftol=1e-8, int verbose=-1, mode='LM'):
        status = self._model.FinalSetupForFitting()
        if status < 0:
            raise Exception('Failure in ModelObject::FinalSetupForFitting().')
        if mode == 'LM':
            if self._model.UsingCashStatistic():
                raise Exception('Cannot use Cash statistic with L-M solver.')
            self._fitStatus = LevMarFit(self._nParams, self._nFreeParams, self._nPixels,
                                        self._paramVect, self._paramInfo,
                                        self._model, ftol, self._paramLimitsExist,
                                        self._fitResult, verbose)
        elif mode == 'DE':
            self._fitStatus = DiffEvolnFit(self._nParams, self._paramVect, self._paramInfo,
                                           self._model, ftol, verbose)
        elif mode == 'NM':
            self._fitStatus = NMSimplexFit(self._nParams, self._paramVect, self._paramInfo,
                                           self._model, ftol, verbose)
        else:
            raise Exception('Invalid fit mode: %s' % mode)

        self._fitMode = mode
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
        if model_image is NULL:
            raise Exception('Error: model image has not yet been computed.')
        output_array = np.empty((self._nRows, self._nCols), dtype='float64')
        memcpy(&output_array[0,0], model_image, imsize)

        return output_array
        
        
    def getFitStatistic(self, mode='none'):
        cdef double fitstat
        if self.fittedLM:
            fitstat = self._fitResult.bestnorm
        else:
            fitstat = self._model.GetFitStatistic(self._paramVect)
        cdef int n_valid_pix = self._model.GetNValidPixels()
        cdef int deg_free = n_valid_pix - self._nFreeParams

        if mode == 'none':
            return fitstat
        elif mode == 'reduced':
            return fitstat / deg_free
        elif mode == 'AIC':
            return AIC_corrected(fitstat, self._nFreeParams, n_valid_pix, 1)
        elif mode == 'BIC':
            return BIC(fitstat, self._nFreeParams, n_valid_pix, 1);
        else:
            raise Exception('Unknown statistic mode: %s' % mode)


    @property
    def fittedLM(self):
        return self._fitted and (self._fitMode == 'LM')


    @property
    def nPegged(self):
        if self.fittedLM:
            return self._fitResult.npegged
        else:
            return -1
    
    
    @property
    def nIter(self):
        if self.fittedLM:
            return self._fitResult.niter
        else:
            return -1
    
    
    @property
    def nFev(self):
        if self.fittedLM:
            return self._fitResult.nfev
        else:
            return -1
    
    
    @property
    def nValidPixels(self):
        return self._model.GetNValidPixels()
    

    @property
    def validPixelFraction(self):
        return self._model.GetNValidPixels() / self._nPixels
    

    @property
    def fitConverged(self):
        if not self._fitted:
            raise Exception('Not fitted yet.')
        return (self._fitStatus > 0) and (self._fitStatus < 5)
    
    
    @property
    def fitError(self):
        if not self._fitted:
            raise Exception('Not fitted yet.')
        return self._fitStatus <= 0
    
    
    @property
    def fitTerminated(self):
        if not self._fitted:
            raise Exception('Not fitted yet.')
        # See Imfit/src/mpfit.cpp for magic numbers.
        return self._fitStatus >= 5
    
    
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
        
        