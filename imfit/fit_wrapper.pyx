'''
Created on Sep 15, 2013

@author: andre
'''

from imfit cimport ModelObject, mp_par, configOptions
from imfit cimport ReadConfigFile, AddFunctions, LevMarFit
from imfit cimport MASK_ZERO_IS_GOOD, WEIGHTS_ARE_SIGMAS

cimport numpy as np
import numpy as np
from os import path

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from stdlib cimport calloc, free 

def fit(np.ndarray[np.double_t, ndim=2] image,
        np.ndarray[np.double_t, ndim=2] mask,
        np.ndarray[np.double_t, ndim=2] noise,
        np.ndarray[np.double_t, ndim=2] psf,
        str configFileName,
        int nCombined,
        double expTime,
        double gain,
        double readNoise,
        double originalSky,
        double ftol,
        int verbose=1,
        int errorType=WEIGHTS_ARE_SIGMAS,
        int maskFormat=MASK_ZERO_IS_GOOD):
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
    
    cdef vector[string] functionList
    cdef vector[double] parameterList
    cdef vector[mp_par] paramLimits
    cdef vector[int] functionSetIndices
    cdef bool paramLimitsExist = False
    cdef mp_par *parameterInfo
    cdef int status
    cdef configOptions userConfigOptions
    cdef int i
     
    if not path.exists(configFileName):
        print '*** WARNING: Unable to find configuration file \"%s\"!' % configFileName
        return
    status = ReadConfigFile(configFileName, True, functionList, parameterList,
                            paramLimits, functionSetIndices, paramLimitsExist, userConfigOptions)
 
    if status != 0:
        print '*** WARNING: Failure reading configuration file!'
        return
     
    image = np.ascontiguousarray(image).copy()
    nColumns = image.shape[0]
    nRows = image.shape[1]
    nPixels_tot = nColumns * nRows
    allPixels = &image[0,0]
     
    mask = np.ascontiguousarray(mask).copy()
    nMaskColumns = mask.shape[0]
    nMaskRows = mask.shape[1]
    allMaskPixels = &mask[0,0]
 
    noise = np.ascontiguousarray(noise).copy()
    nErrColumns = noise.shape[0]
    nErrRows = noise.shape[1]
    allErrorPixels = &noise[0,0]
 
    psf = np.ascontiguousarray(psf).copy()
    nColumns_psf = psf.shape[0]
    nRows_psf = psf.shape[1]
    nPixels_psf = nColumns_psf * nRows_psf
    psfPixels = &psf[0,0]
    
    theModel = new ModelObject()
    theModel.SetMaxThreads(2)
     
    status = AddFunctions(theModel, functionList, functionSetIndices, False)
    if status < 0:
        print '*** WARNING: Failure in AddFunctions!\n\n'
        return
#     
    # Set up parameter vector(s), now that we know total # parameters
    nParamsTot = nFreeParams = theModel.GetNParams();
    print '%d total parameters' % nParamsTot
    if nParamsTot != parameterList.size():
        print '*** WARNING: number of input parameters (%d) does not equal required number of parameters for specified functions (%d)!' % (parameterList.size(), nParamsTot)
        return
 
    theModel.AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels)
    theModel.AddImageDataVector(allPixels, nColumns, nRows)
    theModel.AddImageCharacteristics(gain, readNoise, expTime, nCombined,
                                     originalSky)
    theModel.PrintDescription()
    theModel.AddErrorVector(nPixels_tot, nColumns, nRows, allErrorPixels,
                            errorType)
    theModel.AddMaskVector(nPixels_tot, nColumns, nRows, allMaskPixels,
                           maskFormat)
    theModel.ApplyMask()
 
    '''
    START OF MINIMIZATION-ROUTINE-RELATED CODE */
    Parameter limits and other info:
    First we create a C-style array of mp_par structures, containing parameter constraints
    (if any) *and* any other useful info (like X0,Y0 offset values).  This will be used
    by DiffEvolnFit (if called) and by PrintResults.  We also decrement nFreeParams for
    each *fixed* parameter.
    Then we point the mp_par-array variable mpfitParameterConstraints to this array *if*
    there are actual parameter constraints; if not, mpfitParameterConstraints is set = NULL,
    since that's what mpfit() expects when there are no constraints.
    '''
    print 'Setting up parameter information array ...'
    parameterInfo = <mp_par *> calloc(nParamsTot, sizeof(mp_par))
    for i in xrange(nParamsTot):
        parameterInfo[i].fixed = paramLimits[i].fixed
        if parameterInfo[i].fixed == 1:
            nFreeParams -= 1
        parameterInfo[i].limited[0] = paramLimits[i].limited[0]
        parameterInfo[i].limited[1] = paramLimits[i].limited[1]
        parameterInfo[i].limits[0] = paramLimits[i].limits[0]
        parameterInfo[i].limits[1] = paramLimits[i].limits[1]
 
    nDegFreedom = theModel.GetNValidPixels() - nFreeParams
    print '%d free parameters (%d degrees of freedom)' % (nFreeParams, nDegFreedom)
   
   
    # Copy initial parameter values into C array, correcting for X0,Y0 offsets
    paramsVect = <double *> calloc(nParamsTot, sizeof(double))
    for i in xrange(nParamsTot):
        paramsVect[i] = parameterList[i]
   
    # DO THE FIT!
    print 'Calling Levenberg-Marquardt solver ...'
    status = LevMarFit(nParamsTot, nFreeParams, nPixels_tot, paramsVect, parameterInfo,
                       theModel, ftol, paramLimitsExist, verbose)
 
#     free(paramsVect)
#     free(parameterInfo)
   
    print 'Done!'
