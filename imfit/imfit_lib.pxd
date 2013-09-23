from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from 'imfit/param_struct.h':
    ctypedef struct mp_par:
        bint fixed          # 1 = fixed; 0 = free
        bint limited[2]     # 1 = low/upper limit; 0 = no limit
        double limits[2]    # lower/upper limit boundary value
 
        char *parname       # Name of parameter, or 0 for none
 
        double offset       # Offset to be added when printing/writing output
                            # (e.g., X0  or Y0 offset when an image subsection is
                            # being fitted)


cdef extern from 'imfit/model_object.h':
    cdef cppclass ModelObject:
        # WARNING: calling SetupModelImage and AddImageDataVector in the
        # same ModelObject instance (or any of them more than once) will
        # cause memory leak!
        void SetupModelImage(int nImageColumns, int nImageRows)
        void AddImageDataVector(double *image, int n_columns, int n_rows)
        
        void AddImageCharacteristics(double imageGain, double readoutNoise, double expTime, 
                                     int nCombinedImages, double originalSkyBackground)
        void AddErrorVector(int nDataValues, int nImageColumns, int nImageRows,
                            double *pixelVector, int inputType)
        void AddMaskVector(int nDataValues, int nImageColumns, int nImageRows,
                           double *pixelVector, int inputType)
        void ApplyMask()
        void AddPSFVector(int nPixels_psf, int nColumns_psf, int nRows_psf,
                          double *psfPixels)
        void SetMaxThreads(int maxThreadNumber)
        int GetNParams()
        int GetNValidPixels()
        void PrintDescription()
        double *GetModelImageVector()
        


cdef extern from 'imfit/add_functions.h':
    int AddFunctions(ModelObject *theModel, vector[string] &functionNameList,
                     vector[int] &functionSetIndices, bool subamplingFlag)
    int GetFunctionParameters(string &functionName, vector[string] &parameterNameList)
    void GetFunctionNames(vector[string] &functionNameList)


cdef extern from 'imfit/definitions.h':
    int WEIGHTS_ARE_SIGMAS    = 100 # "weight image" pixel value = sigma
    int WEIGHTS_ARE_VARIANCES = 110 # "weight image" pixel value = variance (sigma^2)
    int WEIGHTS_ARE_WEIGHTS   = 120 # "weight image" pixel value = weight
    int MASK_ZERO_IS_GOOD     =  10 # "standard" input mask format (good pixels = 0)
    int MASK_ZERO_IS_BAD      =  20 # alternate input mask format (good pixels = 1)


cdef extern from 'imfit/levmar_fit.h':
    int LevMarFit(int nParamsTot, int nFreeParams, int nDataVals, double *paramVector, 
                  mp_par *parameterLimits, ModelObject *theModel, double ftol, 
                  bool paramLimitsExist, int verbose)

