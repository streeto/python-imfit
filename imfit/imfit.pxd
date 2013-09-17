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


cdef extern from 'imfit/config_file_parser.h':
    ctypedef struct configOptions:
        vector[string] optionNames
        vector[string] optionValues
        int  nOptions

    int ReadConfigFile(string& configFileName, bool mode2D, vector[string]& functionList,
                       vector[double]& parameterList, vector[mp_par]& parameterLimits,
                       vector[int]& setStartFunctionNumber, bool& parameterLimitsFound,
                       configOptions& configFileOptions)



cdef extern from 'imfit/model_object.h':
    cdef cppclass ModelObject:
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
        


cdef extern from 'imfit/add_functions.h':
    int AddFunctions(ModelObject *theModel, vector[string] &functionNameList,
                     vector[int] &functionSetIndices, bool subamplingFlag)


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


# {
#   public:
#     // Constructors:
#     ModelObject( );
# 
#     void SetDebugLevel( int debuggingLevel );
#     
#     void SetMaxThreads( int maxThreadNumber );
#     
#     // common, not specialized
#     // Adds a new FunctionObject pointer to the internal vector
#     void AddFunction( FunctionObject *newFunctionObj_ptr );
#     
#     // common, but Specialized by ModelObject1D
#     virtual void DefineFunctionSets( vector<int>& functionStartIndices );
#     
#     // 1D only, but needs to be part of base interface
#     virtual void AddDataVectors( int nDataValues, double *xValVector, 
#                             double *yValVector, bool magnitudeData ) { nDataVals = nDataValues; };
# 
#     // Probably 1D only, but might be usable by 2D version later...
#     virtual void SetZeroPoint( double zeroPointValue );
# 
#  
# //     void SetGain( double gainValue );
# // 
# //     void SetSkyBackground( double originalSkyValue );
# 
#     // 2D only
#     void AddImageDataVector( double *pixelVector, int nImageColumns, int nImageRows );
# 
#     // 2D only
#     void AddImageCharacteristics( double imageGain, double readoutNoise, double expTime, 
#                                 int nCombinedImages, double originalSkyBackground );
#     
#     // 2D only
#     void SetupModelImage( int nImageColumns, int nImageRows );
#     
#     // 2D only
#     virtual void AddErrorVector( int nDataValues, int nImageColumns, int nImageRows,
#                          double *pixelVector, int inputType );
# 
#     // 1D only
#     virtual void AddErrorVector1D( int nDataValues, double *pixelVector, int inputType ) { ; };
# 
#     // 1D only
#     virtual void AddMaskVector1D( int nDataValues, double *inputVector, int inputType ) { ; };
#     
#     // 2D only
#     virtual void GenerateErrorVector( );
# 
#     // 2D only
#     virtual void AddMaskVector( int nDataValues, int nImageColumns, int nImageRows,
#                          double *pixelVector, int inputType );
# 
#     // 2D only
#     void AddPSFVector( int nPixels_psf, int nColumns_psf, int nRows_psf,
#                          double *psfPixels );
# 
#     // 1D only
#     virtual void AddPSFVector1D( int nPixels_psf, double *xValVector, double *yValVector ) { ; };
#     
#     // 2D only [1D maybe needs something similar, but with diff. interface]
#     virtual void ApplyMask( );
# 
#     // common, but Specialized by ModelObject1D
#     virtual void CreateModelImage( double params[] );
#     
#     // 2D only
#     void UpdateWeightVector(  );
# 
#     // Specialized by ModelObject1D
#     virtual void ComputeDeviates( double yResults[], double params[] );
# 
#      // common, not specialized
#     virtual void UseModelErrors( );
# 
#      // common, not specialized
#     virtual void UseCashStatistic( );
#  
#      // common, not specialized
#     virtual bool UsingCashStatistic( );
#  
#     // common, not specialized
#     virtual double GetFitStatistic( double params[] );
#     
#     // common, not specialized
#     virtual double ChiSquared( double params[] );
#     
#     // common, not specialized
#     virtual double CashStatistic( double params[] );
#     
#     // common, but Specialized by ModelObject1D
#     virtual void PrintDescription( );
# 
#     // common, not specialized
#     void GetFunctionNames( vector<string>& functionNames );
# 
#     // common, but Specialized by ModelObject1D
#     virtual void PrintModelParams( FILE *output_ptr, double params[], mp_par *parameterInfo,
#                                                                         double errs[] );
# 
#     // 2D only; NOT USED ANYWHERE!
#     void PrintImage( double *pixelVector, int nColumns, int nRows );
# 
#         // 2D only
#     void PrintInputImage( );
# 
#         // 2D only
#     void PrintModelImage( );
# 
#     // 2D only; NOT USED ANYWHERE!
#     void PrintWeights( );
# 
#     // common, but Specialized by ModelObject1D
#     virtual void PopulateParameterNames( );
# 
#     // common, might be specialized...
#     virtual void FinalSetup( );
# 
#     // common, not specialized
#     string& GetParameterName( int i );
# 
#     // common, not specialized
#     int GetNFunctions( );
# 
#     // common, not specialized
#     int GetNParams( );
# 
#     // common, not specialized -- returns total number of data values (e.g., pixels in image)
#     int GetNDataValues( );
# 
#     // common, not specialized -- returns total number of *non-masked* data values
#     int GetNValidPixels( );
# 
#         // 2D only
#     double * GetModelImageVector( );
# 
#         // 2D only
#     double * GetExpandedModelImageVector( );
# 
#         // 2D only
#     double * GetResidualImageVector( );
# 
#         // 2D only
#     double * GetWeightImageVector( );
# 
#         // 2D only
#     double FindTotalFluxes(double params[], int xSize, int ySize, 
#                                                 double individualFluxes[] );
# 
#     // Generate a model image using *one* of the FunctionObjects (the one indicated by
#     // functionIndex) and the input parameter vector; returns pointer to modelVector.
#     double * GetSingleFunctionImage( double params[], int functionIndex );
# 
#     // 1D only
#     virtual int GetModelVector( double *profileVector ) { return -1; };
# 
#     // 1D or 2D
#     virtual void UseBootstrap( );
#     
#     // 1D or 2D
#     virtual void MakeBootstrapSample( );
#     
#     // Destructor
#     virtual ~ModelObject();
