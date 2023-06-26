// wb_compute_cost.cpp
//
// J = wb_compute_cost(tau, S.DCM, fc, fs, h, t0, dt);
//
// mex function used to compute cost for a possible DOA solution.  Called
// by wb_cost_function.m and wb_initialization.m.
//
// To run, first compile wb_compute_cost.cpp from the
// MATLAB command line as follows:
//
// mex -v -largeArrayDims -lmwlapack -lmwblas wb_compute_cost.cpp
//   -v = verbose (not necessary)
//   -largeArrayDims required for LAPACK and BLAS libraries
//
// Then use the following syntax from the MATLAB command line:
// J = wb_compute_cost(tau, S.DCM, fc, fs, h, t0, dt);
//
// Linear Algebra PACKage (LAPACK) Information (used for matrix inverse, zgetrf and zgetri):
//  http://www.netlib.org/lapack/explore-html/d7/d99/_v_a_r_i_a_n_t_s_2lu_2_c_r_2zgetrf_8f.html
//  edit([matlabroot '/extern/include/lapack.h'])
//
// For cross platform compatibility, use size_t and ptrdiff_t for
// referencing pointers.
//
//
// Inputs:
//  tau   = Nc x Nsig matrix of relative two way time delays for each 
//          channel given a possible DOA solution.  Nc is the number of 
//          sensors and Nsig is the user specified dimensionality of the 
//          signal subspace.  Let 1 <= i <= Nc and 1 <= j <= Nsig denote
//          row and column indexes of the tau matrix and let [tau]i,j 
//          denote the (i,j) entry of tau.  [tau]i,j is the two way 
//          propagation between the origin of the SAR flight coordinate 
//          system and the ith sensor for the jith source, 
//  DCM   = space-time data covariance matrix built in array_proc.m,
//  fc    = carrier frequency at passband,
//  fs    = sampling frequency at basedband after CSARP,
//  h     = impulse response time series built in combine_wf_chan_task.m
//  t0    = first time value associated with the first value in h,
//  dt    = sampling interval of h.
//
// Outputs:
//  J     = real valued scalar in linear scale describing the cost of a
//          given DOA solution specified by the input tau.
//
// Authors:  John Paden, Theresa Stumpf
//
// See Also:  array_proc.m, wb_initialization.m, wb_cost_function.m
// ========================================================================

#if defined( _WIN32 ) || defined( _WIN64 )
#define _USE_MATH_DEFINES // For windows
#endif

#include <math.h>
#include <mex.h>
#include <lapack.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  // ----------------------------------------------------------------------
  // Declare Pointers to Input arguments and Initialize
  // ----------------------------------------------------------------------
  // tau  = num_sens x num_src array whose entries correspond to the 
  //        2 way time delays to each sensor for a particular DOA solution.
  //        tau = (2/c) * (y_pc*sin(DOA) - z_pc*cos(DOA)), where y_pc, z_pc 
  //        are phase centers from FCS, and DOA is 1xnum_src vector 
  //        containing a possible solution,
  // DCMr = real part of the space time data covariance matrix,
  // DCMi = imag part of the space time data covariance matrix,
  // Fc   = 1x1 scalar whose value is the center frequency,
  // Fs   = 1x1 scalar whose value is the sampling frequency AFTER SAR 
  //        processing,
  // h    = vector containing samples of the modeled impulse response,
  // t0   = 1x1 scalar whose value corresponds to the initial time of h,
  // dt   = sampling interval of h
  // NOTE:  given t0 and dt, we know the time vector associated with h and
  //        compute these values as needed instead of passing in the time 
  //        vector.
  
  double *tau, *DCMr, *DCMi, *Fc, *Fs, *h, *t0, *dt;
  tau   = mxGetPr(prhs[0]);
  DCMr  = mxGetPr(prhs[1]);
  DCMi  = mxGetPi(prhs[1]);
  Fc    = mxGetPr(prhs[2]);
  Fs    = mxGetPr(prhs[3]);
  h     = mxGetPr(prhs[4]);
  t0    = mxGetPr(prhs[5]);
  dt    = mxGetPr(prhs[6]);
  
  // Input checking
  if (!mxIsComplex(prhs[1]))
  {
    mexPrintf("DCM must be complex\n");
    return;
  }
  
  if (!mxIsComplex(prhs[1]))
  {
    mexPrintf("DCM must be complex\n");
    return;
  }
  
  //  Convert carrier frequency to radians
  double Omega_c = *Fc*2*M_PI; 
  
  // ----------------------------------------------------------------------
  // Setup dimension variables
  // ----------------------------------------------------------------------
  // num_sens = number of sensors in the array
  // num_meas = number of rows in the space time covariance matrix.  
  //            num_meas = num_sens * Widening factor.
  // num_cols  = number of columns of the space time covariance matrix.  
  //            The DCM is square by definition but defining num_cols 
  //            separately allows a subset of columns to be passed in.  
  //            Sometimes it's useful to look at the cost function of each 
  //            column separately to debug.
  // num_src  = number of sources, determined by the number of columns in 
  //            tau
  
  ptrdiff_t num_meas  = (ptrdiff_t)mxGetM(prhs[1]);
  ptrdiff_t num_cols   = (ptrdiff_t)mxGetN(prhs[1]);
  ptrdiff_t num_sens  = (ptrdiff_t)mxGetM(prhs[0]);
  ptrdiff_t num_src   = (ptrdiff_t)mxGetN(prhs[0]);
  
  // declare pointer to output and allocate memory
  double *Jr, *Ji;
  Jr    = (double *)mxMalloc(1 * sizeof(double));
  Ji    = (double *)mxMalloc(1 * sizeof(double));
  
  // allocate memory for A.  A is the model of the space-time covariance 
  // matrix 
  double *bigAr, *bigAi;
  bigAr = (double *)mxMalloc(num_meas *num_meas * num_src * sizeof(double));
  bigAi = (double *)mxMalloc(num_meas *num_meas * num_src * sizeof(double));
  
  // matrix used to prepare for inverse (A' * A) stored in smallA.
  double *smallA;
  smallA = (double *)mxMalloc(num_src * num_src * 2 * sizeof(double));
  
  // First matrix multiply
  double *mat1r, *mat1i;
  mat1r = (double *)mxMalloc(num_src * 1 * sizeof(double));
  mat1i = (double *)mxMalloc(num_src * 1 * sizeof(double));
  
  // Second matrix multiply
  double *mat2r, *mat2i;
  mat2r = (double *)mxMalloc(num_src * 1 * sizeof(double));
  mat2i = (double *)mxMalloc(num_src * 1 * sizeof(double));
  
  // Third matrix multiply
  double *mat3r, *mat3i;
  mat3r = (double *)mxMalloc(num_src * 1 * sizeof(double));
  mat3i = (double *)mxMalloc(num_src * 1 * sizeof(double));
  
  // ----------------------------------------------------------------------
  // Setup dimension for accessing appropriate locations in memory
  // ----------------------------------------------------------------------
  ptrdiff_t ref_index     = 0; //
  ptrdiff_t mSensSkip     = 0; // used to skip num_sens when looping over m
  ptrdiff_t mMeasSkip     = 0; // used to skip num_meas when looping over m
  ptrdiff_t jMeasSrcSkip  = 0; // used to skip num_meas*num_src when 
                               // looping over j
  
  double tArg1      = 0.0;  // time argument used to evaluate phase term
  double tArg2      = 0.0;  // time argument used to evaluate impulse resp
  double shift      = 0.0;  // integer value 
  double t1         = 0.0;
  double h1         = 0.0;
  double h2         = 0.0;
  double H          = 0.0;
  ptrdiff_t n1 = 0;
  
  // J(tau(theta)) = 
  //    sum over num_cols{DCM(:,col)'*DCM(:,col) - 
  //                        DCM(:,col)'*Pa(col;tau)*DCM(:,col)} 
  // where Pa(col;tau) is the projection matrix of a specific column given 
  // tau.  Pa(col;tau) = A * inv(A'*A) * A' with A being the matrix whose 
  // columns model a specific column of the DCM model in response to tau.
  // A is unique for each column of the DCM model so in contrast to MLE, 
  // this method uses a unique projection matrix for each column/tau 
  // combination.
  
  // ======================================================================
  // PART 1:  SETUP MODEL OF COLUMNS OF THE SPACE TIME COVARIANCE MATRIX
  // FOR THE GIVEN tau AND STORE IN bigA
  // ======================================================================
  // Create the bigAr (expectation of the columns of the data covariance 
  // matrix, DCM)
  // bigA = [num_meas x colsDCM x num_src]
  
  ptrdiff_t MeasSrcSkip       = num_meas*num_src;
  for (ptrdiff_t j=0; j<num_cols; j++) { // loop over columns of DCM model
    
    // Index to access appropriate memory location for the jth column of 
    // bigA (have to skip over num_meas*num_src)
    jMeasSrcSkip       = j*MeasSrcSkip; 
    
    for (ptrdiff_t m=0; m<num_src; m++) { 
      mSensSkip = m*num_sens;        
      mMeasSkip = m*num_meas;
      ref_index = j%num_sens + mSensSkip;
      
      for (ptrdiff_t i=0; i<num_meas; i++) {
        // Setup the shift variable used to evaluate the space-time 
        // covariance at different lags.  ie let X(t) be num_sens by 1 
        // vector of array data.  We assume the array samples a WSS random
        // process.  Then Rxx(0) = E(X(t)*X'(t)) is the more common spatial
        // covariance matrix, defined as the multichannel covariance matrix
        // zero lag.  We use the shift variable to evaluate the 
        // multichannel covariance matrix at different lags.  Each lag is a 
        // multiple of 1/sampling frequency or Ts.  So the case of shift=1
        // corresponds to Rxx(Ts) = E(X(t)*X'(t-Ts)).  The shift variable
        // takes integer values over the interval 
        // -floor(W-1)/2 : floor(W-1)/2, where W is the widening factor.
        // For the Basler, W = 3 and shift takes a value from the set
        // {-1,0,1}.
        shift   = (double)(j/num_sens - i/num_sens); 
//         mexPrintf("shift = %f \n",shift);
        
        // Setup time argument used to evaluate the phase of the (i,j) 
        // entry of the space time covariance model for the mth source.
        // NOTE: tArg1 must agree with the convention of FCS. The sign of 
        // tArg1 may need to be flipped to agree.
        tArg1   = tau[i%num_sens + mSensSkip] - tau[ref_index]; // for tau = y_pc*sin(theta) - z_pc*cos(theta) (!)
                
        // Setup time argument used to evaluate the magnitude term of the 
        // (i,j) entry of the space time covariance model of the mth 
        //source.
        tArg2   = tArg1 - (shift / *Fs);  
        // Use linear interpolation to evaluate the impulse response h at
        // tArg2.
        n1      = (ptrdiff_t)((tArg2 - *t0) / (*dt));
        t1      = *t0 + *dt * n1;
        h1      = h[n1];
        h2      = h[n1 + 1];
        H       = h1 + (h2-h1) * ((tArg2 - t1) / (*dt));
        
        // Compute corresponding entry in bigA
        bigAr[i + mMeasSkip + jMeasSrcSkip] = H*cos(Omega_c*tArg1);
        bigAi[i + mMeasSkip + jMeasSrcSkip] = H*sin(Omega_c*tArg1);      
      }
    }
  }
  
  // ======================================================================
  // PART 2:  SETUP  A'*A TERM OF THE PROJECTION MATRIX AND STORE IN smallA 
  // FOR INVERSION.  AT SAME TIME, COMPUTE THE PORTION OF J THAT IS THE SUM
  // OF INNER PRODUCTS OF EACH COLUMN OF THE DCM
  // ======================================================================
  
  // Initialize cost function to zero
  *Jr = 0.0;
  *Ji = 0.0;
  
  // Reinitialize measurement skip variable to 0
  mMeasSkip = 0;
  
  // Declare sumr and sumi, variables used to store real and imaginary 
  // parts of the A'*A matrix multiply
  double sumr = 0, sumi = 0;
  ptrdiff_t iMeasSkip = 0, jMeasSkip = 0, jSrcSkip = 0, mMeasSrcSkip = 0;
  
  for (ptrdiff_t m = 0; m < num_meas; m++) { // loop over columns of DCM
    mMeasSkip     = m*num_meas;    // index for accessing memory location 
    mMeasSrcSkip  = m*MeasSrcSkip; // index for accessing memory location 
    
    // sum up inner product of each column of the DCM (first term in cost
    // function equation)
    for (int m2 = 0; m2 < num_meas; m2++) {
      *Jr = *Jr + DCMr[m2 + mMeasSkip]*DCMr[m2 + mMeasSkip] + DCMi[m2 + mMeasSkip]*DCMi[m2 + mMeasSkip];
    }
    
    // perform matrix multiply A'*A for each column modeled
    for (int i = 0; i < num_src; i++) { // loop over columns of A'
      iMeasSkip = i*num_meas;
      for (int j = 0; j < num_src; j++) { // loop over columns of A
        sumr = 0;
        sumi = 0;
        jMeasSkip = j*num_meas;
        jSrcSkip = j*num_src;
        // multiply and sum appropriate entries
        for (int k = 0; k < num_meas; k++) {
          
          sumr += bigAr[k + iMeasSkip  + mMeasSrcSkip] * bigAr[k + jMeasSkip + mMeasSrcSkip]
                  + bigAi[k + iMeasSkip + mMeasSrcSkip] * bigAi[k + jMeasSkip + mMeasSrcSkip];
          
          sumi += bigAr[k + iMeasSkip + mMeasSrcSkip] * bigAi[k + jMeasSkip + mMeasSrcSkip]
                  - bigAi[k + iMeasSkip + mMeasSrcSkip] * bigAr[k + jMeasSkip + mMeasSrcSkip];
        }
        // store in smallA and interleave real and imaginary components 
        // (necessary for zgetrf and zgetri)
        smallA[2*(i+jSrcSkip)]      = sumr;
        smallA[2*(i+jSrcSkip) + 1]  = sumi;
      }
    }
    
    // declare and initialize variable needed to do matrix inverse
    ptrdiff_t *ipiv, lwork, info;
    double *wk;
    lwork   = 10*num_src;
    ipiv    = (ptrdiff_t *)mxMalloc(2*num_src*sizeof(*ipiv));
    wk      = (double *)mxMalloc(2*lwork*sizeof(*wk));
    // Do the matrix inversion of smallA and store back in the same memory 
    // location
    zgetrf(&num_src, &num_src, smallA, &num_src, ipiv, &info);
    zgetri(&num_src, smallA, &num_src, ipiv, wk, &lwork, &info);
    
    //=====================================================================
    // PART 3: COMPUTE REMAINING MATRIX MULTIPLIES IN COST FUNCTION 
    // EXPRESSION.  THIS PART OF THE COST FUNCTION CORRESPONDS TO THE SUM
    // OF THE L2 NORMS OF THE PROJECTIONS OF EACH COLUMN OF THE DCM ONTO 
    // ITS CORRESPONDING MODEL (GIVEN BY bigA)
    //=====================================================================
    //
    //  DCM(:,column)'*bigA  * smallA * bigA' * DCM(:,column)
    // |______mat1_________|  |         |_______mat2________|
    //                        |__________mat3_______________|
    //
    // 1) create mat1 = DCM(:,column)'*bigA
    for (ptrdiff_t i = 0; i < num_src; i++) {
      sumr = 0;
      sumi = 0;
      iMeasSkip = i*num_meas;
      for (ptrdiff_t j = 0; j< num_meas; j++) {
        
        sumr += DCMr[j + mMeasSkip] * bigAr[j + iMeasSkip + mMeasSrcSkip]
                + DCMi[j + mMeasSkip] * bigAi[j + iMeasSkip + mMeasSrcSkip];
        
        sumi += DCMr[j + mMeasSkip] * bigAi[j + iMeasSkip + mMeasSrcSkip]
                - DCMi[j + mMeasSkip] * bigAr[j + iMeasSkip + mMeasSrcSkip];
      }
      mat1r[i] = sumr;
      mat1i[i] = sumi;
    }
    
    // 2) create mat2 = bigA'*DCM(:,column)
    for (ptrdiff_t i = 0; i < num_src; i++) {
      sumr = 0;
      sumi = 0;
      iMeasSkip = i*num_meas;
      for (ptrdiff_t j = 0; j < num_meas; j++) {
        sumr += bigAr[j + iMeasSkip + mMeasSrcSkip] * DCMr[j + mMeasSkip]
                +bigAi[j + iMeasSkip + mMeasSrcSkip] * DCMi[j + mMeasSkip];
        sumi += bigAr[j + iMeasSkip + mMeasSrcSkip] * DCMi[j + mMeasSkip]
                -bigAi[j + iMeasSkip + mMeasSrcSkip] * DCMr[j + mMeasSkip];
      }
      mat2r[i] = sumr;
      mat2i[i] = sumi;
    }
    
    // 3) create mat3 = smallA * mat2
    for (ptrdiff_t i = 0; i < num_src; i++) {
      sumr = 0;
      sumi = 0;
      for (ptrdiff_t j = 0; j < num_src; j++) {
        jSrcSkip = j*num_src;
        sumr += smallA[2*(i + jSrcSkip)] * mat2r[j]
                -smallA[2*(i + jSrcSkip) + 1] * mat2i[j];
        sumi += smallA[2*(i + jSrcSkip)] * mat2i[j]
                +smallA[2*(i + jSrcSkip) + 1] * mat2r[j]; 
      }
      mat3r[i] = sumr;
      mat3i[i] = sumi;
    }
    
    // 4) Compute mat1*mat3 and sum
    for (ptrdiff_t i = 0; i < num_src; i++) {
//       *Jr = *Jr + mat1r[i]*mat3r[i] - mat1i[i]*mat3i[i];
      *Jr = *Jr - mat1r[i]*mat3r[i] + mat1i[i]*mat3i[i];
      
    }
  }
  
  *Jr = *Jr / num_cols;
  
  plhs[0] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxCOMPLEX);
  mxSetPr(plhs[0], Jr);
  mxSetPi(plhs[0], Ji);
  mxSetM(plhs[0], 1);
  mxSetN(plhs[0], 1);
  
}

