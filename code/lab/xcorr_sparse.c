#include <stdio.h>
#include <math.h>
#include "mex.h"

/* 
   % XCORR_SPARSE  Compute cross-correlation between discrete time events
   %
   % usage:  T = dxcorr(T1, T2, min, max, delta);
   %
   % Cross-correlation between discrete events. Far faster than XCORR.
   %
   % Inputs:      T1, T2 - column vector of events
   %            min, max - range of correlation window
   %               delta - step size in output histogram
   %
   % All arguments are assumed to be in the same time units.
   %
   % Mex file.
   %
*/



int histogram_size(int start, int end, int delta) {
  return (int) ceil((end - start) / delta) + 1;
}

void increment_histogram(int* histogram, int N, int value, 
			 int start, int end, int delta) {
  int index = (value - start) / delta;
  //int index = (int) floor((value - start) / delta);  // CHECK OUT
  histogram[index]++;
} 



/* XCORR_DISCRETE
   Compute cross-correlation between discrete time events. Replaces Matlab function
   XCORR when the process is mainly sparse. Useful for spike trains.
*/

void xcorr_sparse(int* t1, int N1, int* t2, int N2, int* histogram, int N, 
	    int start, int end, int delta) {
  
  int i=0, j_init=0, j=0;      // indices
  int min_value, max_value; // range in histogramming
  int value;                // identified value

  // create and initialize histogram
  for (i=0;i< N; i++) histogram[i] = 0;
  
  // for each spike in t1
  // Note: t1[i], t2[j]
  for (i=0; i < N1; i++) {

    // determine acceptable range for correlation
    min_value = t1[i] + start;
    max_value = t1[i] + end + delta;

    // start t2 index at previous spike
    j = j_init;

    // increment t2 spikes until find candidate
    // Note: order of test is important!
    while ((j < N2) && (t2[j] < min_value)) j++;
    
    // determine starting point for t2 on the next t1 spike
    j_init = j - 1;
    if (j_init < 0) j_init = 0;

    // grab all spikes in t2 within interval
    // Note: order of test is important!
    while ((j < N2) && (t2[j] < max_value)) {

      // determine candidate value to include in histogram
      value = t2[j] - t1[i];
      // increment histogram value
      increment_histogram(histogram, N, value, start, end, delta);
      j++;

    }
  }
}


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {

  int N1, N2, N;          // length of spike trains 1, 2, correlation histogram
  int *T1, *T2;           // spike trains 1, 2
  int* hist;              // cross-correlation histogram
  int i;
  int min, max, delta;  // window for detecting correlations


  /* PRELIMINARY DATA CHECKS */
  /* ======================= */
  
  /* 1. Number of arguments */
  if ( nrhs != 5 ) 
    mexErrMsgTxt("xcorr_sparse: incorrect number of arguments");

  if (nlhs > 1) 
    mexErrMsgTxt("xcorr_sparse: only outputs 2 item at max.");
    
  /* 2. Check arguments are appropriate flavor */
  for (i=0; i<nrhs; i++)
    if ( mxIsComplex(prhs[i]) ||
         mxIsClass(prhs[i],"sparse") || 
         mxIsChar(prhs[i]) )
      mexErrMsgTxt("xcorr_sparse: all arguments must be real, full and non-string");
  
  /* 3. Check data sizes */
  for (i=0; i<nrhs; i++) {
    if (mxGetN(prhs[i]) != 1) 
      mexErrMsgTxt("xcorr_sparse: all RHS args must have 1 column.");
  }

  for (i=2; i<nrhs; i++) {
    if (mxGetM(prhs[i]) != 1) 
      mexErrMsgTxt("xcorr_sparse: all RHS must be scalars.");
  }
  
  /* 4. Check integers */
  for (i=0;i<2;i++) {
    if (!mxIsInt32(prhs[i])) {
      mexErrMsgTxt("xcorr_sparse: double/float not implemented. Must be int32");
    }
  }

  /* GRAB INFORMATION */
  /* ================ */

  // get length of vector
  N1 = mxGetNumberOfElements(prhs[0]); 
  N2 = mxGetNumberOfElements(prhs[1]);

  // grab spike train
  T1 = (int*) mxGetData(prhs[0]);
  T2 = (int*) mxGetData(prhs[1]);

  // grab paramters  
  min = (int) mxGetScalar(prhs[2]);
  max = (int) mxGetScalar(prhs[3]);
  delta = (int) mxGetScalar(prhs[4]);

  // check arguments are correct
  if (min >= max)
    mexErrMsgTxt("xcorr_sparse: min >= max");
  
  if (delta > (max-min))
    mexErrMsgTxt("xcorr_sparse: delta > window");

  /* ALLOCATE MEMORY */
  /* =============== */

  // determine histogram size
  N = histogram_size(min, max, delta);

  // controlled by Matlab
  hist = (int *) mxMalloc( N * sizeof(int));
 
  if (hist == NULL)
    mexErrMsgTxt("xcorr_sparse: unable to allocate memory");


  /* RUN COMPUTE-COINCIDENCE */
  /* ======================= */
  xcorr_sparse(T1, N1, T2, N2, hist, N, min, max, delta);

  /* SAVE OUTPUT DATA */
  /* ================ */

  // return correlation window
  plhs[0] = mxCreateNumericMatrix(N, 1, mxINT32_CLASS, mxREAL);
  mxSetData(plhs[0], hist);

}

