/*=================================================================
 *
 * stim_design.c - Construct the stimulus portion of the design matrix  - this is a simplified version of convmtx in MATLAB
 * X = stim_design(stim,mk)
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  double *X, *stim;
  int n,T,m;
  
  int i,j,k,numtoupdate;
  
  
  if (nrhs != 2) {
      mexErrMsgTxt("need 2 input args: stim_design(X,m)\n");
  }
  
  /*  Read inputs into appropriate local variables */
  n = mxGetM(prhs[0]);
  T = mxGetN(prhs[0]);
  stim = mxGetPr(prhs[0]);
  m = (int)(*mxGetPr(prhs[1]));
  
  
  /* Initialize the output */
  plhs[0] = mxCreateDoubleMatrix(n*m,T-m+1,mxREAL);
  X = mxGetPr(plhs[0]);
  
  /* printf("T=%d,n=%d,m=%d\n",T,n,m); */
  
 
  for (i=0;i<T-m+1;i++) { /* Iterate through cols */
      for (j=0; j < m; j++) {
        memcpy((void *)(X+i*n*m+j*n),(void *)(stim+(i+m-1-j)*n),n*sizeof(double));
      }
  }
}
