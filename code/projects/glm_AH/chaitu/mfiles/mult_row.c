/*=================================================================
 *
 * mult_row.c - Multiply rows of matrix A by weights supplied by a vector k
 * M = mult_row(A,k)
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  double *A, *k, *M, *stim;
  int m,n,nm;
  
  int i,j,numtoupdate;
  
  
  if (nrhs != 2) {
      mexErrMsgTxt("need 2 input args: mult_col(A,k)\n");
  }
  
  /*  Read inputs into appropriate local variables */
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  nm = n*m;
  A = mxGetPr(prhs[0]);
  if (mxGetM(prhs[1]) != m || mxGetN(prhs[1]) != 1) {
      mxErrMsgTxt("2nd argument k needs to be a col vector of length equal to num rows of first argument A!\n");
  }
  k = mxGetPr(prhs[1]);
  
  
  /* Initialize the output */
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  M = mxGetPr(plhs[0]);
  
    /* Copy entries of A to M */
  memcpy((void *)(M),(void *)(A),nm*sizeof(double));
  
  for (i=0;i<nm;i++) {
    M[i] = M[i]*k[(i%m)];
  }

}
