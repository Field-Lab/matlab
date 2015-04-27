/*=================================================================
 *
 * mult_col_mult.c - Multiply columns of matrix A by weights supplied by a vector k then multiply by the original matrix transposed
 * M = A*diag(k)*A'
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

#include "gsl_cblas.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  double *A, *k, *M;
  double sum;
  int m,n,nm;
  
  int i,j,r,numtoupdate;
  
  
  if (nrhs != 2) {
      mexErrMsgTxt("need 2 input args: mult_col(A,k)\n");
  }
  
  /*  Read inputs into appropriate local variables */
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  nm = n*m;
  A = mxGetPr(prhs[0]);
  if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != 1) {
      mxErrMsgTxt("2nd argument k needs to be a col vector of length equal to num columns of first argument A!\n");
  }
  k = mxGetPr(prhs[1]);
  
  
  /* Initialize the output */
  plhs[0] = mxCreateDoubleMatrix(m,m,mxREAL);
  M = mxGetPr(plhs[0]);
  
  
  
  /* Only compute upper triangle */
  for (i=0; i<m; i++) {
      for (j=i; j < m; j++) {
          sum = 0;
          for (r=0; r < n; r++) {
            /* A_ir*k_r*A_jr */
            sum += A[m*r+i]*k[r]*A[m*r+j];
          }
          M[m*j+i] = sum;  
      }
  }

}
