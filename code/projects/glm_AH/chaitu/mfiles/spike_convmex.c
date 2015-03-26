/*=================================================================
 * spike_convmex.c - Convolve a binary spike vector with a matrix 
 * whose columns are the postspike basis functions
 * A = spike_convmex(T,splist,B)
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  double *ps_basis, *h, *result, *splist;
  int nsp, T, hlen, nbasis;
  
  int i,j,k,numtoupdate,sptime,offset,basisoffset;
  
  
  if (nrhs != 3) {
      mexErrMsgTxt("need 2 input args: spike_convmex(T,spikes,postspike_basis)\n");
  }
  
  /*  Read inputs into appropriate local variables */
  T = (int)(*mxGetPr(prhs[0])); 
  nsp = mxGetM(prhs[1])*mxGetN(prhs[1]);
  splist = mxGetPr(prhs[1]);
  hlen = mxGetM(prhs[2]);
  nbasis = mxGetN(prhs[2]);
  ps_basis = mxGetPr(prhs[2]);
  
  /* Initialize the output */
  plhs[0] = mxCreateDoubleMatrix(T,nbasis,mxREAL);
  result = mxGetPr(plhs[0]);
  h = mxGetPr(plhs[0]);
  
  
  /* printf("T=%d,hlen=%d,nbasis=%d\n",T,hlen,nbasis); */
  
  /* Iterate through spikes */
  for (i=0;i<nsp;i++) {
    sptime = (int)(splist[i]) -1; /* time index of spike */
    numtoupdate = hlen;
    if (hlen >= T-sptime) {
        numtoupdate = T-sptime-1;
    }
    
    /* Convolve with each basis function */
    for (j=0; j < nbasis; j++) {
        /* memcpy((void *)(result+ j*T+sptime+1),(void *)(ps_basis + j*hlen),sizeof(double)*numtoupdate); */
        offset = j*T + sptime+1;
        basisoffset = j*hlen;
        for (k = 0; k<numtoupdate; k++) {
            result[offset+k] += ps_basis[basisoffset+k];            
        }
    }
  }
}
