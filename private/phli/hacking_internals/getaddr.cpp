#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if (nrhs < 1) 
    mexErrMsgTxt("One input required.");
  
  plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  unsigned long *out = static_cast<unsigned long *>(mxGetData(plhs[0]));
  out[0] = (unsigned long) prhs[0];
}