#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if (nrhs < 1) 
    mexErrMsgTxt("One input required.");
  printf("%p\n", prhs[0]);
}