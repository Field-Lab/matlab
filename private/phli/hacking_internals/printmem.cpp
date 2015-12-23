#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if (nrhs < 1 || !mxIsUint64(prhs[0]) || mxIsEmpty(prhs[0]))
    mexErrMsgTxt("First argument must be a uint64 memory address.");
  unsigned long *addr = static_cast<unsigned long *>(mxGetData(prhs[0]));
  unsigned char *mem = (unsigned char *) addr[0];
  
  if (nrhs < 2 || !mxIsDouble(prhs[1]) || mxIsEmpty(prhs[1]))
    mexErrMsgTxt("Second argument must be a double-type integer byte size.");      
  unsigned int nbytes = static_cast<unsigned int>(mxGetScalar(prhs[1]));
      
  for (int i = 0; i < nbytes; i++) {
      printf("%.2x ", mem[i]);
      if ((i+1) % 16 == 0) printf("\n");
  }
  printf("\n");
}