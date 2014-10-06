#include "mex.h"
#include "string.h"

const char *byte_to_binary(unsigned char c)
{
    char b[9];
    b[0] = '\0';

    int z;
    for (z = 128; z > 0; z >>= 1) {
        strcat(b, ((c & z) == z) ? "1" : "0");
    }

    return b;
}


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
      printf("%s ", byte_to_binary(mem[i]));
      if ((i+1) % 8 == 0) printf("\n");
  }
  printf("\n");
}