#include "mex.h"
#include <fstream>

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if (nrhs < 1)
    mexErrMsgTxt("Require an array address to dump from");
  char *dumpaddress = (char *) prhs[0];
  
//  if (nrhs > 1 && mxIsScalar(prhs[1]) {
//  }
  
  if (nrhs > 2 && mxIsChar(prhs[2])) {
    filename = mxArrayToString(prhs[1]);  
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    file.write(dumpaddress, 104);
    file.close();
  } else {
  }
}