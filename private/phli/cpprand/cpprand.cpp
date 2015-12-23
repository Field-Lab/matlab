/* cpprand.cpp
 * Ver 0.9.9
 * Peter H. Li 2011 FreeBSD License 
 * See cpprand.m for documentation
 */

#include "mex.h"
#include "cpprand_struct2dist.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check inputs
  if (nrhs < 3)
    mexErrMsgIdAndTxt("cpprand:nrhs", "Arguments should be M and N dimensions, a starting seed/state, and an optional distribution");
  if (!mxIsNumeric(prhs[0]) || !mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[0]) != 1 || mxGetNumberOfElements(prhs[1]) != 1)
    mexErrMsgIdAndTxt("cpprand:prhs", "First and second arguments should be scalars giving the M and N dimensions");
  
  // 3rd argument can either be a uint32 seed, or a char array of the full state
  CppRandFunct<> cpprandfunct;
  switch (mxGetClassID(prhs[2])) {
    case mxUINT32_CLASS:
      if (mxGetNumberOfElements(prhs[2]) != 1)
        mexErrMsgIdAndTxt("cpprand:prhs", "3rd argument should be a scalar uint32 seed, or a char array of previous state");
      cpprandfunct = CppRandFunct<>(*static_cast<unsigned int *>(mxGetData(prhs[2])));
      break;
    case mxCHAR_CLASS:
      cpprandfunct = CppRandFunct<>(mxArrayToString(prhs[2]));
      break;
    default:
      mexErrMsgIdAndTxt("CPPRAND:prhs", "3rd argument should be uint32 seed or state string");
  }

  // Create output matrix
  plhs[0] = mxCreateNumericMatrix(mxGetScalar(prhs[0]), mxGetScalar(prhs[1]), mxDOUBLE_CLASS, mxREAL);

  // Get random numbers, save end state of RNG
  std::string outstate;
  if (nrhs > 3) { 
    // Distribution specified as Matlab struct; decode into C++ distribution
    CppRandDistFunct<> crdf(cpprandfunct, MWArrayContainer<double>(plhs[0]));
    outstate = diststruct_apply<CppRandDistFunct<>, std::string>(prhs[3], crdf);

  } else {
    // Use default distribution
    CPPRAND_RANDOM_NAMESPACE::uniform_01<> dist;
    outstate = cpprandfunct(MWArrayContainer<double>(plhs[0]), dist);
  }
  
  // Copy ending state string to output if requested
  if (nlhs > 1) plhs[1] = mxCreateString(outstate.c_str());
}
