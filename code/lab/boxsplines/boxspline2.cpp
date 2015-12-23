/* boxspline2_parallel.cpp
 * Version 0.8
 *
 * Taken directly from:
 *   Three-Directional Box-Splines: Characterization and Efficient
 *   Evaluation, Laurent Condat and Dimitri Van De Ville, IEEE Signal
 *   Processing Letters, 2006
 * with minor optimizations.
 *
 * 2011-09 Peter H. Li
 * FreeBSD License
 */


#include "mex.h"
#include "math.h"
        
float boxspline2(float x, float y) {
  float u = fabs(x) - fabs(y)/sqrt(3.0);
  float v = fabs(x) + fabs(y)/sqrt(3.0);
  
  if (u < 0) {
    u = -u;
    v = v + u;
  }
  if (2*u < v) u = v - u;
  
  if (v > 2.0) return 0.0; // Outside the support
  
  float g = u - v/2.0;
  if (v < 1.0) return 0.5 + ((5/3.0 - v/8.0)*v - 3)*v*v/4.0 + ((1 - v/4.0)*v + g*g/6.0 - 1)*g*g;
  if (u > 1.0) return (v - 2) * (v - 2) * (v - 2) * (g - 1) / 6.0;
  return 5/6.0 + ((1 + (1/3.0 - v/8.0)*v)*v/4.0 - 1)*v + ((1 - v/4.0)*v + g*g/6.0 - 1)*g*g;
}
        

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check inputs
  if (nrhs != 2)
    mexErrMsgIdAndTxt("boxspline2:nrhs", "Arguments should be the x and y coordinate matrices");
  if (!mxIsDouble(prhs[0]))
    mexErrMsgIdAndTxt("boxspline2:prhs", "First argument must be a double matrix.");
  if (!mxIsDouble(prhs[1]))
    mexErrMsgIdAndTxt("boxspline2:prhs", "Second argument must be a double matrix.");
  
  // Check that x and y have same number of elements
  const mwSize numelx = mxGetNumberOfElements(prhs[0]);
  const mwSize numely = mxGetNumberOfElements(prhs[1]);
  if (numelx != numely)
    mexErrMsgIdAndTxt("boxspline2:prhs", "First and second arguments must have same number of arguments");
  
   double *inx = mxGetPr(prhs[0]);
   double *iny = mxGetPr(prhs[1]);
   plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
   double *out = mxGetPr(plhs[0]);
   for (int i = 0; i < numelx; i++) out[i] = boxspline2(inx[i], iny[i]);
}