/**
 * To compile:
 * mex -I/home/peterli/cuda-4.0.17/include thrust_sort.cpp -L/home/peterli/cuda-4.0.17/lib64 -lcudart thrust_sort_lib.o
 *
 * Note, CUDA toolkit version must match that included in Matlab.
 *
 */

#include "mex.h"
#include "thrust/device_vector.h"
#include "thrust_sort_lib.h"
#include "mx_functional/mx_functional.h"

class ThrustSortFunct
{
public:
  template <typename D>
  void operator()(MWArrayContainer<D> a) const {
    thrust_device_sort(a.begin(), a.end());
  }
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  plhs[0] = mxDuplicateArray(prhs[0]);
  ThrustSortFunct *tsf = new ThrustSortFunct();
  ColApplyFunct<ThrustSortFunct> f(*tsf);
  numarr_apply<void>(plhs[0], f);
}
