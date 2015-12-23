/**
 * To compile:
 * /home/peterli/cuda-4.0.17/bin/nvcc -I/home/peterli/cuda-4.0.17/include -c thrust_sort_lib.cu -Xcompiler -fPIC
 *
 * Note, CUDA toolset version must match that included in Matlab if this is to
 * be called from Matlab Mex.
 *
 * -fPIC needed for call from Matlab Mex as well.
 *
 */

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include "mx_functional/skip_ptr.h"


/**
 * thrust_device_sort
 */
template<typename T>
void thrust_device_sort(T first, T last) {
  typedef typename T::value_type V;
  thrust::device_vector<V> dv(first, last);
  thrust::sort(dv.begin(), dv.end());
  thrust::copy(dv.begin(), dv.end(), first);
}

// Explicit instantiations
template void thrust_device_sort<skip_ptr<double> >(skip_ptr<double>, skip_ptr<double>);
template void thrust_device_sort<skip_ptr<float> >(skip_ptr<float>, skip_ptr<float>);
template void thrust_device_sort<skip_ptr<int> >(skip_ptr<int>, skip_ptr<int>);
template void thrust_device_sort<skip_ptr<unsigned int> >(skip_ptr<unsigned int>, skip_ptr<unsigned int>);
template void thrust_device_sort<skip_ptr<short> >(skip_ptr<short>, skip_ptr<short>);
template void thrust_device_sort<skip_ptr<unsigned short> >(skip_ptr<unsigned short>, skip_ptr<unsigned short>);
template void thrust_device_sort<skip_ptr<char> >(skip_ptr<char>, skip_ptr<char>);
template void thrust_device_sort<skip_ptr<unsigned char> >(skip_ptr<unsigned char>, skip_ptr<unsigned char>);



/**
 * thrust_sort
 */
template<typename T>
void thrust_sort(T first, T last) {
  thrust::sort(first, last);
}

// Explicit instantiations
template void thrust_sort<skip_ptr<double> >(skip_ptr<double>, skip_ptr<double>);
template void thrust_sort<skip_ptr<float> >(skip_ptr<float>, skip_ptr<float>);
template void thrust_sort<skip_ptr<int> >(skip_ptr<int>, skip_ptr<int>);
template void thrust_sort<skip_ptr<unsigned int> >(skip_ptr<unsigned int>, skip_ptr<unsigned int>);
template void thrust_sort<skip_ptr<short> >(skip_ptr<short>, skip_ptr<short>);
template void thrust_sort<skip_ptr<unsigned short> >(skip_ptr<unsigned short>, skip_ptr<unsigned short>);
template void thrust_sort<skip_ptr<char> >(skip_ptr<char>, skip_ptr<char>);
template void thrust_sort<skip_ptr<unsigned char> >(skip_ptr<unsigned char>, skip_ptr<unsigned char>);
