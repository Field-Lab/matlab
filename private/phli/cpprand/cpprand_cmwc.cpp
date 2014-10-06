#include <limits>
#include "mex.h"
#include "mx_functional/mw_array_container.h"

class cmwc4096 {
public:
  typedef unsigned int result_type;

  explicit cmwc4096(result_type x0 = std::numeric_limits<int>::max()) {
    seed(x0);
  }

  void seed(result_type x0) {
    _Q[0] = x0;
    _Q[1] = x0 + _phi;
    _Q[2] = x0 + _phi + _phi;
    for (int i = 3; i < _state_size; ++i) _Q[i] = _Q[i-3] ^ _Q[i-2] ^ _phi ^ i;
  }

  result_type operator()() {
    unsigned long t, a = 18782LL;
    static unsigned int i = _state_size-1;
    result_type x, r = 0xfffffffe;
    i = (i + 1) & 4095;
    t = a * _Q[i] + c;
    c = (t >> 32);
    x = t + c;
    if (x < c) {
      x++;
      c++;
    }
    return (_Q[i] = r - x);
  }

  // friend functions
  template<typename CharT, typename Traits> 
  friend std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT,Traits>& os, const cmwc4096& cmwc) {
    for (int i = 0; i < _state_size; ++i)
      os << cmwc._Q[i] << " ";
    return os;
  }

  template<typename CharT, typename Traits> 
  friend std::basic_istream<CharT,Traits>& 
  operator>>(std::basic_istream<CharT,Traits>& is, cmwc4096& cmwc) {
    for (int i = 0; i < _state_size; ++i)
      is >> cmwc._Q[i] >> std::ws;
    return is;
  }

  static result_type min() { return __min; }
  static result_type max() { return __max; }

private:
  static const result_type __min = std::numeric_limits<result_type>::min();
  static const result_type __max = std::numeric_limits<result_type>::max();
  static const unsigned int _state_size = 4096;
  result_type _Q[_state_size];
  result_type _c = 362436;
  result_type _phi = 0x9e3779b9;
};
 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //if (nrhs > 2 && mxGetClassID(prhs[2]) == mxUINT32_CLASS && mxGetNumberOfElements(prhs[2]) == 1) {
      //init_rand(*static_cast<unsigned int *>(mxGetData(prhs[2])));
    //} else {
      //unsigned int seed = 2147483647;
      //init_rand(seed);
    //}

    //plhs[0] = mxCreateNumericMatrix(mxGetScalar(prhs[0]), mxGetScalar(prhs[1]), mxUINT32_CLASS, mxREAL);
    //MWArrayContainer<unsigned int> out(plhs[0]);
    //for (int i = 0; i < out.nelem(); ++i) {
      //out[i] = rand_cmwc();
    //}
}
