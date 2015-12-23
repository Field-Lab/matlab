/* cpprandpar.cpp
 * Ver 0.9.9
 * Peter H. Li 2011 FreeBSD License 
 * See cpprandpar.m for documentation.  Compare with cpprand.cpp for basic logic.
 */

#include "mex.h"
#include "cpprand_struct2dist.hpp"
#include "mx_functional/mx_functional_parallel.h"


// Adapter between BlockApplyFunct and CppRandFunct
template <class D = CPPRAND_RANDOM_NAMESPACE::uniform_01<> >
class CppRandBlockFunct
{
public:
  CppRandBlockFunct(MWArrayContainer<unsigned int> seeds, MWArrayContainer<double> out, std::string *states, D dist) :
    _seeds(seeds), _out(out), _states(states), _dist(dist) {}
    
  void operator()(int blocknum, int nblocks) const {
    MWArrayContainer<double> outblock = _out.get_block(blocknum, nblocks);
    CppRandFunct<> crf;

    // If seeds were given, use them, otherwise assume full states were given
    if (_seeds.nelem() > 0) crf = CppRandFunct<>(_seeds[blocknum]);
    else crf = CppRandFunct<>(_states[blocknum]);
    
    // Call rand generate functor and save ending state strings
    D dist = _dist;
    _states[blocknum] = crf(outblock, dist);
  }
  
  std::string *get_states() const { return _states; }
  
private:
  MWArrayContainer<unsigned int> _seeds;
  MWArrayContainer<double> _out;
  std::string *_states;
  D _dist;
};


// Adapter between StructToDistFunct and CppRandBlockFunct
class CppRandDistBlockFunct
{
public:
  CppRandDistBlockFunct(MWArrayContainer<unsigned int> seeds, MWArrayContainer<double> out, std::string *states, size_t n) :
    _seeds(seeds), _out(out), _states(states), _n(n) {}
  
  template <typename D>
  void operator()(D dist) {
    CppRandBlockFunct<D> crbf(_seeds, _out, _states, dist);
    BlockApplyParallelFunct<CppRandBlockFunct<D> > bapf(crbf);
    bapf(_n);
  }
      
  std::string *get_states() const { return _states; }
  
private:
  MWArrayContainer<unsigned int> _seeds;
  MWArrayContainer<double> _out;
  std::string *_states;
  size_t _n;
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check arguments
  if (nrhs < 3)
    mexErrMsgIdAndTxt("cpprandpar:nrhs", "Arguments should be M and N dimensions and a rowvector of seeds");
  if (!mxIsNumeric(prhs[0]) || !mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[0]) != 1 || mxGetNumberOfElements(prhs[1]) != 1)
    mexErrMsgIdAndTxt("cpprand:nrhs", "First and second arguments should be scalars giving the M and N dimensions");

  // 3rd argument gives the starting state for RNGs, so its length 
  // determines how many RNGs to initialize, i.e. how many parallel threads
  size_t n = mxGetNumberOfElements(prhs[2]);
  
  // Initialize seeds and states; last arg will specify one or the other
  std::string instates[n];
  const mxArray *seeds;
  
  // 3rd argument can be either a uint32 vector of seed values, or a cell array of state strings
  switch (mxGetClassID(prhs[2])) {
    case mxUINT32_CLASS:
      // Set up with seeds and blank instates
      seeds = prhs[2];
      break;
      
    case mxCELL_CLASS:
      // Set up with blank seeds and filled instates
      seeds = mxCreateNumericMatrix(0, 0, mxUINT32_CLASS, mxREAL);
      for (mwIndex cell = 0; cell < n; ++cell)
        instates[cell] = mxArrayToString(mxGetCell(prhs[2], cell));
      break;
    
    default:
      mexErrMsgIdAndTxt("CPPRANDPAR:prhs", "3rd argument should be a vector of uint32 seeds or a cell array of state strings");
  }
  
  // Create output matrix
  plhs[0] = mxCreateDoubleMatrix(mxGetScalar(prhs[0]), mxGetScalar(prhs[1]), mxREAL);
  MWArrayContainer<double> out(plhs[0]);
  if (n > out.nelem()) mexErrMsgIdAndTxt("CPPRANDPAR:prhs", "You specified more RNGs/threads than output elements???");
  
  std::string *outstates;
  if (nrhs > 3) {
    // Distribution specified as Matlab struct; decode into C++ distribution
    CppRandDistBlockFunct crdbf(MWArrayContainer<unsigned int>(seeds), out, instates, n);
    diststruct_apply<CppRandDistBlockFunct,void>(prhs[3], crdbf);
    outstates = crdbf.get_states();
    
  } else {
    // Use default distribution
    CPPRAND_RANDOM_NAMESPACE::uniform_01<> dist;

    // Setup main logic functor and pass to parallelized block apply
    CppRandBlockFunct<> crbf(MWArrayContainer<unsigned int>(seeds), out, instates, dist);
    BlockApplyParallelFunct<CppRandBlockFunct<> > bapf(crbf);
    bapf(n);
    outstates = crbf.get_states();
  }

  // Done?
  if (nlhs < 2) return;

  // Copy end state string to output
  plhs[1] = mxCreateCellMatrix(n,1);
  for (int i = 0; i < n; ++i) {
    mxSetCell(plhs[1], i, mxCreateString(outstates[i].c_str()));
  }
}
