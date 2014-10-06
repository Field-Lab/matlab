/* mx_functional.h
 * Functional library for Matlab MEX functions
 * Ver 0.3.1
 * Peter H. Li 2011 FreeBSD License
 */

#ifndef MX_FUNCTIONAL
#define MX_FUNCTIONAL

#include <mex.h>
#include "mw_array_container.h"


// Functional approach to decoding the underlying type of numeric mxArrays.
// We call this with an mxArray and a functor giving the desired 
// post-decode processing.  Therefore the functor must have templated 
// operator() to be able to handle the different numeric types.
template <typename R, typename F>
R numarr_apply(const mxArray *inarr, F f) {
  switch (mxGetClassID(inarr)) {
    case mxDOUBLE_CLASS:
      return (R) f(MWArrayContainer<double>(inarr));

    case mxSINGLE_CLASS:
      return (R) f(MWArrayContainer<float>(inarr));

    case mxINT8_CLASS:
      return (R) f(MWArrayContainer<char>(inarr));

    case mxUINT8_CLASS:
      return (R) f(MWArrayContainer<unsigned char>(inarr));

    case mxINT16_CLASS:
      return (R) f(MWArrayContainer<signed short>(inarr));

    case mxUINT16_CLASS:
      return (R) f(MWArrayContainer<unsigned short>(inarr));
      
    case mxINT32_CLASS:
      return (R) f(MWArrayContainer<signed int>(inarr));

    case mxUINT32_CLASS:
      return (R) f(MWArrayContainer<unsigned int>(inarr));

    // Uncomment these if int64 is needed, but note that on some compilers
    // it's called "__int64" instead of "long"
    //case mxINT64_CLASS:
      //return (R) f(MWArrayContainer<signed long>(inarr);
      
    //case mxUINT64_CLASS:
      //return (R) f(MWArrayContainer<unsigned long>(inarr);
      
    default:
      mexErrMsgIdAndTxt("Functional:numarr_apply:prhs", "Unrecognized numeric array type.");
  }
}




// Dummy class for template simplification
class MxFunctionalNoOutput {};


// Further abstraction that is commonly useful.  This wraps another Functor
// that accepts input and output MWArrayContainer<D>s.  This can then be 
// passed to numarr_apply to decode the numeric type and the wrapped function 
// will be applied to each column of data and each column of output.
template <typename F, typename R = MxFunctionalNoOutput>
class ColApplyFunct
{
public:
  ColApplyFunct(F f, MWArrayContainer<R> out) : _f(f), _out(out) {}
  
  template <typename D>
  void operator()(MWArrayContainer<D> a) {      
    for (mwIndex col = 0; col < a.ncols(); col++)
      _f(a.get_col(col), _out.get_col(col));
  }
  
private:
  F _f;
  MWArrayContainer<R> _out;
};


// Often, if we don't know /a priori/ the type of the input data, we also don't
// know /a priori/ the type that the output data should be.  This version 
// assumes that what is wanted is for the output data array to be cast to the
// same type as the input data once the input data type has been decoded.
template <typename F>
class ColApplyFunct<F,mxArray *>
{
public:
  ColApplyFunct(F f, const mxArray *out) : _f(f), _out(out) {}

  template <typename D>
  void operator()(MWArrayContainer<D> a) {
      
    // Assume cast output to D is wanted
    MWArrayContainer<D> out(_out);    
    for (mwIndex col = 0; col < a.ncols(); col++)
      _f(a.get_col(col), out.get_col(col));
    
  }

private:
  F _f;
  const mxArray *_out;
};


// Simple version for no output needed; all side-effect methods
template <typename F>
class ColApplyFunct<F,MxFunctionalNoOutput>
{public:
  ColApplyFunct(F f) : _f(f) {}

  template <typename D>
  void operator()(MWArrayContainer<D> a) {
    for (mwIndex col = 0; col < a.ncols(); col++) 
      _f(a.get_col(col));
  }

private:
  F _f;
};



// This version simply calls the passed functor with an iterating block num
// and leaves everything else (usually splitting the operations into blocks)
// to the passed functor.  This is really really pointless; the reason it's 
// here is as a mirror to the BlockApplyParallelFunct class in 
// mx_functional_parallel.  That version is useful because each call of the
// passed functor is executed on a different thread.  So if you can get your
// algorithm working with this stupid stand-in method, you should be able to 
// get out-of-the box multithreading by substituting in BlockApplyParallelFunct
// once you have the necessary parallel libraries (currently Theron Actors
// library).
template <typename F>
class BlockApplyFunct
{
public:
  BlockApplyFunct(F f) : _f(f) {}
  
  void operator()(int nblocks) const {
    for (int blk = 0; blk < nblocks; ++blk) {
      _f(blk, nblocks);
    }
  }
private:
  F _f;
};

#endif
