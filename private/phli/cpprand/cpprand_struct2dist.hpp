/* cpprand_struct2dist.hpp
 * Version 0.3.3
 * Peter H. Li 2011 FreeBSD License
 * Converting from Matlab struct arrays to generic distribution classes
 * 
 * See cpprand_dist2struct.hpp for general background on the design.
 * Basically these two headers allow us to bring flexible C++ probability
 * distributions into our implementation of C++ RNGs in Matlab.
 *
 * This file deals with converting the Matlab structs back into various
 * possible C++ distribution types.  This is the more difficult half of 
 * the equation, as we are going from a single Matlab struct array type 
 * to many possible C++ distribution types.  
 *
 * There are various possible designs to handle this, but I've found 
 * the functor approach is simpler than inheritance based methods.  
 * Basically there is a single method, diststruct_apply, that you can 
 * call on a Matlab struct representing the C++ distribution.  This 
 * converts the Matlab struct into the appropriate C++ distribution and 
 * then applies a passed functor.  The passed functor should be 
 * templated to handle every possible C++ distribution type.  As long as 
 * you can wrap all the logic you want to operate on the C++ 
 * distribution inside a functor, this approach allows you to avoid the
 * issue of coming up with a way to store the decoded Matlab struct that can
 * be any of a number of classes.
 */

#ifndef CPPRAND_STRUCT2DIST_HPP
#define CPPRAND_STRUCT2DIST_HPP

#include <string.h>
#include "cpprand.hpp"


// Just decodes distclass string and calls Functor argument with templated 
// operator() with correct class.  This is used by createdist to parse the
// dist type string originally passed by the user.  It's also used to decode
// Matlab diststructs, although that functionality is more nicely wrapped in
// diststruct_apply.
template <typename F, typename R>
R decode_distclass_string_and_call(const char *distclass, F f) {
  if (strcmp(distclass, "uniform_01") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::uniform_01<> >();
  else if (strcmp(distclass, "uniform_smallint") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::uniform_smallint<> >();
  else if (strcmp(distclass, "uniform_int_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::uniform_int_distribution<> >();
  else if (strcmp(distclass, "uniform_real_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::uniform_real_distribution<> >();

  else if (strcmp(distclass, "bernoulli_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::bernoulli_distribution<> >();
  else if (strcmp(distclass, "binomial_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::binomial_distribution<> >();
  else if (strcmp(distclass, "geometric_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::geometric_distribution<> >();
  
  else if (strcmp(distclass, "poisson_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::poisson_distribution<> >();
  else if (strcmp(distclass, "exponential_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::exponential_distribution<> >();
  
  else if (strcmp(distclass, "normal_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::normal_distribution<> >();
  else if (strcmp(distclass, "lognormal_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::lognormal_distribution<> >();
  else if (strcmp(distclass, "cauchy_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::cauchy_distribution<> >();
  
  else if (strcmp(distclass, "triangle_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::triangle_distribution<> >();
//   else if (strcmp(distclass, "uniform_on_sphere") == 0)
//     return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::uniform_on_sphere<> >();

#ifdef BOOST_RANDOM_GAMMA_DISTRIBUTION_HPP
  else if (strcmp(distclass, "gamma_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::gamma_distribution<> >();
#endif
  
#ifdef CPPRAND_USE_PIECEWISE_LINEAR_DISTRIBUTION
  else if (strcmp(distclass, "piecewise_linear_distribution") == 0)
    return (R) f.template operator()<CPPRAND_RANDOM_NAMESPACE::piecewise_linear_distribution<> >();
#endif
          
  else
    mexErrMsgIdAndTxt("cpprand:decode_dist_string_and_call", "Unrecognized distribution class");
}


// Convert diststruct into the actual distribution and apply another functor;
// assumes diststruct type has already been decoded; i.e. this should only be
// called from decode_distclass_string_and_call
template <typename F, typename R>
class DistStructApplyFunct
{
public:
  DistStructApplyFunct(const mxArray *diststruct, F f) : _diststruct(diststruct), _f(f) {}
  
  template <typename D>
  R operator()() {
    D dist;

    mxArray *statefield = mxGetField(_diststruct, 0, "state");
    if (statefield == NULL) return (R) _f(dist);

    const char *state_cstr = mxArrayToString(statefield);
    if (state_cstr == NULL) return (R) _f(dist);
    
    std::string statestr = state_cstr;
    std::stringstream statestream(std::stringstream::in);
    statestream.str(statestr);
    statestream >> dist;
    
    return (R) _f(dist);
  }

private:
  F _f;
  const mxArray *_diststruct;
};


// Assuming you have a Matlab struct representing a C++ distribution, this
// nicely wraps everything for you so you can put your logic in a templated
// functor, pass the diststruct and functor into this, and get the correct
// template call from your functor running.
//
// Contract is that this will check that diststruct is a struct array and
// is not empty.  Will also check that a 'class' field exists set to 
// 'CppRandDist' and a 'type' field exists and has a string value.  Via
// underlying functions, it then checks that the type is a valid
// distribution and that any 'state' field is valid.
template <typename F, typename R>
R diststruct_apply(const mxArray *diststruct, F f) {
  if (!mxIsStruct(diststruct) || mxGetNumberOfElements(diststruct) < 1)
    mexErrMsgIdAndTxt("cpprand:diststruct_apply", "Argument must be a non-empty struct array");

  mxArray *classfield = mxGetField(diststruct, 0, "class");
  if (classfield == NULL) 
    mexErrMsgIdAndTxt("cpprand:diststruct_apply", "Struct array must have a 'class' field specifying CppRandDist");
  const char *classname = mxArrayToString(classfield);
  if (classname == NULL || strcmp(classname, "CppRandDist") != 0)
    mexErrMsgIdAndTxt("cpprand:diststruct_apply", "'class' must be 'CppRandDist'");

  mxArray *disttypefield = mxGetField(diststruct, 0, "type");
  if (disttypefield == NULL) 
    mexErrMsgIdAndTxt("cpprand:diststruct_apply", "Distribution struct must have a 'type' field specifying a C++ distribution type");
  const char *disttype = mxArrayToString(disttypefield);
  if (disttype == NULL)
    mexErrMsgIdAndTxt("cpprand:diststruct_apply", "Failed to read string from 'type' field");

  DistStructApplyFunct<F,R> dsaf(diststruct, f);
  return (R) decode_distclass_string_and_call<DistStructApplyFunct<F,R>,R>(disttype, dsaf);
}

#endif
