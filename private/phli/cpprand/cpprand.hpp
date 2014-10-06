/* cpprand.hpp
 * Ver 0.9.9
 * Peter H. Li 2011 FreeBSD License 
 * See cpprandpar.m for more documentation.
 */

#ifndef CPPRAND_HPP
#define CPPRAND_HPP

#include <sstream>
#include "mx_functional/mx_functional.h"

// Detect C++11 compiler?  Not working yet... (copied from ../include/boost/config/compiler/gcc.hpp)
// #if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 2)) && defined(__GXX_EXPERIMENTAL_CXX0X__)
// #define CPPRAND_RANDOM_NAMESPACE std

#include "boost/random.hpp"
#include "boost/version.hpp"

// If this is defined, then we used some kind of Boost Random >= 1.47!
#ifdef BOOST_RANDOM_PIECEWISE_LINEAR_DISTRIBUTION_HPP_INCLUDED
#define CPPRAND_RANDOM_NAMESPACE boost::random
#define CPPRAND_PARAMS_SERIALIZABLE
#else
#error You probably should not try to compile this without Boost Random 1.47
#error Try compiling/mexing with the -Iboost_random_1_47 option to kludge Boost Random 1.47 into place.
#error Or you can try taking out these lines and going ahead with your system Boost Random, but that is not recommended
#define CPPRAND_RANDOM_NAMESPACE boost
#define uniform_real_distribution uniform_real
#define uniform_int_distribution uniform_int
#define CPPRAND_OLD_GAMMA
#endif

// Even if we kludged Boost Random >= 1.47, if the real Boost on the system 
// isn't >= 1.47 then we may have problems with piecewise_linear
#if BOOST_VERSION>=104700
#define CPPRAND_USE_PIECEWISE_LINEAR_DISTRIBUTION
#endif


// Functor handling the core logic of generating pseudorandom samples.  Also
// handles on construction various ways of seeding the Mersenne Twister RNG.
// The main methods are to provide a seed, or provide an entire state string.
// 
// In the future I would like to expand this to flexibly handle different 
// RNGs.  It appears that all the generators in Boost Random implement the
// stream operators for serializing state, and they also all seem to be able
// to take some kind of int as a seed; nice!
template <class RNG = CPPRAND_RANDOM_NAMESPACE::mt19937>
class CppRandFunct
{
public:  
  CppRandFunct() {}
  CppRandFunct(RNG rng) : _rng(rng) {}
  CppRandFunct(uint32_t seed) : _rng(seed) {}
  CppRandFunct(std::string rngstatestr) {
    std::stringstream rngstate(std::stringstream::in);
    rngstate.str(rngstatestr);
    rngstate >> _rng;
  }
  
  // Must be const for parallel operation
  template <class D>
  std::string operator()(MWArrayContainer<double> out, D dist) const {
    // Must copy RNG so that operator() can be const
    RNG rng = _rng;
    for (mwSize i = 0; i < out.nelem(); ++i)
      out[i] = dist(rng);
    
    std::stringstream state(std::stringstream::out);
    state << rng;
    return state.str();
  }
  
private:
  RNG _rng;
};


// Simple adapter between StructToDistFunct, which must be used to 
// decode Matlab structs holding C++ distribution classes, and 
// CppRandFunct
template <class RNG = CPPRAND_RANDOM_NAMESPACE::mt19937>
class CppRandDistFunct
{
public:
  CppRandDistFunct(CppRandFunct<RNG> crf, MWArrayContainer<double> out) : _crf(crf), _out(out) {}

  template <typename D>
  std::string operator()(D dist) const {
    return _crf(_out, dist);
  }
  
private:
  CppRandFunct<RNG> _crf;
  MWArrayContainer<double> _out;
};


#endif
