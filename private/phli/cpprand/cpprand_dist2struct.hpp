/* cpprand_dist2struct.hpp
 * Version 0.3.3
 * Peter H. Li 2011 FreeBSD License
 * Converting from generic distribution classes to Matlab struct arrays
 * 
 * C++ probability distribution classes allow a single RNG to flexibly generate
 * samples from many underlying distributions, including the user-definable 
 * piecewise_linear_distribution and discrete_distribution.  These classes are
 * required to override the stream operators (<< and >>) to implement string
 * based serialization of their internal state.  Thus it is not too painful to
 * bring the flexibility of distribution classes into Matlab by allowing the
 * mex file createdist.cpp to build a C++ ditstribution object and return it
 * to Matlab in its serialized state, and then allow Matlab to pass the 
 * serialized string back into C++ to get random samples.
 *
 * There was of course still a little futzing required to get this working.  
 * I had to decide on a Matlab data structure to store the serialized 
 * distribution objects.  It wasn't easy to return a Matlab OO object from 
 * C++, so in the end I decided to stick with a simple struct array.
 * 
 * The logic here is used primarily to convert from C++ distributions to
 * a Matlab struct, while cpprand_struct2dist deals with the reverse
 * conversion.
 */

#ifndef CPPRAND_DIST2STRUCT_HPP
#define CPPRAND_DIST2STRUCT_FPP

#include <sstream>
#include "cpprand.hpp"


// Essentially just a map from distribution types to strings.  There are
// various ways to implement this but this seemed cleanest in the end and
// there are not so many types as to make this a total pain.
template <typename D> std::string get_disttype();
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::uniform_01<>                    >() { return "uniform_01"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::uniform_smallint<>              >() { return "uniform_smallint"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::uniform_int_distribution<>      >() { return "uniform_int_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::uniform_real_distribution<>     >() { return "uniform_real_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::bernoulli_distribution<>        >() { return "bernoulli_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::binomial_distribution<>         >() { return "binomial_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::geometric_distribution<>        >() { return "geometric_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::poisson_distribution<>          >() { return "poisson_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::exponential_distribution<>      >() { return "exponential_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::normal_distribution<>           >() { return "normal_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::lognormal_distribution<>        >() { return "lognormal_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::cauchy_distribution<>           >() { return "cauchy_distribution"; }
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::triangle_distribution<>         >() { return "triangle_distribution"; }
// template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::uniform_on_sphere<>             >() { return "uniform_on_sphere"; }

#ifdef BOOST_RANDOM_GAMMA_DISTRIBUTION_HPP
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::gamma_distribution<>            >() { return "gamma_distribution"; }
#endif

#ifdef CPPRAND_USE_PIECEWISE_LINEAR_DISTRIBUTION
template <> std::string get_disttype<CPPRAND_RANDOM_NAMESPACE::piecewise_linear_distribution<> >() { return "piecewise_linear_distribution"; }
#endif

// Convert from the distribution D to a Matlab struct.  The struct
// identifies itself as a CppRandDist, in place of real OO, and also
// stores the specific C++ distribution type (e.g. "uniform_01") and
// of course the state string that represents the distribution.
template <typename D>
mxArray *dist_to_struct(D dist) {
  mxArray *classname = mxCreateString("CppRandDist");
  mxArray *type = mxCreateString(get_disttype<D>().c_str());

  std::string randversionstr = "Boost ";
  randversionstr += BOOST_LIB_VERSION;
  mxArray *randversion = mxCreateString(randversionstr.c_str());

  std::stringstream s(std::stringstream::out);
  s << dist;
  mxArray *state = mxCreateString(s.str().c_str());

  const char *fieldnames[] = {"class", "randversion", "type", "state"};
  mxArray *diststruct = mxCreateStructMatrix(1, 1, 4, fieldnames);
  mxSetField(diststruct, 0, "class", classname);
  mxSetField(diststruct, 0, "randversion", randversion);
  mxSetField(diststruct, 0, "type", type);
  mxSetField(diststruct, 0, "state", state);

  return diststruct;
}


#endif
