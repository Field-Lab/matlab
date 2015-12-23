/* createdist.cpp
 * Ver 0.9.9
 * Peter H. Li 2011 FreeBSD License 
 * See createdist.m for documentation
 */

#include "mex.h"
#include "cpprand_dist2struct.hpp"
#include "cpprand_struct2dist.hpp"

class CreateDistFunct
{
public: 
  CreateDistFunct() {}
  CreateDistFunct(const mxArray *params) : _params(params) {}  
  
  template <typename D>
  mxArray *operator()() {    
    D dist;
    return this->build_dist(dist);
  }

private:
  const mxArray *_params;
  
  /*** build_dist overloaded for every known dist type ***/
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::uniform_01<> dist) {
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::uniform_smallint<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      int min = data[0];
      int max = data[1];
      dist = CPPRAND_RANDOM_NAMESPACE::uniform_smallint<>(min, max);
    }
    return dist_to_struct(dist);
  }

  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::uniform_int_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      int min = data[0];
      int max = data[1];
      dist = CPPRAND_RANDOM_NAMESPACE::uniform_int_distribution<>(min, max);
    }
    return dist_to_struct(dist);
  }

  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::uniform_real_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      double min = data[0];
      double max = data[1];
      dist = CPPRAND_RANDOM_NAMESPACE::uniform_real_distribution<>(min, max);
    }
    return dist_to_struct(dist);
  }

  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::bernoulli_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 0) {
      double *data = mxGetPr(_params);
      double p = data[0];
      if (0 <= p && p <= 1) 
        dist = CPPRAND_RANDOM_NAMESPACE::bernoulli_distribution<>(p);
    }
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::binomial_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      double t = data[0];
      double p = data[1];
      if (t >= 0 && 0 <= p && p <= 1)
        dist = CPPRAND_RANDOM_NAMESPACE::binomial_distribution<>(t, p);
    }
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::geometric_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 0) {
      double *data = mxGetPr(_params);
      double p = data[0];
      if (0 <= p && p <= 1) 
        dist = CPPRAND_RANDOM_NAMESPACE::geometric_distribution<>(p);
    }
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::poisson_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 0) {
      double *data = mxGetPr(_params);
      double mean = data[0];
      if (mean > 0)
        dist = CPPRAND_RANDOM_NAMESPACE::poisson_distribution<>(mean);
    }
    return dist_to_struct(dist);
  }

  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::exponential_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 0) {
      double *data = mxGetPr(_params);
      double lambda = data[0];
      if (lambda > 0)
        dist = CPPRAND_RANDOM_NAMESPACE::exponential_distribution<>(lambda);
    }
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::normal_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      double mean = data[0];
      double sigma = data[1];
      dist = CPPRAND_RANDOM_NAMESPACE::normal_distribution<>(mean, sigma);
    }
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::lognormal_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      double m = data[0];
      double s = data[1];
      dist = CPPRAND_RANDOM_NAMESPACE::lognormal_distribution<>(m, s);
    }
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::cauchy_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      double median = data[0];
      double sigma = data[1];
      dist = CPPRAND_RANDOM_NAMESPACE::cauchy_distribution<>(median, sigma);
    }
    return dist_to_struct(dist);
  }
  
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::triangle_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 2) {
      double *data = mxGetPr(_params);
      double a = data[0];
      double b = data[1];
      double c = data[2];
      dist = CPPRAND_RANDOM_NAMESPACE::triangle_distribution<>(a, b, c);
    }
    return dist_to_struct(dist);
  }

  // Appears that uniform_on_sphere is supposed to return a vector giving 
  // coordinates in dim dimensions, so this doesn't work with current framework.
//   mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::uniform_on_sphere<> dist) {          
//     if (mxIsChar(_params)) return serialized_params_dist(dist);
//     else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 0) {
//       double *data = mxGetPr(_params);
//       double dim = data[0];
//       if (dim >= 0)
//         dist = CPPRAND_RANDOM_NAMESPACE::uniform_on_sphere<>(dim);
//     }
//     return dist_to_struct(dist);
//   }

#ifdef BOOST_RANDOM_GAMMA_DISTRIBUTION_HPP
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::gamma_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);

#ifdef CPPRAND_OLD_GAMMA
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 0) {
      double *data = mxGetPr(_params);
      double alpha = data[0];
      if (alpha > 0)
        dist = CPPRAND_RANDOM_NAMESPACE::gamma_distribution<>(alpha);
#else
    else if (mxIsDouble(_params) && mxGetNumberOfElements(_params) > 1) {
      double *data = mxGetPr(_params);
      double alpha = data[0];
      double beta = data[1];
      if (alpha > 0 && beta > 0)
        dist = CPPRAND_RANDOM_NAMESPACE::gamma_distribution<>(alpha, beta);
#endif

    }
    return dist_to_struct(dist);
  }
#endif  
    
#ifdef CPPRAND_USE_PIECEWISE_LINEAR_DISTRIBUTION
  mxArray *build_dist(CPPRAND_RANDOM_NAMESPACE::piecewise_linear_distribution<> dist) {          
    if (mxIsChar(_params)) return serialized_params_dist(dist);
    else if (mxIsDouble(_params)) {
      MWArrayContainer<double> params(_params);
      MWArrayContainer<double> intervals;
      MWArrayContainer<double> weights;
      if (params.ncols() == 2) {
        intervals = params.get_col(0);
        weights = params.get_col(1);
      } else if (params.nrows() == 2) {
        intervals = params.get_row(0);
        weights = params.get_row(1);
      } else {
        mexErrMsgIdAndTxt("createdist:prhs", "piecewise_linear_distribution should be specified by either a 2xN, or an Mx2 numeric matrix");
      }
      dist = CPPRAND_RANDOM_NAMESPACE::piecewise_linear_distribution<>(intervals.begin(), intervals.end(), weights.begin());
    }
    return dist_to_struct(dist);
  }
#endif

  // General handler for params passed as serialized string (but this probably shouldn't be used in general)
  // Before Boost 1.47, this doesn't work, so just gets the default distribution.
  template <typename D>
  mxArray *serialized_params_dist(D dist) {
#ifdef CPPRAND_PARAMS_SERIALIZABLE
    std::string params_str = mxArrayToString(_params);
    std::stringstream params_stream(std::stringstream::in);
    params_stream.str(params_str);

    typename D::param_type params;
    params_stream >> params;

    dist = D(params);
#endif
    return dist_to_struct(dist);
  }

};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check inputs
  if (nrhs < 1)
    mexErrMsgIdAndTxt("createdist:nrhs", "Argument should be distribution type name, e.g. 'uniform_01'");
  
  CreateDistFunct f;
  if (nrhs > 1) f = CreateDistFunct(prhs[1]);
  plhs[0] = decode_distclass_string_and_call<CreateDistFunct, mxArray*>(mxArrayToString(prhs[0]), f);
}