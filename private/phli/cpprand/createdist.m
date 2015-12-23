% CREATEDIST    Brings flexible C++ random distributions into Matlab
% usage: d = createdist('poisson_distribution');
%        d = createdist('poisson_distribution', 3);
%        d = createdist('normal_distribution', [0 1.5]);
%
% then: A = cpprand(10000, 100, uint32(5489), d) % Generate random samples on distribution d
%
% They can all be called without parameters for a default version, or with 
% parameters (usually listed in a numeric vector).
%
% Implemented so far:
%   uniform_01
%   uniform_smallint,           parameters: [min max]
%   uniform_int_distribution,   parameters: [min max]
%   uniform_real_distribution,  parameters: [min max]
%   bernoulli_distribution,     parameter: 0 <= p <= 1
%   binomial_distribution,      parameters: t >= 0, 0 <= p <= 1
%   **geometric_distribution,   parameter: 0 <= p <= 1
%   poisson_distribution,       parameter: mean > 0
%   exponential_distribution,   parameter: lambda > 0
%   **gamma_distribution,       parameters: alpha > 0, beta > 0
%   normal_distribution,        parameters: [mean sigma]
%   **lognormal_distribution,   parameters: [m s]
%   cauchy_distribution,        parameters: [median sigma]
%   triangle_distribution,      parameters: [a b c]
%
%   piecewise_linear_distribution, for this give a 2xN or Mx2 matrix where
%   the first row/column is the values and the second is the probability
%   weights
%
% **Note that gamma distribution has changed parameters between Boost 1.46
% and Boost 1.47; in < 1.47 it takes only an alpha parameter.  Geometric
% and lognormal distributions may also have changed implementation between
% Boost 1.42 and Boost 1.47.
%
% Note that piecewise_linear is only available on Boost >= 1.47.
%
% Note that, on Boost < 1.47 systems, cpprand is generally slower than
% Matlab's rand; even cpprandpar with multithreading is slower.
%
% In general, probabilities do not need to be normalized.
%
%
% Available in Boost Random:
%   http://www.boost.org/doc/libs/1_47_0/doc/html/boost_random/reference.html#boost_random.reference.distributions 
% (If you need one that I haven't wired in yet, shoot me an email).
%

% Version 0.9.9
% Peter H. Li 27-Sep-2011
% As required by MatLab Central FileExchange, licensed under the FreeBSD License