% CPPRANDPAR    Parallelized wrapper for Mersenne Twister from Boost Random
% usage: R = cpprandpar(M,N,seeds)
%        [R,states] = cpprandpar(M,N,seeds)
%        [R,states] = cpprandpar(M,N,states)
%
% See createdist.m for flexible distribution usage:
%        [R, states] = cpprandpar(M,N,states,dist)
%
% See within cpprandpar.m below for details on the package in general,
% including a little build advice.  CPPRANDPAR is a parallelized wrapper
% for the boost::random (soon to be C++11 std::random) Mersenne Twister
% implementation.  This is significantly faster than Matlab's default
% random number generator (which is also a Mersenne Twister, in all recent
% Matlabs for the last several years), provided you have an up-to-date
% Boost library.
%
% Currently it is only implemented to create 2D matrices given M and N
% dimensions.  It also does not maintain a seed in the background for you,
% as Matlab's implementation does.  So you must seed the RNGs every time
% you call it.  The seeds should be uint32 type (I may add logic to handle
% other types, but it's really better to stick to uint32).  For example, to
% generate a 10x10 random matrix single-threaded you could call:
%   R = cpprandpar(M,N,uint32(5489));
%
% The number of seeds you give is the number of threads that are used; each
% thread has its own independent RNG.  So you should make sure to give
% different seeds.  Two possible approaches are:
%   R = cpprandpar(M,N,uint32(5489:5492))
% to run with 4 threads and defined seeds, or:
%   seeds = randi(intmax, 1, 4, 'uint32');
%   R = cpprandpar(M,N,seeds);
% to have random seeds but still keep track of them (useful if you want to
% be able to rerun later with the same pseudo-random sequence.
%
% If you call this multiple times with the same seed, you will get the same
% sequences.  This is often not what you want.  If you need to keep getting
% different random sequences on subsequent calls, one approach is to
% provide new random seeds on each call with randi:
%   for i = 1:10
%     A = cpprandpar(M,N, randi(intmax, 1, 4, 'uint32'));
%     ... do something with A
%   end
%
% I have also provided a way to get the end state of the RNG back out at
% the end of generating a sequence, and you can feed this state back in to
% the next call to keep the returns changing pseudo-randomly; this is nice
% because you only have to keep track of the initial seed values in order
% to be able to recreate the entire sequence later.  For the Mersenne
% Twister, it is apparently not possible to simply pull a seed back out
% from a used generator, you must instead pull out an entire state vector.
% Boost provides this state as a string:
%   initseeds = randi(intmax, 1, 4, 'uint32'); % 4 seeds, therefore 4 RNGs, therefore 4 threads
%   [~,states] = cpprandpar(0,0, initseeds);   % states is a cell array of strings representing the state of each RNG
%   for i = 1:10
%     [A,states] = cpprandpar(M,N,states);
%     ... do something with A
%   end
%


% BUILD HELP
% I've tested this on Boost 1.42 and Boost 1.47 systems.  Performance is
% much better on 1.47, and in fact things are somewhat unstable on Boost
% 1.42.  Rather than trying to use this on a Boost 1.42 system, it appears
% that better results can be had by kludging some of the 1.47 headers in
% for Random.
%
% To simplify this kludging, I've included a directory, boost_random_1_47.
% If you mex/compile with -Iboost_random_1_47, you will override your
% system's Boost Random with Boost Random 1.47, leaving the rest of your
% Boost intact.  This works surprisingly well, but caused problems for one
% distribution, piecewise_linear_distribution, so I removed this
% distribution unless you have a real Boost >= 1.47 system.  If you want to
% try forcing piecewise_linear back in, you can pass
% -DCPPRAND_USE_PIECEWISE_LINEAR_DISTRIBUTION and see how it goes; YMMV; on
% my machine it builds okay but then running the function fails.
%
% Multithreading (i.e. cpprandpar) relies on the Theron Actors library,
% which in turn relies on either Windows threads or Boost Threads.  Please
% see my blog posts for help getting this built:
%   http://absurdlycertain.blogspot.com/2011/09/preamble-what-follows-is-guide.html
%
% This is fairly straightforward on Mac, but on Linux can be tricky; try
% getting just createdist and cpprand working first before tackling
% cpprandpar.


% CPPRAND PACKAGE
% The CppRand package wraps standard C++ Random Number Generators (RNGs).
% Although it currently relies on the Boost Random library, this
% functionality is included in the new C++11 standard, so as C++11
% compilers become available we should be able to drop the Boost dependency.  
% Therefore I named the package CppRand.
%
% Currently I have only implemented wrappers for the Mersenne Twister
% implementation.  This is the default RNG for recent versions of Matlab.
% Importantly, the Boost 1.42 implementation is actually slower than the
% latest Matlab implementations, but the Boost 1.47 implementation is
% faster.  So this will only be helpful to you if you have a recent Boost 
% (I have not determined exactly when the speedup occurred).  Happily, 
% Boost Random is a header-only library (at least the part used here), 
% so it is not difficult to include Boost Random 1.47 in place of something 
% else that might be on your system.
%
% Benchmarks on a 4-core Linux machine:
%   Boost 1.47:
%       >> tic; for i = 1:10, rand(1000000,10); clear ans; end; toc; clear i
%       Elapsed time is 1.151257 seconds.
%       >> tic; for i = 1:10, cpprand(1000000,10,uint32(5489)); clear ans; end; toc; clear i
%       Elapsed time is 0.729212 seconds.
%       >> tic; for i = 1:10, cpprandpar(1000000,10,uint32(5489:5492)); clear ans; end; toc; clear i
%       Elapsed time is 0.335498 seconds.
%
%   Boost 1.42:
%       >> tic; for i = 1:10, cpprandpar(1000000,10,uint32(5489:5492)); clear ans; end; toc; clear i
%       Elapsed time is 1.405573 seconds.
%   SLOWER than Matab rand even with 4 threads!  Make sure you have up to
%   date Boost.
%


% CHANGELOG (for whole package):
%   0.9.9 - Accepting that it just doesn't work well with Boost 1.42,
%           trying to provide a way to kludge Boost Random 1.47 in instead
%   0.9.8 - Added more distributions; uniform_on_sphere excluded
%   0.9.7 - fixes to get to at least build on Boost < 1.47 systems (still 
%           slower on these systems though, so not much point)
%   0.9.6 - fixing/adding documentation
%   0.9.5 - added major new functionality for flexible distributions
%   0.8.0 - initial release (previously called 0.80)


% Version 0.9.9
% Peter H. Li 27-Sep-2011
% As required by MatLab Central FileExchange, licensed under the FreeBSD License