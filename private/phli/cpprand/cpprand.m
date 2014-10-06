% CPPRAND    Wrapper for Mersenne Twister from Boost Random
% usage: R = cpprand(M,N,seed)
%        [R,state] = cpprandpar(M,N,seed)
%        [R,state] = cpprandpar(M,N,state)
%
% See createdist.m for flexible distribution usage:
%        [R, states] = cpprandpar(M,N,states,dist)
%
% See CPPRANDPAR for more details on implementation of the package.
%
% CPPRAND is a wrapper for the boost::random (soon to be C++11 std::random)
% Mersenne Twister implementation.  This is ~2x faster than Matlab's
% default implementation (which is also a Mersenne Twister, in all
% recent Matlabs for the last several years), provided you have an
% up-to-date Boost library.
%
% Currently it is only implemented to create 2D matrices given M and N
% dimensions.  It also does not maintain a seed in the background for you,
% as Matlab's implementation does.  So you must seed the RNG every time
% you call it.  The seeds should be uint32 type (I may add logic to handle
% other types, but it's really better to stick to uint32).  For example, to
% generate a 10x10 random matrix you could call:
%   R = cpprandpar(M,N,uint32(5489));
%
% If you call this multiple times with the same seed, you will get the same
% sequence.  This is often not what you want.  If you need to keep getting
% different random sequences on subsequent calls, one approach is to
% provide a new random seed on each call with randi:
%   for i = 1:10
%     A = cpprand(M,N, randi(intmax, 'uint32'));
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
%   initseed = randi(intmax, 'uint32'); % Save this seed value to recreate sequence later
%   [~,state] = cpprand(0,0, initseed); % state is a char array representing the ending state of the RNG
%   for i = 1:10
%     [A,state] = cpprand(M,N,state);
%     ... do something with A
%   end

% Version 0.9.9
% Peter H. Li 27-Sep-2011
% As required by MatLab Central FileExchange, licensed under the FreeBSD License