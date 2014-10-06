% BOXSPLINE2_PARALLEL
% usage: val = boxspline2_parallel(x, y)
%
% Parallelized version of BOXSPLINE2.  Relies on Theron Actors library for
% cleaner multithreading.  See http://absurdlycertain.blogspot.com/2011/09/preamble-what-follows-is-guide.html
% for a guide to getting Theron building with mex.
%
% See also: BOXSPLINEN, BOXSPLINE2_PARALLEL
%
% 2011-08, phli
%