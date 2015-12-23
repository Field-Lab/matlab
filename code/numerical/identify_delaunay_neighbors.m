function neighbors = identify_delaunay_neighbors(centers, scale)
% MY_FUNCTION     This template function does nothing.  Ha!
%
% usage:  result = my_function(arg1, <params>)
%
% arguments:     arg1 - argument 1
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
%
% 2008-10  gauthier
%



% get NND

% get all pairwise distances in matrix
dists = squareform(pdist(centers));

% put large values on diagonal
dists = dists + realmax*eye(size(dists));

% get median nearest neighbor distance
nnd = median(min(dists));

neighbors = nnd;
