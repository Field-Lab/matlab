function vwa = variance_weighted_average(estimates, variances)
% variance_weighted_average     combine estimates based on inverse variance weighting
%
%   The estimates are specified in a matrix of any dimension, with the last dimension being an index
%   of the estimates.  Thus this code can compute the variance weighted average for 
%   scalars, vectors, matrices, RGB images, etc.
%
% usage:  vwa = variance_weighted_average(estimates, variances)
%
% arguments:   estimates - D1 x D2 x ... x Dn x N matrix (N = # of estimates)
%              variances - Nx1 vector of the variance of each estimate
%
% outputs:     vwa - D1 x D2 x ... x Dn matrix, the inverse variance weighted average
%
%
%
% 2009-12  gauthier
%


% note size of estimates matrix
orig_size = size(estimates);

% ensure variances is a column vector
variances = reshape(variances,[],1);

% ensure dimensions match
if orig_size(end) ~= length(variances)
    error('Number of estimates (%d) must match number of variances (%d)',orig_size(end),length(variances))
end

% compute weight for each estimate
weights = 1./variances;
weights = weights/sum(weights);

% reshape so each estimate is a column
estimates_reshaped = reshape(estimates,[],orig_size(end));

% compute weighted sum of estimates
vwa = estimates_reshaped*weights;

% reshape
vwa = reshape(vwa,[orig_size(1:end-1) 1]);
