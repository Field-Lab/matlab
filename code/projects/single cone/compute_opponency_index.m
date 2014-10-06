function opponency_index = compute_opponency_index(weights, labels, varargin)
%
% compute_opponency_index   Computes the L versus M cone opponency from a weight matrix
%
% usage: oppenency_index = compute_opponency_index(weights, labels)
%
% arguments:         weights - a matrix [P x Q] of weights for P cones inputs to Q RGCs
%                     labels - character ('L' and 'M') labels to cones. S is ignored
%
% outputs:   opponency_index - a vector of opponency indices for the RGCs
%
% optional params, their default values, and what they specify:
%
% none at this time
%
% GDF 2009-02
%

% initialize and check inputs
p = inputParser;

p.addRequired('weights', @isnumeric)
p.addRequired('labels', @ischar)
p.parse(weights, labels, varargin{:})

% check inputs
[num_cones, num_RGCs] = size(weights);
if length(labels) ~= num_cones
    error('Input dimension mismatch: number of cone labels must equal number of cones')
end

% body of function

% get indices to l and m cones
L_indices = find(labels == 'L');
M_indices = find(labels == 'M');

L_weights = weights(L_indices,:);
M_weights = weights(M_indices,:);

opponency_index = (sum(L_weights,1) - sum(M_weights,1)) ./ (sum(abs(L_weights),1) + sum(abs(M_weights),1));
