function cones = get_cones_by_weight_rank(datarun, rgc_id, ranks, varargin)
% 2011-07 phli

% Get weights
[mosaic_weights, selection, extras] = select_cone_weights(datarun, rgc_id, varargin{:});
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% Rank cones by weight
[~,sorted_cone_indices] = sort(connectivity, 1, 'descend');
cones = sorted_cone_indices(ranks);