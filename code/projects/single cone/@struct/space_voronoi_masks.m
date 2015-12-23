function spaced = space_voronoi_masks(datarun, neighbor_min_dist, self_max_dist)
% SPACE_VORONOI_MASKS   Cut Voronoi mask pixels that impinge too closely on neighboring cone centers or are too far from their own cone centers
% usage: spaced = space_voronoi_masks(datarun, neighbor_min_dist, self_max_dist)
% 
% inputs: DATARUN (should run MAKE_VORONOI_MASKS first)
%         NEIGHBOR_MIN_DIST      If empty, set to 0
%         SELF_MAX_DIST          If empty, set to Inf
%
% outputs: SPACED masks cell array of length equal to the number of cones,
%          where each element is logical pixel array.
%
% See also: MAKE_VORONOI_MASKS
%
% phli 2011-01
%

if nargin < 2 || isempty(neighbor_min_dist)
    neighbor_min_dist = 0;
end

if nargin < 3 || isempty(self_max_dist)
    self_max_dist = Inf;
end

masks = datarun.cones.mosaic.voronoi_masks;
scale = datarun.cones.mosaic.voronoi_masks_scale;
centers = datarun.cones.centers .* scale;
near_cones_list = delaunay_neighbors(datarun.cones.mosaic.delaunay);

spaced = space_voronoi_masks(masks, centers, near_cones_list, neighbor_min_dist, self_max_dist);