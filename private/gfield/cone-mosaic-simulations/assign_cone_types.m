function types = assign_cone_types(locations, l_proportion, varargin)
% assign_cone_types     assign the cone types in a mosaic
%
% usage:  assign_cone_types(locations, varargin)
%
% arguments:  locations - an N by 2 matrix of cone locations
%
% outputs:     types - a vector of length N assigning types to the cones
%
%
% optional fields in params, their default values, and what they specify:
%
%  seed                     []              random number seed to make the "random" sampling reproducible
%  clumping                false            if true spatial filtering is used to clump cones
%  filter_sigma            1.0              clumping parameter: spatial scale of clumping
%  delta_pL                0.75             clumping parameter: clumping strength
%  clumping_seed            []              clumping parameter: seeds the random pattern that is filtered
%
% 2009-03 gdf
%

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addRequired('locations', @isnumeric);
p.addRequired('l_proportion', @isnumeric);
p.addParamValue('clumping', false, @islogical);
p.addParamValue('filter_sigma', 1.0, @isnumeric);
p.addParamValue('delta_pL', 0.75, @isnumeric);
p.addParamValue('seed', [], @isnumeric);

p.parse(locations, l_proportion, varargin{:});

seed = p.Results.seed;

clumping = p.Results.clumping;

% check to see if seed is specified

if ~clumping
    if ~isempty(seed)
        rand('twister', seed);
    end
    temp_cone_types = rand(length(locations(:,1)),1);
    l_cone_ids = find(temp_cone_types <= l_proportion);
    m_cone_ids = find(temp_cone_types >= l_proportion);
    cone_types(l_cone_ids) = 'L';
    cone_types(m_cone_ids) = 'M';
    types = cone_types;
end

% clumped cone mosaic
if clumping
    filter_sigma = p.Results.filter_sigma;
    delta_pL = p.Results.delta_pL;
    [cone_types, extras] = filter_to_clump(locations, filter_sigma, l_proportion, delta_pL, 'seed', seed);
    types = cone_types;
end

