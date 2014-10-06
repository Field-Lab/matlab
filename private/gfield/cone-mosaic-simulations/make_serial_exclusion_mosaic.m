function locations = make_serial_exclusion_mosaic(grid_size, num_points, exclusion_mean, exclusion_sigma, varargin)
% make_serial_exclusion_mosaic     the function returns a set of mosaic
%                                  locations produced by an exclusion zone algorithm for generating 2-D
%                                  spatial mosaics.  The exclusion zone is
%                                  produced by a Guassian distribution
%                                  where for each point a sample is drawn
%                                  from the distribution to determine a
%                                  point falling close to an old point is
%                                  kept or discarded.
%
% usage:  locations = make_serial_exclusion_mosaic(grid_size, num_points,...
%                                                           exclusion_mean, exclusion_sigma, varargin)
%
% arguments:     grid_size - 2D vector specifying the spatial extent of mosaic
%               num_points - number of points in mosaic
%           exclusion_mean - mean of exclusion gaussian
%          exclusion_sigma - sigma of exclusion gaussian
%
% outputs:       locations - locations of points in a mosiac
%
%
%
% GDF 2009-2-17
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addRequired('grid_size', @isnumeric);
p.addRequired('num_points', @isnumeric);
p.addRequired('exclusion_mean', @isnumeric);
p.addRequired('exclusion_sigma', @isnumeric);
p.addParamValue('verbose', false, @islogical)
p.addParamValue('seed', [], @isnumeric)

% resolve user input and default values
p.parse(grid_size, num_points, exclusion_mean, exclusion_sigma, varargin{:});

verbose = p.Results.verbose;
seed = p.Results.seed;

% initialize random number generator
if ~isempty(seed)
    rand('twister',seed);
end

if verbose
    h = waitbar(0, 'generating mosaic locations');
end

temp_locations = zeros(num_points,2);
samples = 1;
% generate intial, unconstrained location
temp_locations(samples, :) = grid_size .* rand(1,2);
% generate the remaining locations
while samples <= num_points
    if verbose 
        waitbar(samples/num_points, h)
    end
    test_location = grid_size .* rand(1,2);
    d_min = normrnd(exclusion_mean, exclusion_sigma);
    distances = pdist([test_location; temp_locations]);
    distances = distances(1:samples);
    if isempty(find(distances < d_min, 1))
       temp_locations(samples,:) = test_location;
       samples = samples + 1;
    end 
end

if verbose
    close(h)
end

locations = temp_locations;
    
    