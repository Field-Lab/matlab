function distances = distance_from_point(matrix_size, point, varargin)
% distance_from_point     at every value in a matrix, identify the distance from a specified point, e.g.
%
%
% > distance_from_point([5 5],[3 3])
%
% ans =
%
%     2.8284    2.2361    2.0000    2.2361    2.8284
%     2.2361    1.4142    1.0000    1.4142    2.2361
%     2.0000    1.0000         0    1.0000    2.0000
%     2.2361    1.4142    1.0000    1.4142    2.2361
%     2.8284    2.2361    2.0000    2.2361    2.8284
%
%
% usage:  distances = distance_from_point(matrix_size, point, varargin)
%
% arguments:  matrix_size - size of the matrix (height, width)
%                   point - x,y coordinates of the point
%                  params - struct or list of optional parameters (see below)
%
% outputs:      distances - distance from the point.  size(distances) = matrix_size
%
%
% optional params, their default values, and what they specify:
%
% radius     	Inf         	radius in which to return the distance
%                                   distance is computed within a square region, sidelength 2*ceil(radius)+1
%                                   points outside the region have the value Inf
%
%
%
% gauthier 2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('radius', Inf, @(x)(x>=0)); % must be non-negative

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% ensure it's an integer
radius = ceil(params.radius);



% BODY OF THE FUNCTION


% identify effective region, i.e. the region in which distances will be computed
eff_x = max(1,round(point(1))-radius) : min(matrix_size(2),round(point(1))+radius);
eff_y = max(1,round(point(2))-radius) : min(matrix_size(1),round(point(2))+radius);
eff_size_x = length(eff_x);
eff_size_y = length(eff_y);

% re-locate center in effective region
eff_ctr = point - [min(eff_x) min(eff_y)] + 1;

% compute distance of each point from the center
x_dist = ([1:eff_size_x] - eff_ctr(1)).^2;
y_dist = ([1:eff_size_y] - eff_ctr(2)).^2;
dist = sqrt(ones(eff_size_y,1)*x_dist + y_dist'*ones(1,eff_size_x));


% make matrix to store result
distances = Inf*ones(matrix_size);

% fill in distances in the desired region
distances(eff_y,eff_x) = dist;

