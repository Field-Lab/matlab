function locations = hexagonal_mosaic(points, varargin)
%
% usage: locations = hexagonal_mosaic(points, varargin)
%
% description: generates a hexagonal array of locations
%
% required input    points        number of points in the hexagonal lattice
%                                 this value must satisfy sqrt(points) = integer
%
% optional inputs        default              Explanation
%   jitter_method    'bounded_uniform'        Selects the method for jittering mosiacs away from a 
%                                             perfectly hexagonal lattice.
%                                             Options:
%                                             - bounded_uniform
%                                                 BU_scale (default = 0.3) sets the width of the jitter
%                                             - Gaussian
%                                                 Gaussian_SD (default = 0.25) sets the width of the jitter
%                                             - exponential
%                                                 exp_mean (default = 0.5) sets the width of the jitter
%
%   spacing               1.0                  Distance between points in multiples of sqrt(3)/2
%   seed                   []                  The random jitter can be seeded so that
%                                              it is fixed from one simulation to the next
%                   
%
%   GDF 2009-01

p = inputParser;
p.addRequired('points', @(x)mod(sqrt(x),1) == 0)

p.addParamValue('spacing', 1.0, @isnumeric)
p.addParamValue('seed', [])
p.addParamValue('jitter_method','bounded_uniform', @ischar)
p.addParamValue('Gaussian_sd', 0.0, @isnumeric)
p.addParamValue('BU_scale', 0.25, @isnumeric)
p.addParamValue('exp_mean', 0.5, @isnumeric)

p.parse(points, varargin{:})

seed = p.Results.seed;
Gaussian_sd = p.Results.Gaussian_sd;
spacing = p.Results.spacing * sqrt(3)/2;
jitter_method = p.Results.jitter_method;
BU_scale = p.Results.BU_scale;
exp_mean = p.Results.exp_mean;

% construct hexagonal lattice
[x y] = meshgrid(0:1:(sqrt(points)-1));
n = size(x,1);
X = x * spacing;
Y = (y*spacing) + repmat([0 0.5], [n,n/2]);

temp_X = reshape(X, [], 1);
temp_Y = reshape(Y, [], 1);

% summary of cone locations and type
locations = [temp_X, temp_Y];

jitter = zeros(points, 2);  % initialize jitter

% make jittered noise
if strcmp(jitter_method, 'Gaussian')
    if ~isempty(seed)
        rand('twister', seed)
    end
    jitter = mvnrnd([0 0], diag([Gaussian_sd Gaussian_sd]), points);
    
end

if strcmp(jitter_method, 'bounded_uniform')
    if ~isempty(seed)
        rand('twister', seed)
    end
    jitter = BU_scale .* ((rand(points,2) .* 2) -1);
    
end
    
if strcmp(jitter_method, 'exponential')
    if ~isempty(seed)
        rand('twister', seed)
    end
    x_jitter = exprnd(exp_mean, points, 1);
    y_jitter = exprnd(exp_mean, points, 1);
    x_signs = rand(points,1);
    y_signs = rand(points,1);
    x_pos_index = find(x_signs > 0.5);
    x_neg_index = find(x_signs <= 0.5);
    y_pos_index = find(y_signs > 0.5);
    y_neg_index = find(y_signs <= 0.5);
    x_signs(x_pos_index) = 1;
    x_signs(x_neg_index) = -1;
    y_signs(y_pos_index) = 1;
    y_signs(y_neg_index) = -1;
    
    x_jitter = x_jitter .* x_signs;
    y_jitter = y_jitter .* y_signs;
    jitter = [x_jitter, y_jitter];
        
    size(jitter);
end
    

locations = locations + jitter;

