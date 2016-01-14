function output_matrix = make_Gaussian_two_d(varargin)
% diff_Gaussians_two_d(varargin)
%
%   This function makes a 2 dimensional Gaussian in a field set by x_dim
%   and y_dim.
%
% Author: GDF 
% Date: 2011-06-11


% DEFINE INPUTS AND PARSE

p = inputParser;
p.addParamValue('center_point_x', 3, @isnumeric);
p.addParamValue('center_point_y', 3, @isnumeric);
p.addParamValue('sd_x', 1, @isnumeric);
p.addParamValue('sd_y', 1, @isnumeric);
p.addParamValue('amp_scale', 1, @isnumeric);
p.addParamValue('rotation_angle', 0, @isnumeric);
p.addParamValue('x_dim', 5, @isnumeric);
p.addParamValue('y_dim', 5, @isnumeric);

p.parse(varargin{:});

% assign the parsed results to simpler (original) names
center_point = [p.Results.center_point_x, p.Results.center_point_y];
sd_scale = [p.Results.sd_x, p.Results.sd_y];
amp_scale = p.Results.amp_scale;
rotation_angle = p.Results.rotation_angle;
x_dim = p.Results.x_dim;
y_dim = p.Results.y_dim;


% BODY OF FUNCTION
% make the Gaussian according the given parameters

% initialize the output matrix
output_matrix = zeros(y_dim, x_dim);

% make an array of points for matrix (STA) values
width_points = 1:1:x_dim;
height_points = 1:1:y_dim;

% calculate the distances of these points from the center of Gauss
width_dists = center_point(1) - width_points;
height_dists = center_point(2) - height_points;

% calculate rotation matrix: couterclockwise rotation with respect to angle
rotation_matrix = [cos(rotation_angle), -1*sin(rotation_angle); sin(rotation_angle), cos(rotation_angle)];

% define covariance matrix given the sd_scale and rotation matrix
covariance_matrix = rotation_matrix * [1/sd_scale(1)^2 0; 0 1/sd_scale(2)^2] * rotation_matrix';

% calculate the value of the Gaussian at each point in output_matrix
for wd = 1:x_dim
    for ht = 1:y_dim
        pt = [height_dists(ht); width_dists(wd)];
        output_matrix(ht,wd) = amp_scale .* exp(-0.5 .* (pt' * covariance_matrix * pt));
    end
end
