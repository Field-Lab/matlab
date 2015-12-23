function [x,y,indices] = rf_cone_profile(cone_weights,cone_centers, center, varargin)
% rf_cone_profile     RF cone weights as a function of distance from the center
%
% usage:  [x,y,cone_types] = rf_profile(cone_weights,cone_centers,center,varargin)
%
% arguments:       
%        cone_weights - N-length vector of cone weights
%        cone_centers - Nx2 matrix, x,y coordinates of each cone center point
%              center - x,y coordinates of the RF center
%            varargin - struct or list of optional parameters (see below)
%
% outputs:          x - N-length vector of distances from center
%                   y - N-length vector of cone strengths at the points in X
%             indices - N-length vector of cone indices
%
%
% optional params, their default values, and what they specify:
%
% radius            Inf                 only return values within the specified radius
%
%
% gauthier 2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('radius', Inf, @(x)(x>=0));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

num_cones = length(cone_weights);

% get distances from center point
distances = sqrt(sum(    ( cone_centers - repmat(center,num_cones,1) ).^2   ,2));

% identify cones within the specified radius
near_cones = find(distances < params.radius);

% get distances at these points
x = distances(near_cones);

% get cone weights at these points
y = cone_weights(near_cones);

% sort the points
[x,sort_order] = sort(x);
y = y(sort_order,:);
indices = near_cones(sort_order);
