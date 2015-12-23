function [x,y] = rf_profile(rf, center, varargin)
% rf_profile     RF strength as a function of distance from the center
%
% usage:  [x,y] = rf_profile(rf, center, varargin)
%
% arguments:       rf - 2D or 3D RF
%              center - x,y coordinates of the center
%            varargin - struct or list of optional parameters (see below)
%
% outputs:          x - N-length vector of distances from center
%                   y - NxC matrix of rf strength at the points in X.  C = size(rf,3)
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

% get distances from center point
distances = distance_from_point([size(rf,1) size(rf,2)], center, 'radius',params.radius);

% identify region in which distances are less than the specified radius
near_stixels = (distances < params.radius);

% get distances at these points
x = distances(near_stixels);

% get RF values at these points
y = zeros(length(x),size(rf,3));
for cc = 1:size(rf,3)
    rf_temp = rf(:,:,cc);
    y(:,cc) = rf_temp(near_stixels);
end

% sort the points
[x,sort_order] = sort(x);
y = y(sort_order,:);
