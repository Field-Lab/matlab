function d = distance_point_to_segment(point,endpoint_1,endpoint_2)
% distance_point_to_segment     compute distance between a point and a line segment
%
% usage: d = distance_point_to_segment(point,endpoint_1,endpoint_2)
%
%
% 2010-03  gauthier
%


% ensure row vectors
if size(point,1) > 1; point = point'; end
if size(endpoint_1,1) > 1; endpoint_1 = endpoint_1'; end
if size(endpoint_2,1) > 1; endpoint_2 = endpoint_2'; end


% ensure 2D vectors
if numel(point) > 2 || numel(endpoint_1) > 2 || numel(endpoint_2) > 2
    error('inputs must be 2 element vectors')
end

% re-assign names
P = point;
Q1 = endpoint_1;
Q2 = endpoint_2;

% compute distance
d = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
