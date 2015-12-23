function simple = simplify_contour(contour, num_keep)
% SIMPLIFY_CONTOUR    Reduce contour polygon with subset of largest area polygons
% usage: simple = simplify_contour(contour, num_keep)
%
% phli 2010-03
%

if nargin < 2
    num_keep = 1;
end


areas = zeros(num_keep, 1);

for i = 1:numel(contour)
    polygon = contour(i);
    if polygon.hole == 1
        continue;
    end
    
    area = calc_polygon_struct_area(polygon);
    [min_area, min_index] = min(areas);
    if area > min_area
        % Replace smallest polygon with current polygon
        min_index = min_index(1);
        areas(min_index) = area;
        simple(min_index) = polygon;
    end
end