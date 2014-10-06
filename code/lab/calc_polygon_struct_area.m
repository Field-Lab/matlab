function area = calc_polygon_struct_area(polygon_struct)
% CALC_POLYGON_STRUCT_AREA    Calculate the area of a polygon struct
%
% usage: area = calc_polygon_struct_area(polygon_struct)
%
% A polygon struct is the format used by PolygonClip and many SNL-E
% analyses.  Basically fields 'x', 'y', and 'hole'.
%
% 2010-02 phli
%

area = 0;
for i = 1:numel(polygon_struct);
    if polygon_struct(i).hole == 0
        area = area + polyarea(polygon_struct(i).x, polygon_struct(i).y);
    elseif polygon_struct(i).hole == 1
        area = area - polyarea(polygon_struct(i).x, polygon_struct(i).y);
    else
        error('Incorrect argument format');
    end
end