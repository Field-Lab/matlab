function polygons = triangles_to_polygons(triples, coordinates)
% TRIANGLES_TO_POLYGONS    Convert triangles in triple/coord format into polygon structs
% 
% usage: polygons = triangles_to_polygons(triples, coordinates)
%
% Triangles form Matlab delaunay triangulation methods typically come out
% in triple/coord format.  Here, we convert this to the polygon struct
% format used by PolygonClip and often used by SNL-E.
%
% 2010-01 phli
%

% This could be more efficient, but doesn't really matter.

polygons = struct;
for i = 1:size(triples, 1);
    polygons(i).x = coordinates(triples(i,:), 1);
    polygons(i).y = coordinates(triples(i,:), 2);
    polygons(i).hole = 0;
end