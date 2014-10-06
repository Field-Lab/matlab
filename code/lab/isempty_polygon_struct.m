function bool = isempty_polygon_struct(pgs)
% ISEMPTY_POLYGON_STRUCT    Check whether an SNL-E standard polygon struct actually defines any polygon
% usage: bool = isempty_polygon_struct(pgs)
%
% 2010-04 phli
%

bool = isempty(pgs) || (length(pgs) == 1 && length(pgs.x) < 1);