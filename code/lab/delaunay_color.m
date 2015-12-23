function coloring = delaunay_color(dt, init_coloring)
% DELAUNAY_COLOR Do a greedy topocoloring based on Delaunay Triangulation DT
%
% usage: coloring = delaunay_color(dt)
%
% inputs:   DT - Either a DelaunayTri MatLab object, or the SNL-E custom struct format (dt.triangles, dt.coordinates, etc.)
%           INIT_COLORING - Initial filling, perhaps from mapped master datarun; we'll only fill in where init is empty
%
% outputs:  coloring: An Nx1 array of integers each representing a color, where N = the number of coordinate points in DT
%
% This is a standard algorithm for coloring a topological collection of
% regions so that no two regions that border (as determined by DT/VT) will 
% be the same color.  Greedy because we just go through the list of regions
% and assign each the first available color based on its neighbors that
% have already been colored.  This does not generally get an optimal
% coloring (fewest colors possible), but it's good enough.
%
% 2010-02 phli
%

switch class(dt)
    case 'struct'
        % SNL-E custom DT format
        tri = dt.triangles;
        coords = dt.coordinates;
    case 'DelaunayTri'
        tri = dt.Triangulation;
        coords = dt.X;
    case 'TriRep'
        tri = dt.Triangulation;
        coords = dt.X;
    otherwise
        error('Unrecognized triangulation format');
end


if nargin < 2
    init_coloring = zeros(size(coords, 1), 1);
end


neighbors = delaunay_neighbors(dt);
coloring = init_coloring;
for i = 1:numel(coloring)
    % If init_coloring already filled this one in, skip ahead
    if coloring(i)
        continue;
    end
        
    neighbor_colors = coloring(neighbors{i});
    
	j = 1;
	while true
		if ~any(neighbor_colors == j)
			coloring(i) = j;
			break;
		end
		j = j + 1;
	end
end