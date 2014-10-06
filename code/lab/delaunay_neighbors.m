function neighbors_cell = delaunay_neighbors(dt)
% DELAUNAY_NEIGHBORS - Calculate list of neighbors from Delaunay Triangulation DT
%
% usage: neighbors_cell = delaunay_neighbors(dt)
%
% inputs: DT    Delaunay triangulation in either Matlab DelaunayTri, TriRep,
%               or SNL-E DT format.  Or can just be numeric Nx3 of
%               triangles.
%
% output:  Cell array, one element per node.  Each element is a vector
%          indicating the neighbor nodes.
%
% 2010-02 phli
%

if isnumeric(dt)
    tri = dt;
    neighbors_cell = cell(max(tri(:)), 1);
else
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

    neighbors_cell = cell(size(coords, 1), 1);
end

for i = 1:numel(neighbors_cell)
    % Find triangles that include this point
    [row_match, col_match] = find(tri == i);

    neighbors_cell{i} = unique(tri(row_match, :));
    neighbors_cell{i} = setdiff(neighbors_cell{i}, i);
    
    if isfield(dt, 'extra_neighbors')
        neighbors_cell{1} = union(neighbors_cell{i}, dt.extra_neighbors{i});
    end
end