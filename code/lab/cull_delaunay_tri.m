function datarun = cull_delaunay_tri(datarun, cell_spec, distance_factor, varargin)
% CULL_DELAUNAY_TRI    Cut out triangles with excessively long edges
%
% usage: datarun = cull_delaunay_tri(datarun, cell_spec, distance_factor, opts)
%
% opts: neighbor_criterion      'nearest nbr median'
%       plot                    false
%
% Pulled over from Jeff's interdigitation work.
%
% 2010-02 phli
%

opts = inputParser;
opts.addParamValue('neighbor_criterion', 'nearest nbr median');
opts.addParamValue('plot', false);
opts.parse(varargin{:});
opts = opts.Results;


cell_type_num = get_cell_type_nums(datarun, cell_spec);
delaunay = datarun.stas.delaunays{cell_type_num};

% Check whether this has already been culled before
if isfield(delaunay, 'culled') && ~isempty(delaunay.culled)
    warning(['Triangulation has already been culled: ' delaunay.culled{1} ' ' num2str(delaunay.culled{2})]);
end

% this is added by GDF on 05/14 for cases where the function "dist" is not
% available
if exist('dist')     
    distances = dist(delaunay.coordinates');
    
    switch opts.neighbor_criterion
        case 'nearest nbr median'
            sorted_distances = sort(distances);
            nearest_nbr_distances = sorted_distances(2,:); % Skip the first row, which has all zeros for distance to self
            max_dist = median(nearest_nbr_distances) * distance_factor;

        % ToDo: add 'all nbrs median', 'nearest nbr robust std'?
        otherwise
            error('maximum neighbor length criterion not recognized!  (%s)', neighbor_criterion)
    end
    
else
    % get nearest neighbor distances
    distances = ipdm(delaunay.coordinates);

    nearest_nbr_distances = zeros(size(distances,1),1);
    for rw = 1:size(distances,1)
        tmp_rw = sort(distances(rw,:));
        nearest_nbr_distances(rw) = tmp_rw(2);
    end
    
    max_dist = median(nearest_nbr_distances) * distance_factor;
end


triangles = delaunay.triangles;
allowed_triangles = [];
for i = 1:size(triangles, 1)
    tri = triangles(i,:);
    
    edge_lens(1) = distances(tri(1), tri(2));
    edge_lens(2) = distances(tri(2), tri(3));
    edge_lens(3) = distances(tri(3), tri(1));

    if ~any(edge_lens > max_dist)
        allowed_triangles(end+1, :) = tri; %#ok<AGROW>
    end
end

datarun.stas.delaunays{cell_type_num}.triangles = allowed_triangles;
datarun.stas.delaunays{cell_type_num}.culled = {opts.neighbor_criterion, distance_factor};

if opts.plot
    plot_delaunay_tri(datarun, cell_spec);
end