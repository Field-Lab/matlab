function d = nearest_neighbor_distances(varargin)

% NEAREST_NEIGHBOR_DISTANCES(cell_type, data_directory) returns the
% distribution of nearest neighbor distances for cells of specified type.
%
% Inputs:   
%   cell_type = 1 for OnP, 2 for OffP, 3 for OnM, 4 for OffM
%       POTENTIAL BUG: assumes that this is the order of cell types in the
%       data set.
%   data_directory = dataset to load (eg '2012-08-09-3/data003')
%       If not specified, will assume 'datarun' exists already in workspace.
%
% Output: Vector of distance to nearest neighbor for all cells of specified
% type.
%
% Malcolm Campbell 2014
% malcolmc@stanford.edu

% load data
if length(varargin) == 2
    datarun = load_data(varargin{2});
    datarun = load_sta(datarun);
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    datarun = set_polarities(datarun);
else
    datarun = evalin('base','datarun');
end

% find cells of specified type and extract rf centers
cell_type = varargin{1};
cell_ids = find(ismember(datarun.cell_ids,datarun.cell_types{1,cell_type}.cell_ids));
num_cells = length(cell_ids);
rf_centers = zeros(num_cells,2);
for i=1:num_cells
    rf_centers(i,:) = datarun.stas.fits{1,cell_ids(i)}.mean;
end

% calculate paired distances between cells
paired_distances = zeros(num_cells,num_cells);
for i=1:num_cells
    for j=1:num_cells
        paired_distances(i,j) = sqrt(sum((rf_centers(i,:)-rf_centers(j,:)).^2));
    end
end

% find nearest neighbors, masking zeros with Inf
paired_distances(paired_distances==0) = Inf;
d = min(paired_distances,[],2);

end