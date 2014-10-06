function [connectivity, returned_datarun, extras] = extract_connectivity(datarun, cell_type_number, varargin)

% extract_connectivity    This function returns a connectivity matrices for cells given by cell_spec
%
% usage:  connectivity = extract_connectivity(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct 
%    cell_type_number - which cell type to extract
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    connectivity - A connectivity matrix
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false               show output
% remove_cones       []                 String the specifies cone types for removal from connectivity matrix
% min_radius         0                  number of sigmas from RF fits to include cones
% max_radius         2                  number of sigmas from RF fits to include cones
% required_sign      []                 specify the required sign of the cone weights           
%                                           'positive' -- i.e. center
%                                           'negative' -- i.e. surround
% min_convergence    0                  A RGC much have at least this number to be included in connectivity
% max_convergence    200                A RGC can have no more than this number of cones to be included in connectivity
% normalize          true               determines whether the weights are normalized
%
% 2009-05-20 GDF
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('remove_cones', [], @ischar);
p.addParamValue('min_radius', 0, @isnumeric);
p.addParamValue('max_radius', 2, @isnumeric);
p.addParamValue('required_sign', []);

p.addParamValue('min_convergence', 0, @isnumeric);
p.addParamValue('max_convergence', 200, @isnumeric);
p.addParamValue('normalize', true, @islogical);
p.addParamValue('cell_list', [], @isnumeric);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% BODY OF THE FUNCTION

% removed specified cones
if ~isempty(params.remove_cones)
    for tp = 1:length(params.remove_cones)
        temp_type(tp) = params.remove_cones(tp);
        temp_exclude_cone_indices = find(datarun.cones.types == temp_type(tp));
        keep_cone_indices = setdiff([1:length(datarun.cones.types)], temp_exclude_cone_indices);
        
        datarun.cones.centers = datarun.cones.centers(keep_cone_indices,:);
        datarun.cones.types = datarun.cones.types(keep_cone_indices);
        datarun.cones.rgb = datarun.cones.rgb(keep_cone_indices,:);
        datarun.cones.likelihoods = datarun.cones.likelihoods(keep_cone_indices);
        datarun.cones.types_em = datarun.cones.types_em(keep_cone_indices);
        datarun.cones.types_kmeans = datarun.cones.types_kmeans(keep_cone_indices);
        datarun.cones.roi = datarun.cones.roi(keep_cone_indices);
        datarun.cones.weights = datarun.cones.weights(keep_cone_indices,:);
    end
end

cone_locations = datarun.cones.centers;
cone_types = datarun.cones.types;
num_cones = length(cone_types);


%%%%% RGC extraction
if isempty(params.cell_list)
    RGC_indices = get_cell_indices(datarun, {cell_type_number});
else
    RGC_indices = get_cell_indices(datarun, params.cell_list);
end
num_RGCs = length(RGC_indices);
RGC_locations = zeros(num_RGCs,2);
for RGC = 1:num_RGCs
    temp_location = datarun.cones.rf_fits{RGC_indices(RGC)}.center;
    RGC_locations(RGC,:) = temp_location;
end

% restrict cell inclusion in connectivity_matrix 
new_connectivity_matrix = zeros(num_cones, num_RGCs);
excluded_cell_counter = 0;
keeper_cells = 0;
clear convergences;
for RGC = 1:num_RGCs
    % get number of center cones and associated weights
    temp_RGC_center = RGC_locations(RGC,:);
    temp_distances = ipdm(temp_RGC_center, cone_locations);
    
    % find cones in collecting radius
    temp_fits = datarun.cones.rf_fits{RGC_indices(RGC)};
    min_RGC_radius = temp_fits.center_radius * params.min_radius; 
    max_RGC_radius = temp_fits.center_radius * params.max_radius;
    temp_cone_indices_min = find(temp_distances >= min_RGC_radius);
    temp_cone_indices_max = find(temp_distances <= max_RGC_radius);
    
    collected_cone_indices = intersect(temp_cone_indices_min, temp_cone_indices_max);
    
    if ~isempty(params.required_sign)
        if strcmp(params.required_sign, 'positive')
            signed_cone_indices = find(datarun.cones.weights(:, RGC_indices(RGC)) > 0);
        elseif strcmp(params.required_sign, 'negative')
            signed_cone_indices = find(datarun.cones.weights(:, RGC_indices(RGC)) < 0);
        end
        collected_cone_indices = intersect(collected_cone_indices, signed_cone_indices);
    end
    
    temp_num_center_cones = length(collected_cone_indices);
    
    % get only keep RGCs within a region of interest
    clear temp_cone_weights
    if temp_num_center_cones < params.min_convergence || temp_num_center_cones > params.max_convergence
        excluded_cell_counter = excluded_cell_counter + 1;
        temp_cone_weights(1:temp_num_center_cones) = 0;
        % reporters
        %datarun.cell_ids(RGC_indices(RGC))
        %RGC_indices(RGC)
    else
        temp_cone_weights = datarun.cones.weights(collected_cone_indices, RGC_indices(RGC));
        keeper_cells = 1+keeper_cells;
        convergences(keeper_cells) = temp_num_center_cones;
        orig_index_tracker(keeper_cells) = RGC;
    end
    
    if excluded_cell_counter == num_RGCs
        warning('no cell met inclusion criterion')
    end
    
    % assign to new connectivity matrix and normalize
    if params.normalize == true
        new_connectivity_matrix(collected_cone_indices, RGC) = temp_cone_weights ./ sum(temp_cone_weights);
    else
        new_connectivity_matrix(collected_cone_indices, RGC) = temp_cone_weights; 
    end
end
% remove the zero filled RGC weight vectors from the connectivity matrix.
% These zero filled weight vectors correspond to RGCs with a cone
% convergence less than min_cone_convergence or greater than
% max_cone_convergence
temp_keeper_indices = find(max(new_connectivity_matrix, [], 1) > 0);
new_connectivity_matrix = new_connectivity_matrix(:, temp_keeper_indices);
roi_RGC_locations = RGC_locations([temp_keeper_indices],:);
rgc_fit_info = datarun.cones.rf_fits(RGC_indices(temp_keeper_indices));

[num_roi_cones, num_roi_RGCs] = size(new_connectivity_matrix);

extras.mean_convergence = mean(convergences);
extras.stdv_convergence = std(convergences);
extras.convergences = convergences;
extras.original_RGC_indices = orig_index_tracker;
extras.keeper_ids = datarun.cell_ids(RGC_indices(temp_keeper_indices));

returned_datarun.stas.rf_coms = datarun.stas.rf_coms(RGC_indices(temp_keeper_indices));
returned_datarun.cones = datarun.cones;
returned_datarun.cones.weights = new_connectivity_matrix;
returned_datarun.cones.rf_fits = rgc_fit_info;
returned_datarun.cell_ids = datarun.cell_ids(RGC_indices(temp_keeper_indices));
temp_ids = datarun.cell_types{cell_type_number}.cell_ids;
temp_ids = temp_ids(temp_keeper_indices);
temp_name = datarun.cell_types{cell_type_number}.name;
temp_struct.name = temp_name;
temp_struct.cell_ids = temp_ids;
returned_datarun.cell_types{cell_type_number} = temp_struct;

% display how long it took
if params.verbose
    fprintf('done (%0.1f seconds)\n',etime(clock,start_time));
end

connectivity = new_connectivity_matrix;