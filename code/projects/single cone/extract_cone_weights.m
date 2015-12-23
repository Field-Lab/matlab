function [cone_weights,cell_ids,Wc] = extract_cone_weights(datarun, cell_spec, cone_spatial_profiles, cone_types, cone_rgb_observed, params)
% extract_cone_weights     Identify strength of each cone in each RF using regression of the cone spatial profiles
%
% usage:  cone_weights = extract_cone_weights(datarun, cell_spec, cone_spatial_profiles,...
%                            cone_types, cone_rgb_observed, params)
%
% arguments:      datarun - datarun struct
%               cell_spec - which cells to use
%   cone_spatial_profiles - MxN matrix: pixel maps of the ideal shape of each cone
%              cone_types - N-length char vector of cone colors ('L','M','S', or 'U' for unknown)
%       cone_rgb_observed - struct with fields 'L','M','S', giving mean RGB values for each cone type
%                  params - struct of optional parameters (see below)
%
% outputs:     cone_weights - NxM matrix of cone weights (cones = rows, cells = columns)
%                  cell_ids - list of cell ids which are in the weight
%                  matrix
%                        Wc - cone weights of each RGC
%
%
% optional fields in params, their default values, and what they specify:
%
% verbose           true             	show output
%
%
% gauthier  2008-10
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.verbose = true;

% combine user and default parameters
params = default_params( defaults, params);





% BODY OF THE FUNCTION


% get list of cell numbers
cell_nums = get_cell_indices(datarun,cell_spec);


if params.verbose
    fprintf('\nComputing cone weights in %d RFs',length(cell_nums));
    start_time = clock; % note when it started
end



% set cone weights for regression

% set RGB weights for unknown cones as the mean of L and M cones
if ~isfield(cone_rgb_observed, 'U') || isempty(cone_rgb_observed.U)
    cone_rgb_observed.U = mean([cone_rgb_observed.L; cone_rgb_observed.M]);
end

% normalize cone weights to have the same L1 norm (sum of absolute values)
cone_names = fieldnames(cone_rgb_observed);
for cc = 1:length(cone_names)
    rgb = cone_rgb_observed.(cone_names{cc});
    cone_rgb_observed.(cone_names{cc}) = rgb/sum(abs(rgb));    
end



% make a matrix of cone spatial profiles, including the appropriate RGB sensitivity for each cone

% get number of cones
num_cones = size(cone_spatial_profiles,2);

% initialize
Wc = sparse(3 * datarun.stimulus.field_height * datarun.stimulus.field_width, num_cones);
    
% for each cone
for nn = 1:num_cones
    % get cone profile
    cn = cone_spatial_profiles(:,nn);

    % get cone RGB weights of this particular cone
    rgb = cone_rgb_observed.(cone_types(nn));
    
    % apply these to the spatial profile
    Wc(:,nn) = [ rgb(1)*cn; rgb(2)*cn; rgb(3)*cn ];
end



% regress to get the cone weights for each RGC

cone_weights = zeros(num_cones,length(cell_nums));

% go through list of cells
for cc = 1:length(cell_nums)
    
    fprintf('.')

    % get summary frame
    rf = get_rf(datarun,datarun.cell_ids(cell_nums(cc)));

    % check if rf is BW
    if strcmp(datarun.stimulus.independent, 'nil')
        rf = repmat(rf, [1 1 3]);
    end        
    
    if isempty(rf)
        continue
    end

    % reshape for the regression
    rf = reshape(rf,[],1);

    % put in units of SNR
    rf = rf / robust_std(rf);

    % regress to get the cone weights
    cone_weights(:,cc) = Wc\rf;

end



% display how long it took
if params.verbose
    fprintf('\n   done (%0.1f seconds)\n',etime(clock,start_time));
end


% return list of cell ids
cell_ids = datarun.cell_ids(cell_nums);



