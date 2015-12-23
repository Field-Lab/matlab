function  [initial_cone_centers,cone_sources] = find_cones_in_each_rf(datarun, cell_spec, all_sig_stixels, varargin)
% find_cones_in_each_rf     Find cones in each of many RFs, and pool across cells
%                               to find a single collection
%
% usage:  [cones_labeled,initial_cone_centers] = ...
%           find_cones_in_each_rf(datarun, cell_spec, all_sig_stixels, <params>)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_indices for options)
%     all_sig_stixels - see output of compute_spatial_sensitivity
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    cones_labeled - YxX matrix, equal to size of a single STA.  Mostly 0s, and at each cone location
%                               indicates the cone id
%      initial_cone_centers -  Nx2 matrix: x,y locations describing approximate cone centers
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
% foo               'bar'               how to activate the rings
%                                           'bar' - activate on site
%                                           'bore' - activate remotely
%
%
% parameters passed on to 'find_cones_in_rf'. if not specified by the user, these parameters are not passed.
%
%       strength                how to summarize RGB values in the RF
%       filter                  how to filter RF before looking for cones
%       selection_params      	struct of params to pass to function find_cones_in_rf
%       fig_single_cell         (passed as fig_rfs)   figure number to plot each RF.
%
%
% See also: FIND_CONES_IN_RF, COMPUTE_SPATIAL_SENSITIVITY
%
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('foo','bar', @(x)any(strcmpi(x,{'bar','bore'})));

% parameters to be passed on to find_cones_in_rf
p.addParamValue('strength', 'default value');
p.addParamValue('filter', 'default value');
p.addParamValue('selection_params', 'default value');
p.addParamValue('fig_single_cell', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% make params struct to pass to find_cones_in_rf
find_params = make_struct_to_pass(p.Results,{'selection','selection_params',...
    'filter','filter','fig_rfs','fig_single_cell','strength','strength'});





% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% get cell indices and ids
cell_indices = get_cell_indices(datarun,cell_spec);

% error check
if length(cell_indices) ~= size(all_sig_stixels,2)
    error('variable ''all_sig_stixels'' contains %d cells, but should have %d cells',...
        size(all_sig_stixels,2),length(cell_indices))
end

% show output
if params.verbose
    fprintf('\nFinding cones in %d cells',length(cell_indices));
    start_time = clock; % note when it started
end

% initialize variable storing cone locations
cone_locs = [];
%cone_locs = zeros(datarun.stimulus.field_height,datarun.stimulus.field_width);
cone_sources = [];

% loop through cells
for cc = 1:length(cell_indices)
    
    if params.verbose;fprintf('.');end
    
    % get cell index and id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % get rf
    rf = get_rf(datarun,cell_id);
    
    % if empty, skip
    if isempty(rf);continue;end
    
    % get sig stixels and rf strength
    [sig_stixels,rf_strength] = find_cones_in_rf(rf,find_params);
    
    % expand sig_stixels to include a slightly larger area
    %sig_stixels = imfilter(sig_stixels,ones(3,3));
    
    % find local maxima among the significant stixels
    temp_locs = find_local_maxima(rf_strength.*sig_stixels,'return','coordinates');
    
    % get COM of each one
    for tt=1:size(temp_locs,1)
        % get only stixels around the max
        temp = zeros(size(sig_stixels));temp(temp_locs(tt,2),temp_locs(tt,1))=1;
        local_stixels = rf_strength.*imfilter(temp,ones(3,3));
        % find COM
        [ii,jj] = ait_centroid(local_stixels);
        % add to growing list of cone locations
        cone_locs = [cone_locs; ii jj];
        % note source
        cone_sources = [cone_sources cell_id];
    end
        
    % find local maxima among the significant stixels
    %cone_locs = [cone_locs; find_local_maxima(rf_strength,'return','indices')];
    %cone_locs = cone_locs + find_local_maxima(rf_strength);

    if ~isempty(cone_locs)
        figure(1);clf;plot(cone_locs(:,1),cone_locs(:,2),'.');drawnow
    end
    %figure(1);clf;imagesc(cone_locs);drawnow;axis image
end


% display how long it took
if params.verbose
    fprintf('\n   done (%0.1f seconds)\n',etime(clock,start_time));
end

initial_cone_centers = cone_locs;

