function datarun = get_rfs_normalized(datarun, cell_spec, varargin)
% get_rfs_normalized     normalize the amplitude of RFs
%
% datarun = get_rfs_normalized(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field X
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false          	show output
% fig_or_axes       []             	figure or axes to plot in. if 0, make new figure. if empty, don't plot
%
% method        	'variance'      how to get the RF
%                                       'variance' - set each RF to have unit variance within a ROI
%                       
%                         
%  if method == 'variance'
%       center_type 'rf com'        center point of the ROI
%       radius    	3               radius of the ROI, in units
%
%
%
% parameters passed on to other functions. if not specified by the user, these parameters are not passed.
%
%
%       passed to 'compute_temporal_matrix'
%
%           time_bins       bins     	number of bins
%           frames          frames     	which sta frames to use 
%
%
%
%
%
%
% 2008-12 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('save_field', 'rfs');

% specify list of optional parameters
p.addParamValue('method','variance', @(x)any(strcmpi(x,{'variance','peak'})));

% add other parameters based on 'method'
p.KeepUnmatched = true;
p.parse(varargin{:});
switch p.Results.method
    case 'variance'
        p.addParamValue('roi_center', 'com');
        p.addParamValue('roi_radius', 3);
end
p.KeepUnmatched = false;


% parameters to be passed on
%    rf_normalized
p.addParamValue('roi', 'default value');
p.addParamValue('peak', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
norm_params = make_struct_to_pass(p.Results,{'method','method','peak','peak'});





% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);



% ensure field to store normalized RFs exists...
if ~isfield(datarun.stas,params.save_field)
    datarun.stas.(params.save_field) = cell(length(datarun.cell_ids),1);
else
    % ... and is the correct length
    if length(datarun.stas.(params.save_field)) < length(datarun.cell_ids)
        datarun.stas.(params.save_field){length(datarun.cell_ids)} = [];
    end
end



% show output
if params.verbose
    fprintf('\nNoralizing RF for %d cells...',length(cell_indices));
    start_time = clock; % note when it started
end


% loop through cells
for cc = 1:length(cell_indices)
    
    % get cell index and id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % get the rf
    rf = get_rf(datarun,cell_id);
    
    % get normalized rf
    switch params.method
        case 'variance' % normalize by the variance within a ROI
            
            % identify the ROI
            ctr = rf_center(datarun,cell_id,params.roi_center);
            roi = distance_from_point([size(rf,1) size(rf,2)], ctr) <= params.roi_radius;
            
            % normalize by variance within this ROI
            rf_norm = rf_normalized(rf,norm_params,'roi',roi);
            
        case 'peak' % normalize by the peak
            rf_norm = rf_normalized(rf,norm_params);
    end
    
    % save it in datarun
    datarun.stas.(params.save_field){cell_index} = rf_norm;
    
end


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

