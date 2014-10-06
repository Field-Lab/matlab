function datarun = get_rfs_filtered(datarun, cell_specification, varargin)
% get_rfs_filtered     Filter STA spatial summary and save
%
% usage:  datarun = get_rfs_filtered(datarun, cell_specification, params)
%
% arguments:          datarun - datarun struct with field dataset.stas.summaries
%          cell_specification - see get_cell_indices for options
%                      params - struct of optional parameters (see below)
%
% outputs:            datarun - datarun struct with fields:
%
%   datarun.summaries_filt{}    or another name
%
%
% optional fields in params, their default values, and what they specify:
%
% verbose           false               display info
% filt_params     	[]                  type of filter, see 'rf_filtered' for options
% save_filt_params  []                  save parameters to datarun.stas.(save_filt_params)
%                                           if empty, params are not saved
% save_name         'summaries_filt'    name of field in which to save summaries
%
%
%   gauthier 2008-03
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('filt_params', []);
p.addParamValue('save_filt_params', []);
p.addParamValue('save_name', 'summaries_filt');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% get list of cells
cell_nums = get_cell_indices( datarun, cell_specification);


% display how many cells will be loaded
if params.verbose;
    fprintf('\nFiltering spatial STA summary for %d cells ',length(cell_nums));
    start_loading = clock; % note when it started
end


% go through through list of cells
for cell_num = cell_nums
    
    % show a tick of progress
    if params.verbose;fprintf('.');end
    
    % compute filtered summary frame
    [rf_filtered_, filt_params] = rf_filtered(get_rf(datarun,datarun.cell_ids(cell_num)),params.filt_params);
    
    % save it
    datarun.stas.(params.save_name){cell_num} = rf_filtered_;
    
    % save params, if desived
    if ~isempty(params.save_filt_params)
        datarun.stas.(params.save_filt_params){cell_num} = filt_params;
    end
        
end

% if not all cells have summaries, ensure datarun.stas.(params.save_name) is the right length
if length(datarun.stas.(params.save_name)) < length(datarun.cell_ids)
    datarun.stas.(params.save_name){length(datarun.cell_ids)} = [];
end

% display how long it took
if params.verbose;
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
end

