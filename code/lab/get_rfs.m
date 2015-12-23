function datarun = get_rfs(datarun, cell_specification, varargin)
% get_rfs     Compute STA spatial summaries based on datarun.stas.stas{} and save them in datarun
%
% usage:  datarun = get_rfs(datarun, cell_specification, params)
%
% arguments:  datarun - datarun struct, with fields datarun.cell_ids and datarun.stas.stas
%  cell_specification - e.g. [861], 'ON parasol', {1}, 'all' (see get_cell_indices)
%              params - struct of optional parameters (see below)
%
% outputs:    datarun - datarun struct with field datarun.sta.summaries
%
%
% optional fields in params, their default values, and what they specify:
%
% rf_params             []      what type of spatial summary to compute
%                                   see rf_from_sta for options
%                                   if empty, the default in rf_from_sta is used
% verbose               false
%
%
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addParamValue('method',     'project', @ischar)
p.addParamValue('rf_params',   struct,   @isstruct);
p.addParamValue('sig_stixels', cell(0),  @iscell);
p.addParamValue('verbose',     false,    @islogical)

p.parse(varargin{:});
params = p.Results;

%%% commented out by GDF ( seems unnecessary )  %%%
% parse out params for function rf_from_sta
%rf_params.method = params.method;
%rf_method_fields = fieldnames(params.method_params);
%for field = 1:length(rf_method_fields);
%    rf_params = setfield(rf_params, rf_method_fields{field}, getfield(params, rf_method_fields{field}));
%end


% get list of cell IDs for which to compute STA summary
cell_nums = get_cell_indices(datarun, cell_specification);


% display how many cells will be loaded
if params.verbose;
    fprintf('\nComputing spatial summary for %d cells ',length(cell_nums));
    start_loading = clock; % note when it started
end


% loop through cells
for cell_num = cell_nums
    % show a tick of progress
    if params.verbose;fprintf('.');end
    
    % if sig stixels provided for each cell, place in params.rf_params
    if ~isempty(params.sig_stixels)
        params.rf_params.sig_stixels = params.sig_stixels{cell_num};
    end

    % get the STA
    sta = datarun.stas.stas{cell_num};
    
    % compute the spatial summary
    datarun.stas.rfs{cell_num} = rf_from_sta(sta, params.rf_params);
end


% display how long it took
if params.verbose;
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
end

