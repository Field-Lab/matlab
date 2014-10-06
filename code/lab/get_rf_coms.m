function datarun = get_rf_coms(datarun, cell_specification, varargin)
% get_rf_coms     Compute and save center of mass of the STA spatial summary
%
% usage:  datarun = get_rf_coms(datarun, cell_specification, params)
%
% arguments:          datarun - datarun struct with field dataset.stas.summary
%          cell_specification - see get_cell_indices for options
%                      params - struct of optional parameters (see below)
%
% outputs:            datarun - datarun struct with these fields:
%
%   datarun.com             center of mass
%   datarun.com_params      parameters of computation
%
%
% optional fields in params, their default values, and what they specify:
%
% verbose           false   show what's going on
% com_params        []      parameters of the computation (see rf_com)
%                               if empty, reverts to defaults in rf_com
%                               
%   gauthier 2008-03
%




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('com_params', struct);
p.addParamValue('sig_stixels', cell(0), @iscell);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% BODY OF FUNCTION

% get list of cells
cell_nums = get_cell_indices( datarun, cell_specification);


% display how many cells will be loaded
if params.verbose;
    fprintf('\nComputing COM for %d cells ',length(cell_nums));
    start_loading = clock; % note when it started
end


% go through through list of cells
for cell_num = cell_nums

    % show a tick of progress
    if params.verbose;fprintf('.');end
    
    % if sig stixels provided for each cell, place in params.rf_params
    if ~isempty(params.sig_stixels)
        params.com_params.sig_stixels = params.sig_stixels{cell_num};
    end

    % compute COM
    com = rf_com(get_rf(datarun,datarun.cell_ids(cell_num)), params.com_params);
    
    % save it
    datarun.stas.rf_coms{cell_num} = com;
        
end


% display how long it took
if params.verbose;
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
end

