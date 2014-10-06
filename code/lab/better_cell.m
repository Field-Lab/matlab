function [better_cell_id, worse_cell_id] = better_cell(datarun,varargin)
% better_cell     identify which cell is better based on various criteria
%                       in case of a tie, cell id a is used
%
% usage:  [better_cell_id, worse_cell_id] = better_cell(datarun,cell_id_a,cell_id_b,varargin)
%
% arguments:      datarun - datarun struct
%    cell_id_a, cell_id_b - cell ids of the cells to be compared
%                varargin - struct or list of optional parameters (see below)
%
% outputs:     better_cell_id - cell which was better
%               worse_cell_id - cell which was worse
%
%
% optional params, their default values, and what they specify:
%
% criteria     	{'spikes','contam'}
%                               cell array listing which parameters to consider
%                                 	'spikes' - spike rate (higher is better)
%                                  	'contam' - contamination (lower is better)
%                                  	'sta snr' - SNR of the STA (higher is better)
%
%
%  parameters relating how to compare contamination
%
% window        2               temporal window in which a spike is considered contamination (msec)
%
%
%
%
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% check input
p.addRequired('cell_id_a',@(x)length(x)==1)
p.addRequired('cell_id_b',@(x)length(x)==1)

% specify list of optional parameters
p.addParamValue('criteria', {'spikes','contam'},@iscell);
p.addParamValue('window', 2);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% get required
cell_id_a = params.cell_id_a;
cell_id_b = params.cell_id_b;


% BODY OF THE FUNCTION

% go through each criterion
for cc = 1:length(params.criteria)

    % based on the criterion, decide which cell is better
    % if there's a tie, none will be assigned and the code will continue to the next criterion

    switch params.criteria{cc}
        case 'spikes'
            % get spike counts
            count_a = length(datarun.spikes{get_cell_indices(datarun,cell_id_a)});
            count_b = length(datarun.spikes{get_cell_indices(datarun,cell_id_b)});
            % compare
            switch sign(count_a - count_b)
                case 1; better_cell_id = cell_id_a; worse_cell_id = cell_id_b; return
                case -1; better_cell_id = cell_id_b; worse_cell_id = cell_id_a; return
            end


        case 'contam'
            error('contamination not implemented yet')

        case 'sta snr'
            error('STA SNR not implemented yet')

    end

    % if the better cell was found, stop looping through criteria
    if exist('better_cell_id','var')
        break
    end
    
end

% if better_cell_id was not set yet, it's a tie, so just pick cell id a
if ~exist('better_cell_id','var')
    better_cell_id = cell_id_a;
    worse_cell_id = cell_id_b;
end



