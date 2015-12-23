function datarun = set_polarities(datarun, varargin)
% SET_POLARITIES     Set STA polarity for all cells within a type
%
% usage:  datarun = set_polarities(datarun, params)
%
% arguments:  datarun - datarun struct with cell types defined
%              params - struct of optional parameters (see below)
%
% outputs:    datarun - datarun struct with entries in datarun.stas.polarities{}
%
%
% optional fields in params, their default values, and what they specify:
%
% cell_specs     	{{1,3,5},{2,4}} 	cell array of cell specifications (see get_cell_indices for options)
% polarities        [1 -1]              vector of respective polarities to assign to these groups of cells
% guess             false               guess the polarity of each type based on its name
%                                           if true, the arguments 'cell_specs' and 'polarities' are ignored
%
%
% NOTE: to change the polarity of a single cell, use
%
%       datarun.stas.polarities{cell_index} = <polarity>
%
%   Created JLG
%   Edited GDF = made it conform to varargin standards (inputParser)



% SET UP OPTIONAL ARGUMENTS
p = inputParser;

% specify list of optional parameters
p.addRequired('datarun', @isstruct);
p.addParamValue('cell_specs', {{1,3,5},{2,4}});
p.addParamValue('polarities', [1 -1] ,@isnumeric);
p.addParamValue('guess', false);
p.parse(datarun, varargin{:})

params = p.Results;

cell_specs = p.Results.cell_specs;
polarities = p.Results.polarities;



% ensure datarun.cell_types exists
if ~isfield(datarun,'cell_types')
    error('cell types not loaded, use datarun = load_params(datarun)')
end

% ensure datarun.stas.polarities has the right length
cell_count = length(datarun.cell_ids);
if ~isfield(datarun,'stas') || ~isfield(datarun.stas,'polarities')
    datarun.stas.polarities = cell(cell_count,1);
else
    % ... and is the correct length
    if length(datarun.stas.polarities) < cell_count
        datarun.stas.polarities{cell_count} = [];
    end
end


if params.guess

    % go through each cell type and fill in the polarity based on its name
    for tt = 1:length(datarun.cell_types)
        % check the name
        if regexpi(datarun.cell_types{tt}.name,'on')
            [datarun.stas.polarities{get_cell_indices(datarun,{tt})}] = deal(1);
            
        elseif regexpi(datarun.cell_types{tt}.name,'off')
            [datarun.stas.polarities{get_cell_indices(datarun,{tt})}] = deal(-1);
        end

    end

else


    % ensure lengths are correct
    if length(cell_specs) ~= length(polarities)
        error('set_polarities: length of cell_specs must match length of polarities')
    end

    % for each group of cells
    for ss = 1:length(cell_specs)

        % get list of cell IDs
        cell_nums = get_cell_indices(datarun, cell_specs{ss});

        % set each cell's polarity
        for cell_num = cell_nums
            datarun.stas.polarities{cell_num} = polarities(ss);
        end
    end
end
