function cell_types = remove_duplicate_cells(datarun, cell_types, varargin)
% remove_duplicate_cells     Identify duplicates from a cell type struct
%
% usage:  cell_types = remove_duplicates(datarun, cell_types, <params>)
%
% arguments:  datarun - datarun struct
%          cell_types - cell types struct
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field X
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false           show output
% types             1:6             which cell types to look through
% toss_unsure       true            toss cells for which duplicate status can not be determined?
% fig               []              where to plot mosaics
%
%
%
% parameters passed on to other functions. if not specified by the user, these parameters are not passed.
%
%
%       passed to 'identify_duplicate_rfs'
%
%           radius       radius     	cutoff radius for considering cells to be duplicates
%
%
%
%
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('types', 1:6);
p.addParamValue('toss_unsure', true);
p.addParamValue('fig', []);

% parameters to be passed on
%    identify_duplicate_rfs
p.addParamValue('radius', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
dup_params = make_struct_to_pass(p.Results,{'radius','radius'});





% BODY OF THE FUNCTION


% initialize output with input
cell_types_out = cell_types;

% note when it started
if params.verbose
    start_time = clock; 
    fprintf('\nRemoving duplicates... ');
end

% set up figure
if ~isempty(params.fig)
    set_up_fig_or_axes(params.fig);
    plot_axes = subplot_axes(params.fig,[0 0 1 1],0.05,0.1,3,2);
end
    
% go through each cell type
for cc=1:length(params.types)

    % note the cell type
    ct = params.types(cc);
    
    % plot mosaic of all RFs
    if ~isempty(params.fig)
        plot_rf_summaries(datarun,cell_types{ct}.cell_ids,'foa',plot_axes{ct})
        title(cell_types{ct}.name)
        drawnow
    end
    
    % skip if no cells are present
    if isempty(cell_types{ct}.cell_ids);continue;end
    
    % show output
    if params.verbose
        fprintf('\n\t%s: ',cell_types{ct}.name);
    end

    % identify duplicates
    [unique_cells, duplicates, unsure] = ...
        identify_duplicate_rfs(datarun, cell_types{ct}.cell_ids);
    
    % keep only the unique cells
    cell_types_out{ct}.cell_ids = unique_cells;
    
    % optional: keep unsure cells
    if ~params.toss_unsure
        cell_types_out{ct}.cell_ids = [cell_types_out{ct}.cell_ids unsure];
    end
    
    % plot in red RFs which were removed
    if ~isempty(params.fig)
        removed_cells = setdiff(cell_types{ct}.cell_ids,cell_types_out{ct}.cell_ids);
        plot_rf_summaries(datarun,removed_cells,'foa',plot_axes{ct},'fit_color','r','clear',0)
        title(cell_types{ct}.name)
        drawnow
    end
    
    % show output
    if params.verbose
        fprintf('kept %d/%d',length(cell_types_out{ct}.cell_ids),length(cell_types{ct}.cell_ids));
    end
    
end

% display how long it took
if params.verbose
    fprintf('\n    done (%0.1f seconds)\n',etime(clock,start_time));
end

% return output
cell_types = cell_types_out;

