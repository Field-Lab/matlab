function datarun = clean_up_cell_types(datarun, cell_types, varargin)
% CLEAN_UP_CELL_TYPES     modify cell type definitions to remove bad cells
%
% usage:  datarun = clean_up_cell_types(datarun, cell_types, <params>)
%
% arguments:  datarun - datarun struct with field specifying cell types
%          cell_types - list of cell types to clean up (numbers only)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in datarun.cell_types
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           true            show output
% fig               []              where to plot mosaics before and after
% remove            []              list of cell ids to remove
% snr               []              minimum SNR.  if empty, SNRs are plotted
%                                       in a new figure and the user is prompted to specify the threshold 
%
%
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('fig', 0);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



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

    
    
    % REMOVE SPECIFIED CELLS


    
    
    

    % GATHER QUALITY INFORMATION ABOUT EACH CELL


    
    
    
    
    % REMOVE CELLS WHICH ARE LOW QUALITY


    
    
    
    % REMOVE DUPLICATES

    % identify duplicates
    [unique_cells, duplicates, unsure] = ...
        identify_duplicate_rfs(datarun, cell_types{ct}.cell_ids);
    
    % toss duplicates
    good_cells = set_diff(good_cells,duplicates);
    
    % optional: toss unsure cells
    if params.toss_unsure
        good_cells = set_diff(good_cells,unsure);
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








%function datarun = clean_up_cell_types(datarun, cell_types, params)
% CLEAN_UP_CELL_TYPES     modify cell type definitions to remove bad cells
%
%   Cells are removed based on multiple criteria, each of which is optional.
%   First, each cell is identified as "good" or "bad" based on a user criterion.
%   "Bad" cells are removed from the cell type.  
%   Second, duplicates are removed, based purely on RF center location.
%   
%
%
%
% usage:  datarun = clean_up_cell_types(datarun, cell_types, params)
%
% arguments:  datarun - datarun struct with field specifying cell types
%          cell_types - list of cell types to clean up (numbers only)
%              params - struct of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional fields in params, their default values, and what they specify:
%
% criterion         'snr'       criterion for defining cell as "good" or "bad"
%                                   'snr'          - remove cells with low snr
%                                   'remove_cells' - remove specific cells specified by number
%
%
% the following fields depend on the criterion
%
%   for criterion 'snr'
% min_snr           []          minimum SNR to be considered "good"
%                                   if empty, SNRs are plotted in a new figure and the user is prompted
%                                   to specify the threshold
%
%   for criterion 'remove_cells'
% cells             []
%
% 
% gauthier 2008-10
%
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.criterion = 'snr';
defaults.min_snr = [];
defaults.cells = [];

% combine user and default parameters
params = default_params( defaults, params);





% BODY OF THE FUNCTION

% go through each of the listed cell types
for ct_num = 1:length(cell_types)
    
    % get list of cell IDs
    curr_cell_ids = datarun.cell_types{cell_types(ct_num)}.cell_ids;
    
    % get cell numbers
    cell_indices = get_cell_indices(datarun,curr_cell_ids);
    
    
    % GET MEASURE OF QUALITY FOR EACH CELL
    
    switch params.criterion
        case {'snr'}

            % make list of cell qualities
            cell_quality = zeros(length(cell_indices));

            % determine which cells are good/bad
            for cc = 1:length(cell_indices)

                switch params.criterion

                    case 'snr'

                        % get RF
                        rf = get_rf(datarun,datarun.cell_ids(cell_indices(cc)));

                        % if no RF, skip
                        if isempty(rf);
                            cell_quality = 0;
                            continue;
                        end

                        % measure SNR
                        snr = max(max(max(abs(rf)))) / robust_std(reshape(rf,1,[]));

                        %
                        cell_quality(cc) = snr;

                end
            end
            
            
            
            % CHOOSE CELLS

            % clear list of good cells
            good_cells = [];

            % if threshold was specified and SNR is sufficiently high, keep it
            if snr >= params.min_snr
                good_cells = [good_cells cc];
            end
            
            
        case 'remove_cells'

            good_cells = setdiff(curr_cell_ids,params.cells);
    end

    


    
    
    
    
    
    % UPDATE CELL LISTS
    
    % if any cells were removed
    if length(good_cells) < length(curr_cell_ids)

        % store a copy of previous list

        % get number of cell types at present
        curr_num_types = length(datarun.cell_types);

        % store a copy of the cell type
        datarun.cell_types{curr_num_types + 1} = datarun.cell_types{cell_types(ct_num)};

        % with a new name
        datarun.cell_types{curr_num_types + 1}.name = [datarun.cell_types{cell_types(ct_num)}.name ' all'];


        % update the new list
        % set list of cell types to be only cells to keep
        datarun.cell_types{cell_types(ct_num)}.cell_ids = curr_cell_ids(good_cells);
    
    end
    
    % if no cells were removed, do nothing
    
    
    % show how many were kept
    fprintf('%s: kept %d of %d\n',datarun.cell_types{cell_types(ct_num)}.name,length(good_cells),length(curr_cell_ids))
    
end

    