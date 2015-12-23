function [unique_cells, duplicates, unsure] = identify_duplicate_rfs(datarun, cell_spec, varargin)
% identify_duplicate_rfs     Identify RFs which are duplicates, based on several optional criteria
%
%  NOTE: currently the only supported criterion is RF center distance.
%
%
% usage:  result = identify_duplicate_rfs(datarun,cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     unique_cells - cell ids of unique cells
%                duplicates - cell ids of cells which are inferior copies of unique cells
%                    unsure - cell ids of cells which did not have enough info to decide
%                               e.g. no center point, or no RF
%
%
% optional params, their default values, and what they specify:
%
% verbose           false       show output
% foa_scatter       []          figure or axes to plot scatter of dot product vs distance 
%                                   if 0, make new figure. if empty, don't plot
% radius            0.6         maximum distance between center points for cells to be considered distance
%                                   measured in stixel coordinates
%
%
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('foa_scatter', []);
p.addParamValue('radius', 0.6);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
pa_scatter = set_up_fig_or_axes(params.foa_scatter);

% show output
if params.verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end



% GET CENTER POINTS AND FIND NEAR NEIGHBRORS

% get centers of cells which have a center
[centers,cell_ids] = rf_centers(datarun,cell_spec);

% identify cell pairs with neighboring distances less than the radius
dists = ipdm(centers,'subset','max','limit',params.radius,'result','array');



% PLOT (IF DESIRED)
if ~isempty(pa_scatter)
    
    % get dot products and distances of nearby cells

    dlist=[];
    clist=[];

    % cycle through each pair of nearby cells
    for cc=1:length(cell_ids)
        for dd = find(dists(cc,:)>0 & dists(cc,:)<Inf)
            rfa = reshape(get_rf(datarun,cell_ids(cc)),[],1);
            rfb = reshape(get_rf(datarun,cell_ids(dd)),[],1);
            dlist = [dlist dists(cc,dd)];
            clist = [clist sum((rfa/norm(rfa)).*(rfb/norm(rfb)))];
        end
    end

    % plot
    axes(pa_scatter)
    plot(dlist,clist,'.')
end



% DECIDE WHAT CRITERION WILL BE USED FOR DEFINING DUPLICATES

% make function to detect duplicates
switch 1
    case 1  % distance cutoff only
        dups = @(dists,unique_cells) find(dists(unique_cells,unique_cells)<Inf & dists(unique_cells,unique_cells)>0);
end



% ACCUMULATE LIST OF DUPLICATES

% initialize list of unique cells
unique_cells = true(size(dists,1),1);

% find pairs of cells closer than the cutoff distance
while ~isempty(dups(dists,unique_cells))    

    % get one pair of duplicates
    [cc,dd] = dups(dists,unique_cells);
    cc = cc(1); dd = dd(1);
    
    % find the better cell
    unique_cell_ids = cell_ids(unique_cells);
    [better, worse] = better_cell(datarun,unique_cell_ids(cc),unique_cell_ids(dd));
    
    % exclude the worse from the list of unique cells
    unique_cells(worse==cell_ids) = false;
end

% set output
duplicates = cell_ids(~unique_cells);
unique_cells = cell_ids(unique_cells);
unsure = setdiff(datarun.cell_ids(get_cell_indices(datarun,cell_spec)),[unique_cells duplicates]);

