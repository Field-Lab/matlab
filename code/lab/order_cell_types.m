function cell_types = order_cell_types(cell_types, params)
% ORDER_CELL_TYPES     Try to put cell types in the canonical order:
%
% 	1	ON parasol
% 	2	OFF parasol
% 	3	ON midget
% 	4	OFF midget
%   5   SBC
%   6   ac
%
% usage:  cell_types = order_cell_types(cell_types, params)
%
% arguments:  cell_types - cell types in standard cell array form
%                 params - struct of optional parameters (see below)
%
% outputs:    cell_types - reordered version of input
%
%
% optional fields in params, their default values, and what they specify:
%
% req_order    {'on par','off par','on mid','off mid',{'sbc','blue'}}
%                           partial strings of desired cell types in the desired order.
%                               If an entry is a cell array, then the code will search
%                               through the strings in order until a cell type with
%                               the desired string is found.
%                           If a cell type is not found, the element will be an empty type
%
%
% avoid         'contam'	string to avoid when searching for a known string
%                               for example, if there are two types called
%                               'ON parasol' and 'ON parasol contaminated', the
%                               first will be used, because 'contam' will be avoided
%
% Note: to avoid errors, the code operates in two parts:
%   First, cell type names are read in, and the new desired order is generated
%   Second, this desired order is checked for mistakes (e.g. one cell type left out
%       or included twice).  If there are no mistakes, then the new cell types
%       are generated and returned.
%
%
% See also: LOAD_TXT_CELL_TYPES, LOAD_RRS_CELL_TYPES
%
% gauthier 2008-03
% 2010     greschner add AC
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.req_order = {{'on par','on_par','on_y','on y'},...
    {'off par','off_par','off y','off_y'},...
    {'on mid','on_mid','on x','on_x'},...
    {'off mid','off_mid','off x','off_x'},...
    {'sbc','blue'},...
    {'amacrine','ac','off ac','ac off','off-ac','off_ac'},...
    {'lbc','slow blue','large blu','big blu'}};

defaults.avoid = 'contam';

% combine user and default parameters
params = default_params( defaults, params);


% if cell_types is empty, return empty
if length(cell_types) == 0
    return
end


% avoid name confusion
in_types = cell_types;
clear cell_types


% GENERATE THE DESIRED ORDER

% get cell type names
for tt = 1:length(in_types)
    type_names{tt} = in_types{tt}.name;
end

% initialize variable that will store the new ordering of cell types
new_order = [];

% cycle through the desired cell types
for dd = 1:length(params.req_order)
    % determine whether there is a single string to search for, or multiple
    switch class(params.req_order{dd})
        
        case 'char' % a single type
            % look for the specified string in the list of cell types
            new_order(dd) = name_location(type_names,params.req_order{dd},params.avoid);

        case 'cell' % multiple types
            % get list of names to search for
            names = params.req_order{dd};
            
            % loop through them, searching each name
            for nn = 1:length(names)
                % determine whether the string is there or not
                loc = name_location(type_names,names{nn},params.avoid);
                % if the string is there
                if loc
                    % set it in the order
                    new_order(dd) = loc;
                    % and stop looping
                    break
                end
            end
            
            % if the cell type was NOT found
            if length(new_order) < dd
                % note it
                new_order(dd) = 0;
            end
    end
end

% add in to the order any cell types which were not among the desired types
left_over = setdiff( 1:length(type_names), new_order );
if length(left_over)>1
    new_order = [new_order left_over([2:end 1])];%put unclassif at the end
else
    new_order = [new_order left_over];
end


% RETURN A CELL ARRAY WITH THAT ORDER

% ensure new_order is a permutation of the numbers 1 through n (no duplicates,
% no numbers skipped), perhaps with some extra 0's
if ~all(sort(new_order(new_order > 0)) == 1:length(new_order(new_order > 0)))
    error('order_cell_types: Cell types could not be ordered (bug in the code?).  Please order them manualy.')
else
    % if everything checks out, then go ahead
    for nn = 1:length(new_order)

        % if it's not cell type 0
        if new_order(nn) > 0
            % set it to be the desired cell type
            cell_types{nn} = in_types{new_order(nn)};
        else
            % otherwise set it to be an empty type
            cell_types{nn} = struct('name','','cell_ids',[]);
        end

    end
end


function location = name_location(names,find_me,avoid_me)
% Return the ordinality of the name which contains the string find_me
% and does NOT contain the string avoid_me.  If it doesn't exist, return 0.


% go through the names
for nn = 1:length(names)
    % if the find_me is in there AND avoid_me is not
    if ~isempty(regexpi(names{nn},find_me)) && isempty(regexpi(names{nn},avoid_me))
        location = nn;
        return
    end
end

% if it wasn't found, return 0
location = 0;
