function cell_types = load_rrs_cell_types(java_params, params)
% LOAD_RRS_CELL_TYPES     Get cell types from a java params file object
%
% usage:  cell_types = load_rrs_cell_types(java_params, params)
%
% arguments:  java_params - java object paramsFile
%                  params - struct of optional parameters (see below)
%
% outputs:   cell_types - cell types in standard cell array format
%
%
% optional fields in params, their default values, and what they specify:
%
% cell_type_depth       2       how many levels deep to look for cell types in the classification or 'all'
%
% See also: LOAD_TXT_CELL_TYPES, ORDER_CELL_TYPES
%
% gauthier  2008-03
% greschner 2009-02 include key 'cell_type_depth'='all'
% phli      2013-03 fixed so that idmap is not continually regenerated;
%                   much faster for params with many IDs
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.cell_type_depth = 2;

% combine user and default parameters
params = default_params( defaults, params);



% initialize temporary variable
cell_types_temp = struct;

% cycle through cells and get cell type string
idmap = java_params.getClassIDs();
for cell_id = java_params.getIDList'
    % get cell type of this cell
    type_string = idmap.get(uint32(cell_id));
    
    % if there is no class, call it 'unclassified'
    if strcmp(type_string,'All')
        type_string = 'unclassified';
    else

        % cell type names are typically of one of these forms:
        %   All                         unclassified
        %   All/junk                    first level type
        %   All/ON/parasol              second level type
        %   All/ON/parasol/strange      third level type
        %
        %   Cell types are saved out to the level specified by the user.  For example, if
        %   cell_type_depth = 2, then in the example above the last two cells would be in the same class.
        %   A good default value is 2.
        %

        % convert string to be matlab friendly at the appropriate depth

        % find locations of all '/' characters
        slashes = strfind(type_string,'/');
        % delete anything beyond the desired depth
        if all(isequal(params.cell_type_depth,'all'))
            params.cell_type_depth=length(slashes)+1;
        end
        if length(slashes) > params.cell_type_depth
            type_string = type_string(1:slashes(params.cell_type_depth + 1) - 1);
        end
        % remove 'All/', and change to be matlab friendly
        type_string = strrep(strrep(strrep(strrep(strrep(strrep(type_string,'All/',''),'/','_'),' ','_'),'-','_'),'?',''),'.','_');
    end
    
    % store in temporary variable
    if ~isfield(cell_types_temp,type_string)
        cell_types_temp.(type_string) = cell_id;
    else
        cell_types_temp.(type_string) = [cell_types_temp.(type_string) cell_id];
    end
end



% generate final cell types variable

% get cell type names
type_names = fieldnames(cell_types_temp);
% enter each into the final variable
for nn = 1:length(type_names)
    % replace _ with space, and put into growing cell array
    cell_types{nn} =  ...
        struct('name',strrep(type_names{nn},'_',' '),...
        'cell_ids', cell_types_temp.(type_names{nn}) );
end

