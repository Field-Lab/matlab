function cell_types = load_txt_cell_types(file_path, varargin)
% load_txt_cell_types     Get cell types from a classification text file
%
% usage:  cell_types = load_txt_cell_types(file_path, varargin)
%
% arguments:    file_path - path to text file
%               varargin - struct or list of optional parameters (see below)
%
% outputs:   cell_types - cell types in standard cell array format
%
%
% optional arguments, their default values, and what they specify:
%
% cell_type_depth       2       how many levels deep to look for cell types in the classification
% order_cell_types      true    order cell types after loading
%
%
% examples:
%
%   datarun.cell_types = order_cell_types(load_txt_cell_types(...
%       [server_path '2008-06-10-0/data017/data017-jlg.classification.txt']));
%
%
% See also: LOAD_RRS_CELL_TYPES, ORDER_CELL_TYPES,
% LOAD_CELLTYPES_INTO_DATARUN, STRUCT/LOAD_TXT_CELL_TYPES
%
%
% 2008-10 gauthier
%



% SET UP OPTIONAL ARGUMENTS
p = inputParser;

% specify list of optional parameters
p.addParamValue('cell_type_depth', 2, @(x)(x>0 & mod(x,0)==x));
p.addParamValue('order_cell_types', true);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% read in text file
text_file = importdata(file_path,',');

% initialize temporary variable
cell_types_temp = struct;

% cycle through cells and get cell type string
for cc = 1:length(text_file)
    
    % each line of the text file is like this:
    % 578  All/ON/parasol
    %
    % note there are TWO spaces after the cell id


    % get location of the spaces
    space_locs = strfind(text_file{cc},' ');
    
    % get cell id (first characters up to the space
    cell_id = str2num( text_file{cc}(1:strfind(text_file{cc},' ')-1) );
    
    % get cell type of this cell
    type_string = text_file{cc}(space_locs(2)+1:end);
    
    % if there is no class, call it 'unclassified'
    if strcmp(type_string, 'All') || strcmp(type_string, 'All/')
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


% Order?
if params.order_cell_types
    cell_types = order_cell_types(cell_types);
end