function datarun = load_vision_cell_types(datarun, varargin)
% LOAD_CELL_TYPES    Load cell types from the classification.txt file
%
% Provenance not clear, redundant code should be abstracted.
%
% See also: LOAD_TXT_CELL_TYPES, ORDER_CELL_TYPES, LOAD_CELLTYPES_INTO_DATARUN
%

p = inputParser;
% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('sync_cell_ids', true);
p.addParamValue('cell_type_depth', 3);
p.addParamValue('order_cell_types', true);
p.addParamValue('load_into_workinglist', true);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% check arguments
if ~isfield(datarun.names,'rrs_classification_path')
    error('load_vision_cell_types: no rrs_classification_path');
end

if params.verbose;
    fprintf('\nLoading cell types...');
    start_loading = clock; % note when it started
end

% open and load file
fid=fopen(datarun.names.rrs_classification_path);
if (fid == -1)
    error('load_vision_cell_types: file not found');
end

% read entire Index file
file=textscan(fid, '%s', 'Delimiter', '\t');
fclose(fid);
textfile=file{1};

cell_types_temp = struct;
for i=1:length(textfile)
    
    t=strfind(textfile{i},' ');
    cell_id=str2num(textfile{i}(1:t(1)-1));
    type_string=textfile{i}(t(2)+1:end);
    
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


% order cell types, if desired
if params.order_cell_types
    cell_types = order_cell_types(cell_types);
end

% save them in datarun
datarun.vision.cell_types = cell_types;

% load these cell types into datarun.cell_types if none have been defined so far
if params.load_into_workinglist
    datarun.cell_types = datarun.vision.cell_types;
    fprintf('These cell types are loaded into the working list of cell types\n\n')
else
    fprintf('These cell types were NOT loaded into the working list of cell types\n\n')
end

% display how long it took
if params.verbose;
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
end

% display which cell types were loaded
if params.verbose
    fprintf('\nLoaded %d cell types:\n',length(datarun.cell_types))
    show_cell_types(datarun.cell_types)
    fprintf('\n')
end
