function stack = build_stack(varargin)
% BUILD_STACK    Helper for building image stack structs
% usage: stack = build_stack(paths, field_name1, field_value1, field_name2, field_value2, ...)
% usage: stack = build_stack(field_name1, field_value1, field_name2, field_value2, ...)
%
% PATHS is a single path string or a cell array of path strings.  The
% strings can be formatted as, e.g.: '/image1.tif:3' to indicate that the stack
% is from a multi-image tif and this slice is page 3 of the tif.  This will
% be parsed into paths and pages fields automatically.  The path strings
% will also be run through CURLSTR().
%
% See also CURLSTR, GET_STACK_BASEPATH
%
% 2010-08 phli
%

paths = {};

% Odd number of args?
if mod(length(varargin), 2) == 1
    paths = varargin{1};
    varargin = varargin(2:end); 
end


opts = inputParser;
opts.addParamValue('name', '');
opts.addParamValue('desc', '');
opts.addParamValue('paths', {});
opts.addParamValue('pages', {});
opts.addParamValue('data',  {});
opts.addParamValue('basepath', []);
opts.parse(varargin{:});
opts = opts.Results;

if isempty(paths) && isempty(opts.paths) && isempty(opts.data)
    stack = [];
    return
else
    stack = struct(opts);
end


if ~isempty(paths)
    % If single string given, put it in a cell array for processing
    if ischar(paths)
        paths = {paths};
    end
    
    if isempty(stack.basepath) && ~strcmp(paths{1}(1), '/') && ~is_windows_root_path(paths{1})
        stack.basepath = @server_data_path;
    end
    
    % Run curlstr on all given strings; will do nothing if no [] are present
    expanded_paths = {};
    for i = 1:length(paths)
        curlpath = curlstr(paths{i});
        expanded_paths = {expanded_paths{:} curlpath{:}};
    end
else
    expanded_paths = opts.paths;
end


% Get basepath for below
basepath = get_stack_basepath(stack);


% Parse ":" notation for multipage tifs
stack.paths = cell(size(expanded_paths));
stack.pages = cell(size(expanded_paths));
for i = 1:length(expanded_paths)
    pathstr = expanded_paths{i};
    [stack.paths{i}, stack.pages{i}] = parse_pathstr(pathstr);

    % This indicates to load all pages
    if stack.pages{i} == 0
       inf = imfinfo([basepath stack.paths{i}]);
       stack.paths{i} = repmat(stack.paths(i), [1 length(inf)]);
       stack.pages{i} = num2cell(1:length(inf));
    end
end
stack.paths = flatten(stack.paths);
stack.pages = flatten(stack.pages);


% Verify
uniquepaths = unique(stack.paths);
missing = false(1,length(uniquepaths));
for i = 1:length(uniquepaths)
    if ~exist([basepath uniquepaths{i}], 'file')
        missing(i) = true;
    end
end
if any(missing)
    warning(['Files not found for paths: ' basepath join(uniquepaths(missing), [', ' basepath])]);
end