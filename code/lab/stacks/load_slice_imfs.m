function stack = load_slice_imfs(stack, slice_spec, reload)
% LOAD_SLICE_IMFS    Load image info IMF structs for slices into stack
% usage: stack = load_slice_imfs(stack, slice_spec, reload)
%
% inputs:   stack       Image stack struct as described in Proposal.rtf
%           slice_spec  Slice numbers, or a path, or a cell array of paths
%           reload      Whether imf should be reloaded if already cached.
%                           Defaults to false.
%
% outputs: Stores imf in the stack struct
%
% For retrieving, see also slice_imf slice_fullimf
%
% 2010-08 phli
%

if nargin < 2
    slice_spec = 'all';
end

if nargin < 3
    reload = false;
end

switch class(slice_spec)
    case 'double'
        paths = unique(stack.paths(slice_spec));
    case 'char'
        if strcmp(slice_spec, 'all');
            paths = unique(stack.paths);
        else
            paths = {slice_spec};
        end
    case 'cell'
        paths = unique(slice_spec);
end


% Use a HashMap/Associative_Array to store IMFs by filepath
if ~isfield(stack, 'imfs')
    stack.imfs = containers.Map;
end


for i = 1:length(paths)
    path = paths{i};
    fullpath = add_base_path(path, get_stack_basepath(stack));
    if ~stack.imfs.isKey(path) || reload
        stack.imfs(path) = imfinfo(fullpath);
    end
end