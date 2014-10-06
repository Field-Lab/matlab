function [edges_x, edges_y] = map2manhattan(stimmap, varargin)

opts = inputParser();
opts.addParamValue('collapse_output', false);
opts.parse(varargin{:});
opts = opts.Results;

if ischar(stimmap)
    stimmap = load_stimmap(stimmap);
end

indices = unique(stimmap);
indices = setdiff(indices, 0);
indices = sort(indices);
numindices = length(indices);

edges_x = cell(1,numindices);
edges_y = cell(1,numindices);
for i = indices(:)'
    [edges_x{i},edges_y{i}] = mask2manhattan(stimmap == i);
end

if opts.collapse_output
    edges_x = {cell2mat(edges_x)};
    edges_y = {cell2mat(edges_y)};
end


% Can be expanded to load based on parsing either stimulus.lisp or the lisp
% stimulus output (see read_stim_lisp_output);
function map = load_stimmap(filename)
map = dlmread(filename);