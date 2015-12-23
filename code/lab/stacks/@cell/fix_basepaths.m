function stacks = fix_basepaths(stacks)
% FIX_BASEPATHS     Set stack basepaths to standard function handles
% usage: stacks = fix_basepaths(stacks)
%
% Old stack setup had string basepaths, but this was annoying for working
% from different locations.  Makes more sense for basepath to be function
% handle to SERVER_DATA_PATH or the like, which are typically
% customized/shadowed for different locations.
%
% 2013-05, phli
%

% First correct path/basepath for alignment image
[pathstr, name, ext] = fileparts(stacks{1,2}.paths{1});
stacks{1,2}.paths = {fullfile('code', 'projects', 'image alignment', 'images', sprintf('%s%s', name, ext))};
stacks{1,2}.basepath = @matlab_code_path;

% Fix the rest
for r = 2:size(stacks,1)
    for c = 1:size(stacks,2)
        stack = stacks{r,c};
        if isempty(stack),              continue; end
        if ~isfield(stack, 'basepath'), continue; end
        if isempty(stack.basepath),     continue; end
        
        stack.basepath = @server_data_path;
        stacks{r,c} = stack;
    end
end