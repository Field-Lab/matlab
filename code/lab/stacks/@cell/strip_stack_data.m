function stacks = strip_stack_data(stacks)
% STRIP_STACK_DATA      Remove raw data from stacks as long as path to data files exists
% usage: stacks strip_stack_data(stacks)
%
% This is mostly used during SAVE_STACKS to avoid saving large raw data
% that can easily be reloaded from data files on disk.
%
% 2010-09 phli
%

for i = 1:numel(stacks)
    stack = stacks{i};

    if ~isstruct(stack) || ~isfield(stack, 'data')
        continue
    end
    
    % If there is a path to the raw data, strip the data here; it can be
    % reloaded
    for j = 1:length(stack.data)
        if ~isempty(stack.data{j}) && ~isempty(stack.paths{j})
            stack.data{j} = [];
        end
    end
    
    if all(stack.data, @isempty)
        stack.data = {};
    end
    
    stacks{i} = stack;
end