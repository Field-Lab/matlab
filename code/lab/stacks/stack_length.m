function lngth = stack_length(stack)
% STACK_LENGTH
% usage: lngth = stack_length(stack)
%
% 2010-08 phli
%

lngth = 0;

if isempty(stack)
    return
end

if isfield(stack, 'paths')
    lngth = length(stack.paths);
end

if isfield(stack, 'data') && length(stack.data) > lngth
    lngth = length(stack.data);
end