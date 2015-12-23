function stack = get_stack(datarun, varargin)
% GET_STACK     Flexible indexing into stacks cell array
% usage: stack = get_stack(datarun, {2,1})
%        stack = get_stack(datarun, [2,1])
%        stack = get_stack(datarun, 2, 1)
%        stack = get_stack(datarun, 'alive_montage_lores')
%
% 2010-08 phli
%

stack = get_stack(datarun.stacks, varargin{:});