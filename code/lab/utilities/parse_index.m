function index = parse_index(varargin)
% PARSE_INDEX    Flexible indexing utility
% usage: index = parse_index(1, 2, 3)
% usage: index = parse_index([1 2 3])
% usage: index = parse_index({1 2 3})
%
% All three forms return {1 2 3}.  Cell arrays are useful for indexing
% operations because you can do: a = {1 2}; disp(b{a{:}}).  But cell arrays
% can be ugly/clumsy to type and pass around, so it's nice to have indexing
% methods pass through here.
%
% 2010-08 phli
%

index = {};
if length(varargin) == 1
    stack_spec = varargin{1};
    
    switch class(stack_spec)
        case 'cell'
            index = stack_spec;

        case 'double'
            index = num2cell(stack_spec);
    end
else
    index = varargin;
end