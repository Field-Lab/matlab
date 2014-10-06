function index = parse_stack_index(varargin)
% PARSE_STACK INDEX     Convert various stack spec styles into cell indices
% usage: index = parse_stack_index(2, 2)
%        index = parse_stack_index([2 2])
%        index = parse_stack_index('alive_montage_hires')
% all these generate index == {2,2}.  Named indices are converted based on
% valued in STACK_LABELS
%
% See also: STACK_LABELS
%
% 2010-09 phli
%

index = {};
if length(varargin) == 1 && ischar(varargin{1})
    stack_spec = varargin{1};
    
    cell_labels = stack_labels();
    if isfield(cell_labels, stack_spec)
        index = cell_labels.(stack_spec);
    end
else    
    index = parse_index(varargin{:});
end