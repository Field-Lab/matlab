function results = select(cell_or_struct, fhandle)
% FUNCTION_HANDLE/SELECT    Select elements from CELL_OR_STRUCT for which FHANDLE evaluates to TRUE
%
% usage:  results = select(cell_or_struct, fhandle)
%
% arguments: cell_or_struct - The cell array or struct to operate over 
%            fhandle  - Handle for the function to call on each element.
%                       Function can take one arg or two.  If two, the
%                       second arg passed is the index of each element.
%
% outputs: results - Subset of elements from CELL_OR_STRUCT
%
% note: Apparently, function_handle args take precedence over cell args,
% hence this must be in the @function_handle directory...
%
% 2011-01 phli
%

selected = false(size(cell_or_struct));
for i = 1:numel(cell_or_struct)
    if iscell(cell_or_struct)
        elem = cell_or_struct{i};
    elseif isstruct(cell_or_struct)
        elem = cell_or_struct(i);
    else
        error('Invalid input');
    end
    
    if nargin(fhandle) == 1
        selected(i) = fhandle(elem);
    elseif nargin(fhandle) == 2
        selected(i) = fhandle(elem, i);
    end
end

results = cell_or_struct(selected);