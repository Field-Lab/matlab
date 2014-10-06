function results = invertselect(cell_or_struct, fhandle)
% FUNCTION_HANDLE/INVERTSELECT    Select elements from CELL_OR_STRUCT for which FHANDLE evaluates to FALSE
% usage:  results = invertselect(cell_or_struct, fhandle)
%
% 2012-03 phli
%

selected = true(size(cell_or_struct));
for i = 1:numel(cell_or_struct)
    if iscell(cell_or_struct)
        elem = cell_or_struct{i};
    elseif isstruct(cell_or_struct)
        elem = cell_or_struct(i);
    else
        error('Invalid input');
    end
    
    if nargin(fhandle) == 1
        selected(i) = ~fhandle(elem);
    elseif nargin(fhandle) == 2
        selected(i) = ~fhandle(elem, i);
    end
end

results = cell_or_struct(selected);