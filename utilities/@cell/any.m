function bool = any(cell_arr)
% CELL/ANY    Check whether any elements of cell_arr evaluate to true
% usage:  if any(cell_arr)
%
% arguments: cell_arr - The cell array to evaluate over
%
% outputs: bool - LOGICAL value indicating if any call evaluated to true
%
% 2010-09 phli
%

bool = any(cell_arr, @evaltrue);


function bool = evaltrue(x)

bool = false;
if isnumeric(x)
    if x
        bool = true;
    end
else
    if ~isempty(x)
        bool = true;
    end
end
