function [result, found] = detect(cell_arr, fhandle)
% FUNCTION_HANDLE/DETECT    Find element of CELL_ARR for which FHANDLE evaluates to TRUE
%
% usage:  result = detect(cell_arr, fhandle)
%
% arguments: cell_arr - The cell array to operate over 
%            fhandle  - Handle for the function to call on each element.
%                       Function can take one arg or two.  If two, the
%                       second arg passed is the index of each element.
%
% outputs: result - The first element of CELL_ARR for which FHANDLE
%                   evaluates to true.  If no element is found, returns []
%          found - Not generally needed; only useful if it is possible that
%                  the detected element is in fact [].  In this case, FOUND
%                  evaluates to true, otherwise false.
%
% recommendation: If you need to run a function with more complicated
% inputs, suggest wrapping in anonymous function passed into COLLECT.
%
% note: Apparently, function_handle args take precedence over cell args,
% hence this must be in the @function_handle directory...
%
% 2010-01 phli

for i = 1:numel(cell_arr)
    if nargin(fhandle) == 1
        ins = cell_arr(i);
    else
        % ToDo: Which is more efficient?  
        ins = {cell_arr{i}, i}; 
        % ins = vertcat(cell_arr(i), {i})
    end
    
    if fhandle(ins{:})
        result = cell_arr{i};
        found = true;
        return;
    end
end

result = [];
found = false;