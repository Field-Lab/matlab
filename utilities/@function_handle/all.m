function bool = all(cell_arr, fhandle)
% FUNCTION_HANDLE/ALL    Evalute FHANDLE for each element of CELL_ARR and determine if all calls return true
%
% usage:  if all(cell_arr, fhandle) ...
%
% arguments: cell_arr - The cell array to operate over 
%            fhandle  - Handle for the function to call on each element.
%                       Function can take one arg or two.  If two, the
%                       second arg passed is the index of each element.
%
% outputs: bool - LOGICAL value indicating if all calls evaluated to true
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
    
    if ~fhandle(ins{:})
        bool = false;
        return;
    end
end

bool = true;