function bool = any(cell_arr, fhandle)
% FUNCTION_HANDLE/ANY    Pass elements of CELL_ARR to FHANDLE, determine if any evaluates to true
%
% usage:  if any(cell_arr, fhandle) ...
%
% arguments: cell_arr - The cell array to operate over 
%            fhandle  - Handle for the function to call on each element.
%                       Function can take one arg or two.  If two, the
%                       second arg passed is the index of each element.
%
% outputs: bool - LOGICAL value indicating if any call evaluated to true
%
% recommendation: If you need to run a function with more complicated
% inputs, suggest wrapping in anonymous function passed into COLLECT.
%
% note: Apparently, function_handle args take precedence over cell args,
% hence this must be in the @function_handle directory...
%
% 2010-01 phli

[dummy, bool] = detect(cell_arr, fhandle);