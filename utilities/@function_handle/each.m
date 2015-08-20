function each(cell_arr, fhandle)
% FUNCTION_HANDLE/EACH    Run function FHANDLE over each element of CELL_ARR, discarding outputs
%
% usage:  each(cell_arr, fhandle)
%
% arguments: cell_arr - The cell array to operate over 
%            fhandle  - Handle for the function to call on each element.
%                       Function can take one arg or two.  If two, the
%                       second arg passed is the index of each element.
%
% recommendation: If you need to run a function with more complicated
% inputs, suggest wrapping in anonymous function passed into COLLECT.
%
% note: Apparently, function_handle args take precedence over cell args,
% hence this must be in the @function_handle directory...
%
% 2010-01 phli

nout = 0;
collect(cell_arr, fhandle, nout);