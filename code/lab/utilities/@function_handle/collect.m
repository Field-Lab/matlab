function varargout = collect(cell_or_struct, fhandle, nout)
% FUNCTION_HANDLE/COLLECT    Run function FHANDLE over each element of CELL_OR_STRUCT and collect outputs in RESULTS_CELL
%
% usage:  results_cell = collect(cell_or_struct, fhandle, nout)
%
% arguments: cell_or_struct - The cell array or struct to operate over 
%            fhandle        - Handle for the function to call on each element.
%                             Function can take one arg or two.  If two, the
%                             second arg passed is the index of each element.
%            nout           - The number of outputs to collect - defaults to
%                             nargout(fhandle) unless fhandle is varargout.
%
% outputs: results_cell - Cell array collecting the results of the function
%                         calls.  If FHANDLE gives multiple outputs, then
%                         RESULTS_CELL will be a cell array of cell arrays.
%                         If there are no outputs per call, then there will
%                         be no outputs at all.
%
% recommendation: If you need to run a function with more complicated
% inputs, suggest wrapping in anonymous function passed into COLLECT.
%
% note: Apparently, function_handle args take precedence over cell args,
% hence this must be in the @function_handle directory...
%
% TODO: This should be changed to rely on cellfun and arrayfun behind the
% scenes and only be a wrapper to work on either cells or structs.
%
% 2010-01 phli
%

if nargin < 3
    nout = nargout(fhandle);
end

% If FHANDLE is varargout, hack around...
if nout == -1
    nout = 1;
end

% Preallocate results cell array
if nout > 0
    varargout{1} = cell(size(cell_or_struct));
else
    varargout = {};
end

for i = 1:numel(cell_or_struct)
    if iscell(cell_or_struct)
        if nargin(fhandle) == 1
            ins = cell_or_struct(i);
        else
            % ToDo: Which is more efficient?
            ins = {cell_or_struct{i}, i};
            % ins = vertcat(cell_or_struct(i), {i})
        end
    elseif isstruct(cell_or_struct)
        ins = {cell_or_struct(i)};
    else 
        error('Invalid input');
    end

    % Return outputs?
    if nout > 0
        outs{1:nout} = fhandle(ins{:}); %#ok<AGROW>

        % Simplify RESULTS_CELL if only one output per call
        if nout == 1
            varargout{1}{i} = outs{1};
        else
            varargout{1}{i} = outs;
        end
    else
        fhandle(ins{:});
    end
end