function varargout = try_warn(fhandle, args, errid, message, default_output)
% TRY_WARN    Try a function call, converting expected errors into mere warnings
%
% usage: try_warn(fhandle, args, errid, message [, default_output])
%
% arguments: fhandle          - Handle for the function to call
%            args             - Cell array of the arguments to pass to fhandle function
%            errid            - The error id to catch and covert into a warning
%            message          - The message to put in the warning
%            [default_output] - What to return if the function call throws error
%
% example: try_warn(@sin, {1, 2, 3}, 'MATLAB:maxrhs', 'Too many input arguments.  Let this be a warning to you!')
%
% 2010-01 phli
if nargin < 5
    default_output = [];
end

% Determine how many outputs to give
nout = nargout;

try
    % Try running the function call ...
    varargout{1:nout} = fhandle(args{:});
catch err
    % If there was an error matching ERRID, then convert it to a warning
    if strcmp(err.identifier, errid) || strcmp(errid, 'all')
        
        % Show a warning, but not more than once consecutively
        [lastwarnmsg, lastwarnid] = lastwarn();
        if ~strcmp(lastwarnid, err.identifier)
            fprintf(1, '\n');
            warning(err.identifier, message);
        end
        
        % The function call didn't run successfully, but output is still
        % expected, so come up with something
        if iscell(default_output) && size(default_output(:), 1) == nout
            % If the user provided a cell array with the right number of
            % elements as DEFAULT_OUTPUT, use that.
            varargout = default_output(:);
        else
            % Default to assigning the same thing to all outputs
            [varargout{1:nout}] = deal(default_output);
        end

    else
        % Not an error we expected, so pass it along
        rethrow(err)
    end
end