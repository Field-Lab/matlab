function [exit_val, ret_val] = nbinfo(varargin)

if nargout > 0
    [exit_val, ret_val] = system([scripts_path 'nbinfo ' join(varargin, ' ')]);
else
    system([scripts_path 'nbinfo ' join(varargin, ' ')]);
end