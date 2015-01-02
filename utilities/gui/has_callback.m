function has = has_callback(handle, func, callback_name)
% HAS_CALLBACK    Check whether the given callback already includes the given func
% usage: has = has_callback(handle, cbfunc, [callback_name])
%
% Standard MatLab is to set GUI callbacks to single functions, but the
% Image Processing Toolbox adds to this functionality with the
% IPTADDCALLBACK and IPTREMOVECALLBACK methods, which allow loading a
% series of functions onto a single callback.
%
% HAS_CALLBACK checks whether the series of functions assigned to the given
% callback includes the given function.
%
% NOTE: This function uses undocumented/unstable MatLab internal functions
% to pull the function list from iptaddcallback/callbackProcessor.
% Therefore, this may break if new MatLab releases change the backend.
%
% 2010-05 phli
%

if nargin < 3
    callback_name = 'Callback';
end

cbstruct = list_callbacks(handle, callback_name);
for i = 1:length(cbstruct)
    cbfunc = cbstruct(i).func;
    if cbfunc == func
        has = true;
        return;
    end
end

has = false;