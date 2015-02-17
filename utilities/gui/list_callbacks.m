function cbstruct = list_callbacks(handle, callback_name)
% LIST_CALLBACKS    Get a struct of the functions set up on a given callback
% usage: cbstruct = list_callbacks(handle, callback_name)
%
% Standard MatLab is to set GUI callbacks to single functions, but the
% Image Processing Toolbox adds to this functionality with the
% IPTADDCALLBACK and IPTREMOVECALLBACK methods, which allow loading a
% series of functions onto a single callback.
%
% LIST_CALLBACKS gets the list of callback functions, assuming
% IPTADDCALLBACK has been used.  Also returns the ids for each callback,
% which can be used with IPTREMOVECALLBACK.
%
% NOTE: This function uses undocumented/unstable MatLab internal functions
% to pull the function list from iptaddcallback/callbackProcessor.
% Therefore, this may break if new MatLab releases change the backend.
%
% 2010-05 phli
%

if nargin < 2
    callback_name = 'Callback';
end

callback_obj = get(handle, callback_name);

if isempty(callback_obj)
    cbstruct = [];
    return;
end

if ~isa(callback_obj, 'function_handle')
    warning('LIST_CALLBACKS:nofunc', 'Not a function handle property');
    cbstruct = [];
    return;
end

if strcmp(func2str(callback_obj), 'iptaddcallback/callbackProcessor')
    funcs = functions(callback_obj);
    cbstruct = funcs.workspace{1}.callback_list;
else
    cbstruct = struct('func', callback_obj, 'id', []);
end