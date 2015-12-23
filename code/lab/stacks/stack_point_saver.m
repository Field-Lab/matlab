function varargout = stack_point_saver(stack, varargin)
% STACK_POINT_SAVER    Stack point picker that complains if you try to exit without saving the points
% usage: handle = stack_point_saver(stack, opts)
%
% See STACK_POINT_PICKER for opts and details.
%
% NOTE: GUI_SAVE_BEFORE_CLOSE relies on HAS_CALLBACK, which uses
% undocumented/unstable MatLab internal functions. This may break if new
% MatLab releases change the backend.
%
% 2010-05 phli
%

h = stack_point_picker(stack, varargin{:});
if nargout > 0
    varargout{1} = h;
end

fig = getfig(h);
points = getappdata(fig, 'points');
setappdata(fig, 'saved_points', points);

% Record current points before checking whether points need saving
api = getappdata(fig, 'api');
set(fig, 'CloseRequestFcn', api.record_points);

% Warn if closing without saving points
gui_save_before_close(fig, {'points'});
