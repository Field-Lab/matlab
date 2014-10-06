function add_point_picker(ax)
% ADD_POINT_PICKER    Add a point picker to the given axes
% usage: add_point_picker([ax])
%
% Defaults to gca
%
% 2010-05 phli
%

if nargin < 1, ax = gca; end

% Add/enable pointer manager on parent fig
f = getfig(ax);
iptPointerManager(f);

% Add change to cross pointer behavior for ax
pb.enterFcn = @(f, p) set(f, 'Pointer', 'cross');
pb.traverseFcn = [];
pb.exitFcn = [];
iptSetPointerBehavior(ax, pb);

% Add to api
api = getappdata(f, 'api');
api.add_point = @add_point;
setappdata(f, 'api', api);

set(ax, 'ButtonDownFcn', @add_point);




function add_point(ax, ev) %#ok<INUSD>
f = getfig(ax);
points_color = getappdata(f, 'points_color');

% Get click coordinates
cp = get(ax, 'CurrentPoint');
x = cp(1,1);
y = cp(1,2);

% Plot new point
p = impoint(ax, x, y);
if ~isempty(points_color)
    p.setColor(points_color);
end

% Add "Delete" to context menu
make_deleteable(p);

% For older versions of Matlab; should be done automatically on modern versions
set(p, 'Tag', 'impoint');