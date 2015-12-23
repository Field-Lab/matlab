function fig = im_point_picker(im, varargin)
% IM_POINT_PICKER    Simple image plot with added point picking features
% usage: handle = im_point_picker(image, opts)
%
% opts  handle  figure()    Figure to use for GUI; defaults to new figure
%       points  []          Use to pass in previously picked points.
%                           Should be an Nx2 matrix of point XY coordinates.
%       points_color   [0 0 1]   RGB color to use for points
% 
% After picking points, save out to the workspace using the pulldown menu.
% Or they are also available via: getappdata(gcf, 'points').  Or use
% WAITFORDATA.
%
% 2012-09 phli
%

opts = inputParser;
opts.addParamValue('handle', gcf);
opts.addParamValue('points', []);
opts.addParamValue('points_color', [0 0 1]);
opts.parse(varargin{:});
opts = opts.Results;

handle = opts.handle;
fig = figure(handle);
ax = gca;

% Takes a lot of work to get image shown, axes/fig restored to right size,
% and clicks passed through to axes!
fpos = get(fig, 'Position');
axpos = get(ax, 'Position');
hi = imshow(im);
set(hi, 'HitTest', 'off');
set(ax, 'Position', axpos, 'Visible', 'on');
set(fig, 'Position', fpos);
gh = guihandles(fig);

% Add to API object
api = getappdata(fig, 'api');
api.record_points           = @record_points;
api.sort_points             = @sort_points;
api.show_points             = @show_points;
api.points_saved            = @points_saved;
setappdata(fig, 'api', api);

% Set up prepicked points
opts.points = sort_points(opts.points);
setappdata(fig, 'points', opts.points);
setappdata(fig, 'points_color', opts.points_color);

% Add crosshair pointer and point picking callback
add_point_picker(ax);

% Calls to set up initial image
show_points(ax, []);


% Add menu to save points
smenu = uimenu(fig, 'Label', 'Stack', 'Tag', 'smenu');
epmenu = gui_data_export(fig, 'points', smenu);

% Augment the callbacks for point saving menu item
export_points = get(epmenu, 'Callback');
set(epmenu, 'Callback', []);
iptaddcallback(epmenu, 'Callback', @record_points);
iptaddcallback(epmenu, 'Callback', export_points);
iptaddcallback(epmenu, 'Callback', @points_saved);

if nargout < 1, clear fig; end




% Save points on previous section before plotting new image
function record_points(handle, ev) %#ok<INUSD>
fig = getfig(handle);
gh = guihandles(fig);

section_points = [];
ax = getappdata(fig, 'slider_axes');
children = get(ax, 'Children');
for i = 1:length(children)
    child = children(i);
    if ~strcmp(get(child, 'Tag'), 'impoint'), continue; end
    
    symbols = get(child, 'Children');
    x = get(symbols(2), 'XData');
    y = get(symbols(2), 'YData');

    section_points(end+1,1:2) = [x,y];
end
section_points = sort_points(section_points);

points = getappdata(fig, 'points');
if isempty(points), points = {}; end
points{previous_section} = section_points;
setappdata(fig, 'points', points);


% Allows checking whether points have changed before closing figure
function points = sort_points(points)
if isempty(points), return; end
[~,i] = sort(points(:,1));
points = points(i,:);



function show_points(slider, ev) %#ok<INUSD>

% Check whether we already plotted points for this section (handles
% situation where GUI triggers a bunch of events while the slider is still
% moving...)
section = round(get(slider, 'Value'));
plotted_points = getappdata(slider, 'plotted_points');
if plotted_points == section, return; end
setappdata(slider, 'plotted_points', section);

fig = getfig(slider);
ax = getappdata(fig, 'slider_axes');
points = getappdata(fig, 'points');
if isempty(points) || length(points) < section, return; end


% Plot in-plane impoints
points_color = getappdata(fig, 'points_color');
section_points = points{section};
for i = 1:size(section_points, 1)
    x = section_points(i, 1);
    y = section_points(i, 2);
    p = impoint(ax, x, y);
    p.setColor(points_color);

    % Add "Delete" to context menu
    make_deleteable(p);
    
    % For older versions of Matlab; done automatically on modern versions
    set(p, 'Tag', 'impoint');
end


% Plot out-of-plane simple points?
bgpoint_color = getappdata(fig, 'bgpoint_color');
if isempty(bgpoint_color), return; end

x = cell(length(points),1);
y = cell(length(points),1);
axes(ax);
hold on;
for i = 1:length(points)
    if i == section, continue; end
    section_points = points{i};
    if isempty(section_points), continue; end
    
    x{i} = section_points(:, 1);
    y{i} = section_points(:, 2);
end
x = cell2mat(x);
y = cell2mat(y);

bgpoint_style = getappdata(fig, 'bgpoint_style');
p = plot(x,y,bgpoint_style);
set(p, 'Color', bgpoint_color);



function points_saved(handle, ev) %#ok<INUSD>
fig = getfig(handle);
points = getappdata(fig, 'points');
setappdata(fig, 'saved_points', points);