function hout = stack_point_picker(stack, varargin)
% STACK_POINT_PICKER    Stack slider plot with added point picking features
% usage: handle = stack_point_picker(stack, opts)
%
% opts  handle  figure()    Figure to use for GUI; defaults to new figure
%       points  []          Use to pass in previously picked points.
%                           Should be a cell array holding an Nx2 matrix of
%                           point XY coordinates for each section.
%       points_color   [0 0 1]   RGB color to use for points
%       bgpoint_color  []        Color to use for points that are out of plane
%
% 
% After picking points, save out to the workspace using the pulldown menu.
% Or they are also available via: getappdata(gcf, 'points').  Or use
% WAITFORDATA to block until GUI is closed and get points.
%
% See also: STACK_SLIDER_PLOT, WAITFORDATA
%
% 2010-05 phli
%

opts = inputParser;
opts.addParamValue('handle', []);
opts.addParamValue('points', []);
opts.addParamValue('points_color', [0 0 1]);
opts.addParamValue('bgpoint_color', []);
opts.addParamValue('bgpoint_style', '.');
opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;
handle = opts.handle;

fig = stack_slider_plot(stack, handle, unmatched);
if nargout > 0, hout = fig; end
gh = guihandles(fig);
slider = gh.slider;

% Add to API object
api = getappdata(fig, 'api');
api.record_points           = @record_points;
api.sort_points             = @sort_points;
api.sort_all_points         = @sort_all_points;
api.image_click_passthrough = @image_click_passthrough;
api.call_ax                 = @call_ax;
api.show_points             = @show_points;
api.points_saved            = @points_saved;
setappdata(fig, 'api', api);

% Set up prepicked points
opts.points = sort_all_points(opts.points);
setappdata(fig, 'points', opts.points);
setappdata(fig, 'points_color', opts.points_color);
setappdata(fig, 'bgpoint_color', opts.bgpoint_color);
setappdata(fig, 'bgpoint_style', opts.bgpoint_style);

% Add crosshair pointer and point picking callback
ax = getappdata(fig, 'slider_axes');
add_point_picker(ax);

% Calls to set up initial image
image_click_passthrough(slider, []);
show_points(slider, []);


% Add new callbacks both before and after preexisting slider_plot
slider_plot = get(slider, 'Callback');
set(slider, 'Callback', []);
iptaddcallback(slider, 'Callback', @record_points);
iptaddcallback(slider, 'Callback', slider_plot);
iptaddcallback(slider, 'Callback', @image_click_passthrough); % Have to pass click events through for each image to axes :(
iptaddcallback(slider, 'Callback', @show_points);
api = getappdata(fig, 'api');
api.slider_plot = get(slider, 'Callback');
setappdata(fig, 'api', api);


% Add menu to save points
smenu = uimenu(fig, 'Label', 'Stack', 'Tag', 'smenu');
epmenu = gui_data_export(fig, 'points', smenu);

% Augment the callbacks for point saving menu item
export_points = get(epmenu, 'Callback');
set(epmenu, 'Callback', []);
iptaddcallback(epmenu, 'Callback', api.record_points);
iptaddcallback(epmenu, 'Callback', export_points);
iptaddcallback(epmenu, 'Callback', @points_saved);


% Helpful to have a hotkey for this
add_pan_togglekey(fig, 'control');



% Save points on previous section before plotting new image
function record_points(handle, ev) %#ok<INUSD>
fig = getfig(handle);
gh = guihandles(fig);
slider = gh.slider;
section = round(get(slider, 'Value'));
previous_section = round(getappdata(slider, 'PreviousValue'));

if section == previous_section && handle == slider
    return; 
end

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

function points = sort_all_points(points)
for i = 1:length(points)
    points{i} = sort_points(points{i});
end



% Pass clicks on image through to parent axes
% Not sure why just setting HitTest off didn't take care of this...
%
% ToDo - Starting to use this multiple places.  Abstract?  But would like
% to get to bottom of why this is even necessary...
function image_click_passthrough(slider, ev)
fig = getfig(slider);
im = getappdata(fig, 'imhandle');
set(im, 'ButtonDownFcn', @call_ax);

function call_ax(im, ev)
ax = get(im, 'Parent');
ax_bdf = get(ax, 'ButtonDownFcn');
ax_bdf(ax, ev);




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


function wheel_show_points(fig, ev)
gh = guihandles(fig);
slider = gh.slider;
show_points(slider, ev);



function points_saved(handle, ev) %#ok<INUSD>
fig = getfig(handle);
points = getappdata(fig, 'points');
setappdata(fig, 'saved_points', points);