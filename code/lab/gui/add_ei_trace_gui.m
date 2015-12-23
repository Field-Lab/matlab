function add_ei_trace_gui(fig)
% ADD_EI_TRACE_GUI
% usage: add_ei_trace_gui(fig)
%
% Abstracted out of big messy plot_ei_scroll_.  Now it is optional whether
% this GUI to show ei traces is part of the EI scrolling GUI.
%
% See also PLOT_EI_SCROLL, PLOT_EI_SCROLL_
%
% 2012-08 phli
%

% Add callback for the actual EI electrode plotted markers (run by plot_ei_)
setappdata(fig, 'elec_callback', @toggle_electrode_trace_by_handle);

% Add callback for clicking on space between electrode plot markers
ei_slider = getappdata(fig, 'ei_slider');
iptaddcallback(ei_slider, 'Callback', @show_nearest_electrode_callback);

% Call once for initial state
show_nearest_electrode_callback();

% Add plotting markers for traces currently shown
iptaddcallback(ei_slider, 'Callback', @add_markers_callback);



function show_nearest_electrode_callback(handle, event)
set(gca, 'ButtonDownFcn', @show_nearest_electrode);

function add_markers_callback(handle, event)
fig = getfig(handle);
movie_timer = getappdata(fig, 'movie_timer');
if ~isempty(movie_timer) && strcmp(movie_timer.Running, 'off')
    add_markers(fig, []);
end



function e = get_nearest_electrode(click_coords, electrode_positions)
num_electrodes = size(electrode_positions, 1);
distances = electrode_positions - repmat(click_coords, [num_electrodes 1]);
distances = sqrt(distances(:,1).^2 + distances(:,2).^2);
least_distance = min(distances);
closest = find(distances == least_distance);
e = closest(1);


function show_nearest_electrode(handle, event) %#ok<INUSD>
fig = getfig(handle);
figure(fig);
ax = gca;

click_coords = get(ax, 'CurrentPoint');
click_coords = click_coords(1,1:2); % Drop z dimension
electrode_positions = getappdata(fig, 'positions');
e = get_nearest_electrode(click_coords, electrode_positions);
toggle_electrode_trace(fig, e);


% Wrapper for toggle_electrode_trace
function toggle_electrode_trace_by_handle(electrode, event) %#ok<INUSD>
fig = getfig(electrode);
num = getappdata(electrode, 'num');
toggle_electrode_trace(fig, num);


function toggle_electrode_trace(fig, num)
api = getappdata(fig, 'api');
api.stop_movie(fig, []);

trace_plot = getappdata(fig, 'trace_plot');
if isempty(trace_plot) || ~ishandle(trace_plot)
    trace_plot = new_trace_plot(fig);
end
figure(trace_plot);

% If this electrode is already plotted, delete it
traces = getappdata(trace_plot, 'traces');
markers = getappdata(fig, 'markers');
if length(traces) >= num;
    old_trace = traces(num);
    if old_trace ~= 0 && ishandle(old_trace)
        disp(['Deleting electrode trace ' num2str(num)]);
%        axis manual;
        delete(old_trace);
        delete(markers(num));
        markers(num) = 0;
        setappdata(fig, 'markers', markers);
        return
    end
end

% Otherwise plot it
plot_ei_trace(trace_plot, num);


function trace_plot = new_trace_plot(fig)
trace_plot = figure;
hold on;
grid on;

% Move the trace plot to the right of the ei plot.
ei_plot_position = get(fig, 'Position');
ei_right_edge = ei_plot_position(1) + ei_plot_position(3);
set(trace_plot, 'Units', get(fig, 'Units'));
trace_plot_position = get(trace_plot, 'Position');
trace_plot_position(1) = ei_right_edge + (ei_plot_position(3) / 40);
set(trace_plot, 'Position', trace_plot_position);

data = getappdata(fig, 'data');
setappdata(trace_plot, 'data', data);

setappdata(fig, 'trace_plot', trace_plot);
setappdata(trace_plot, 'ei_plot', fig);

% When closing the trace plot, remove electrode markers from ei_plot
set(trace_plot, 'CloseRequestFcn', @remove_markers);
iptaddcallback(trace_plot, 'CloseRequestFcn', 'closereq');


function remove_markers(trace_plot, event) %#ok<INUSD>
ei_plot = getappdata(trace_plot, 'ei_plot');
if ~ishandle(ei_plot), return; end

markers = getappdata(ei_plot, 'markers');
for marker = markers
    if marker ~= 0 && ishandle(marker)
        delete(marker);
    end
end
setappdata(ei_plot, 'markers', []);


function plot_ei_trace(trace_plot, num)
disp(['Plotting electrode trace ' num2str(num)]);
if strcmp(get(gca, 'XLimMode'), 'manual')
    axis auto;
end
ei = getappdata(trace_plot, 'data');
color = next_color;
h = plot(ei(num,:));
set(h, 'Color', color);

traces = getappdata(trace_plot, 'traces');
traces(num) = h;
setappdata(trace_plot, 'traces', traces);


ei_plot = getappdata(trace_plot, 'ei_plot');
marker = plot_marker(ei_plot, num, color);
markers = getappdata(ei_plot, 'markers');
markers(num) = marker;
setappdata(ei_plot, 'markers', markers);
marker_colors = getappdata(ei_plot, 'marker_colors');
marker_colors(num, :) = color;
setappdata(ei_plot, 'marker_colors', marker_colors);


function marker = plot_marker(ei_plot, num, color)
figure(ei_plot);
electrode_positions = getappdata(ei_plot, 'positions');
marker = plot(electrode_positions(num, 1), electrode_positions(num, 2), 'ks');
set(marker, 'MarkerFaceColor', color);
set(marker, 'MarkerSize', 15);
setappdata(marker, 'num', num);
set(marker, 'ButtonDownFcn', @toggle_electrode_trace_by_handle);


function add_markers(handle, event) %#ok<INUSD>
fig = getfig(handle);
electrode_positions = getappdata(fig, 'positions');
markers = getappdata(fig, 'markers');
marker_colors = getappdata(fig, 'marker_colors');
for i = 1:length(markers)
    marker = markers(i);
    if marker == 0
        continue
    end
    color = marker_colors(i, :);
    markers(i) = plot_marker(fig, i, color);
end
setappdata(fig, 'markers', markers);