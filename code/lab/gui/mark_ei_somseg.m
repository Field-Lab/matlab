function [fig markpos] = mark_ei_somseg(datarun, cellid, varargin)

opts = inputParser();
opts.addParamValue('marks', []);
opts.addParamValue('defaultarrowpos', [0 0; 50 50]);
opts.addParamValue('scale', 1);
opts.parse(varargin{:});
opts = opts.Results;

fig = plot_ei_scroll(datarun, cellid, 'scale', opts.scale);
[eipicker csdpicker] = add_ei_csd_toggle(fig, datarun, cellid);
[downer upper] = add_ei_scale_spinner(fig, [0.85 0.04 0.08 0.08]);

% Expand figure to fill vertical space while maintaining AR
set(fig, 'Units', 'normalized');
pos = get(fig, 'Position');
pos(4) = 1;
set(fig, 'Position', pos);
set(fig, 'Units', 'pixels');
pos = get(fig, 'Position');
pos(3) = pos(4);
set(fig, 'Position', pos);
axis fill

% Inactivate the slider until arrow is placed.  A kludge.  Alternative easy
% fix is to place the initial arrow position automatically based on the 
% peak electrode and some other easy to calculate point.
slider = getappdata(fig, 'ei_slider');
set([slider eipicker csdpicker downer upper], 'Enable', 'off');

arrowargs = {};
if ~isempty(opts.marks), arrowargs{end+1} = opts.marks; end
pos = opts.defaultarrowpos;
if isempty(pos)
    [arrow pos] = imarrow(gca, arrowargs{:});
elseif ~isempty(arrowargs)
    [arrow pos] = imarrow(gca, arrowargs{:});
else
    arrow = imarrow(gca, pos);
end
arrow.addNewPositionCallback(@(pos)(setappdata(fig, 'arrowpos', pos)));
setappdata(fig, 'arrow', arrow);
setappdata(fig, 'arrowpos', pos);

iptaddcallback(slider,    'Callback', @show_arrow);
iptaddcallback(eipicker,  'Callback', @show_arrow);
iptaddcallback(csdpicker, 'Callback', @show_arrow);
iptaddcallback(downer,    'Callback', @show_arrow);
iptaddcallback(upper,     'Callback', @show_arrow);
set([slider eipicker csdpicker downer upper], 'Enable', 'on');

if nargout < 1, clear fig; end
if nargout > 1
    dataout = waitfordata(fig, 'arrowpos');
    markpos = dataout{1};
end


function show_arrow(h,~)
fig = getfig(h);
figure(fig);
ax = gca();

pos = getappdata(fig, 'arrowpos');
arrow = imarrow(ax, pos);
arrow.addNewPositionCallback(@(pos)(setappdata(fig, 'arrowpos', pos)));

setappdata(fig, 'arrow', arrow);