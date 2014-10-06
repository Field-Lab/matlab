function hout = gen_slider_plot(varargin)
%GEN_SLIDER_PLOT     Basic set up for a plot with a slider; usually called from another more specific function
%
% 2011-11 phli

opts = inputParser();
opts.addParamValue('handle', []);
opts.addParamValue('data', []);
opts.addParamValue('data_length_func', @data_length_func);
opts.addParamValue('slider_func',      @slider_func);
opts.addParamValue('get_imag_func',    @get_imag_func);
opts.addParamValue('plot_func',        @plot_func);
opts.addParamValue('mouse_wheel_func', @mouse_wheel_func);
opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.handle)
    opts.handle = figure();
end
figure(opts.handle);
fig = getfig(opts.handle);

if nargout > 0
    hout = fig;
end


% Add data to figure
setappdata(fig, 'data', opts.data);


% Create an API to simplify calling GUI methods from outside
api = getappdata(fig, 'api');
if isempty(api)
    api = struct();
end
api.data_length = opts.data_length_func;
api.slider_plot = opts.slider_func;
api.get_imag    = opts.get_imag_func;
api.plot        = opts.plot_func;
setappdata(fig, 'api', api);


% Add image panel to fig
impanel = uipanel('Units', 'normalized', 'Position', [0, 0, 1, 1], 'Parent', opts.handle, 'Tag', 'impanel');


% For slider, switch to Metal Look and Feel if running on Mac OS > 10.6
nativelaf = javax.swing.UIManager.getLookAndFeel();
if aftersnowleopard()
    javax.swing.UIManager.setLookAndFeel('javax.swing.plaf.metal.MetalLookAndFeel');
end

% Add slider to panel
start_index = 1;
index_min = 0.99; 
index_max = api.data_length(opts.data);
slider_step = 1 / (index_max - index_min);
slider_steps = min([slider_step 1]) * [1 5];
shandle = uicontrol(impanel, ...
    'Style'     , 'slider',                        ...
    'Min'       , index_min,                       ...
    'Max'       , index_max,                       ...
    'Units'     , 'normalized',                    ...
    'Position'  , [0.01, 0.025, 0.95, 0.01],       ...
    'Value'     , start_index,                     ...
    'SliderStep', slider_steps,                    ...
    'CallBack'  , api.slider_plot,                ...
    'Tag'       , 'slider'                   ...
);
drawnow;
javax.swing.UIManager.setLookAndFeel(nativelaf);
setappdata(fig, 'slider', shandle);

% Set up mouse wheel listener
set(fig, 'WindowScrollWheelFcn', {opts.mouse_wheel_func, shandle});

add_slider_arrowkeys(fig, shandle);

% Set up figure
set(opts.handle, 'Toolbar', 'figure')

% Set up axes in panel
ax = axes('Parent', impanel, 'Units', 'normalized', 'Position', [0, 0.025, 0.95, 0.9]);
axis image;
setappdata(fig, 'slider_axes', ax);

% Initial plot call
api.slider_plot(shandle, []);

% Now set axes to replace children
set(ax, 'NextPlot', 'replaceChildren');



function l = data_length_func(data)
d = ndims(data);
if d < 3
    l = 1;
else
    l = size(data,d);
end


function slider_func(slider, ev) %#ok<INUSD>
f = getfig(slider);
api = getappdata(f, 'api');

% Determine what value the slider is set to
sval = round(get(slider, 'Value'));
pval = round(getappdata(slider, 'PreviousValue'));
if sval == pval, return; end

% Set axes
ax = getappdata(f, 'slider_axes');
axes(ax);

% Plot image data from stack
data = getappdata(f, 'data');
imag = api.get_imag(data, sval);
i = api.plot(imag, sval, data, f);

% Set handle
setappdata(f, 'imhandle', i);

% Useful for higher level guis that add callbacks before this one, e.g. stack_point_picker
setappdata(slider, 'PreviousValue', sval);


function h = plot_func(imag, sval, data, f)
h = imshow(imag);
title(num2str(sval));


function imag = get_imag_func(data, sval)
d = ndims(data);
if d < 4
    imag = data(:,:,sval);
else
    imag = data(:,:,:,sval);
end


function mouse_wheel_func(fig, ev, slider)
% Determine what value the slider is set to
sval = round(get(slider, 'Value'));
newval = sval + ev.VerticalScrollCount/5;

if newval < get(slider, 'Min') || newval > get(slider, 'Max'), return; end
set(slider, 'Value', newval);

api = getappdata(fig, 'api');
api.slider_plot(slider, []);