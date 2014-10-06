function varargout = plot_ei_scroll_(ei, positions, varargin)
% plot_ei_scroll_     plot an ei with a scroll bar
%
% usage:  fig = plot_ei_scroll_(ei, positions, <params>)
%
% arguments:       ei - ExF matrix, E = # electrodes, F = # frames
%           positions - Ex2 matrix of electrode coordinates
%            <params> - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% spike_frame           []      frame of the spike time
% start_frame           []      first frame to plot
% prefix                ''      text to prepend to figure title
% ei_trace_gui          false   whether to add a bunch of infrastructure to allow clicking to view individual traces
%
%
% all other parameters are passed to plot_ei_
%
% by default, these parameters are passed to plot_ei_:
%
%       cutoff              0.05
%       absolute_cutoff     false
%       neg_color           [0 0 1]
%       pos_color           [1 0 0]
%       disk_points         15
%       
%
%
%
% 2010-03  gauthier
% 2010-06  phli - All hacked up to add movie playing and other interface
%                 stuff.  Could use cleaning...
% 2012-08  phli - Finally cleaned out that stuff and abstraced it off to ADD_EI_TRACE_GUI
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% some provided arguments might be for the low level function and should be ignored.
p.KeepUnmatched = true;

% specify list of optional parameters
p.addParamValue('start_frame', []);
p.addParamValue('spike_frame', []);
p.addParamValue('prefix', '');
p.addParamValue('figure', 0);
p.addParamValue('scale', 3);
p.addParamValue('axon', []);
p.addParamValue('sampling_rate', []);
p.addParamValue('add_slider_plot_callbacks', {});
p.addParamValue('ei_trace_gui', false);
p.addParamValue('output_data', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% collect parameters that are ostensibly arguments to plot_ei_
params.to_pass = p.Unmatched;


% add comma to prefix if nonempty
if ~isempty(params.prefix)
    params.prefix = [params.prefix ', '];
end


% if not specified, set start frame to be spike frame
if isempty(params.start_frame) && ~isempty(params.spike_frame)
    params.start_frame = params.spike_frame;
end

% Default to 1
if isempty(params.start_frame), params.start_frame = 1; end

% make figure
set_up_fig_or_axes(params.figure);
params.plot_fig = gcf;
setappdata(params.plot_fig, 'ei', ei);
setappdata(params.plot_fig, 'data', ei);
setappdata(params.plot_fig, 'positions', positions);
setappdata(params.plot_fig, 'params', params);
setappdata(params.plot_fig, 'elec_spacing', infer_electrode_spacing(positions));

set(params.plot_fig, 'CloseRequestFcn', @stop_movie);
iptaddcallback(params.plot_fig, 'CloseRequestFcn', 'closereq');

% create slider control
ha = make_loop_slider_list(params.start_frame, 1, size(ei,2));
iptaddcallback(ha, 'Callback', @stop_movie);
iptaddcallback(ha, 'Callback', @slider_plot);
setappdata(params.plot_fig, 'ei_slider', ha);


% Add EI trace GUI?
if params.ei_trace_gui, add_ei_trace_gui(params.plot_fig); end


% plot once before any clicks
slider_plot(ha, []);


api = struct;
api.slider_plot   = @slider_plot;
api.advance_frame = @advance_frame;
api.start_stop    = @start_stop;
api.stop_movie    = @stop_movie;
setappdata(params.plot_fig, 'api', api);

movie_timer = timer('TimerFcn', @advance_frame, 'Period', 0.05, 'ExecutionMode', 'fixedRate');
movie_timer.UserData = struct('fig', params.plot_fig);
setappdata(params.plot_fig, 'movie_timer', movie_timer);
eimenu = uimenu(params.plot_fig, 'Label', 'EI', 'Tag', 'eimenu');
psmenu = uimenu(eimenu, 'Label', 'Play/Stop', 'Tag', 'psmenu', 'Callback', @start_stop);

if nargout > 0
    varargout{1} = params.plot_fig;
end


% Some tricks to get data out of the figure as it closes
if ~isempty(params.output_data) && nargout > 1
    dataout = waitfordata(params.plot_fig, params.output_data);
    varargout = [varargout dataout];
end
    

function advance_frame(timer, event) %#ok<INUSD>
ud = timer.UserData;
fig = ud.fig;

slider = getappdata(fig, 'ei_slider');
val = get(slider, 'Value') + 1;
if val > get(slider, 'Max')
    val = 1;
end
set(slider, 'Value', val);

api = getappdata(fig, 'api');
api.slider_plot(slider, []);


function stop_movie(handle, event) %#ok<INUSD>
fig = getfig(handle);
movie_timer = getappdata(fig, 'movie_timer');
stop(movie_timer);


function start_stop(handle, event) %#ok<INUSD>
fig = getfig(handle);
movie_timer = getappdata(fig, 'movie_timer');
if strcmp(movie_timer.Running, 'off')
    start(movie_timer);
else
    stop(movie_timer);
end


% display one frame of the EI
function slider_plot(handle, event) %#ok<INUSD>
fig = getfig(handle);
data      = getappdata(fig, 'data'); % Usually the EI, but can be changed
positions = getappdata(fig, 'positions');
params    = getappdata(fig, 'params');

% get the slider position
frame = round(get(handle,'Value'));

% select the figure
figure(params.plot_fig);

% This is weird, but basically the figure can be closed by other callbacks
% in the middle of slider_plot.  slider_plot is only interruptible when it
% calls figure() just above, so right after that call is the time to check
% that the figure hasn't in the meantime been closed!
if ~ishandle(params.plot_fig), return; end


% note xlim, ylim

% if this is the first time, note it, and don't get the axis limits
if all([xlim ylim]==[0 1 0 1])
    first_time = true;
else
    % if this is NOT the first time, note it, and get the axis limits
    curr_lim = [xlim; ylim];
    first_time = false;
end


% plot ei
params.to_pass.elec_spacing  = getappdata(fig, 'elec_spacing');
params.to_pass.elec_callback = getappdata(fig, 'elec_callback');
params.to_pass.scale = params.scale;
params.to_pass.axon  = params.axon;
cla;
plot_ei_(data,positions,frame, 'cutoff', 0.05, 'absolute_cutoff', false, 'neg_color', [0 0 1], 'pos_color', [1 0 0], params.to_pass);


% if this was NOT the first time, restore the axis limits to their previous values
if ~first_time
    xlim(curr_lim(1,:));
    ylim(curr_lim(2,:));
end


% put prefix, frame number into title text
title_text = sprintf('%sframe %d of %d',params.prefix,frame,size(data,2));

% also include timing relative to the spike, if available
if ~isempty(params.spike_frame)
    if ~isempty(params.sampling_rate)
        title_text = [title_text sprintf(' (%.2f ms relative to spike)', (frame - params.spike_frame) / params.sampling_rate * 1000)];
    else
        title_text = [title_text sprintf(' (%d relative to spike)', frame - params.spike_frame)];
    end
end

% finally add the title
title(title_text);

