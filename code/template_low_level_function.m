function result = my_function(arg1, varargin)
% MY_FUNCTION     This template function does nothing.  Ha!
%
% usage:  result = my_function(arg1, <params>)
%
% arguments:     arg1 - argument 1
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure.
%                                           if empty, don't plot.  if -1, plot in current.
% foo               'bar'               how to activate the rings
%                                           'bar' - activate on site
%                                           'bore' - activate remotely
%
%
% date author
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('foo','bar', @(x)any(strcmpi(x,{'bar','bore'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% show output
if params.verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end


% ...

if ~isempty(plot_axes)
    axes(plot_axes)
    % plot the thing
end


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

