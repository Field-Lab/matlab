function result = my_function(datarun, arg1, params)
% MY_FUNCTION     This template function does nothing.  Ha!
%
% usage:  datarun = my_function(datarun, arg1, params)
%
% arguments:  datarun - datarun struct with field specifying X
%                arg1 - argument 1
%              params - struct of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional fields in params, their default values, and what they specify:
%
% param_1           'default_value_1'	does something
% param_2           false               does something else
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
%
%
% author date
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.param_1 = 'default_value_1';
defaults.param_2 = false;
defaults.verbose = false;
defaults.figure = [];

% combine user and default parameters
params = default_params( defaults, params);



if params.verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end


% BODY OF THE FUNCTION

% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);



% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

