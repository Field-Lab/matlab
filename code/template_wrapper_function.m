function datarun = my_wrapper_function(datarun, cell_spec, arg1, varargin)
% MY_FUNCTION     This template function does nothing.  Ha!
%
% usage:  datarun = my_function(datarun, arg1, <params>)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_indices for options)
%                arg1 - argument 1
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field X
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot.  if -1, plot in current.
% foo               'bar'               how to activate the rings
%                                           'bar' - activate on site
%                                           'bore' - activate remotely
%
%
%
% parameters passed on to other functions. if not specified by the user, these parameters are not passed.
%
%
%       passed to 'compute_temporal_matrix'
%
%           time_bins       bins     	number of bins
%           frames          frames     	which sta frames to use 
%
%
%       passed to 'compute_spatial_matrix'
%
%           space_bins     	bins        number of bins
%           space_color 	frames      which color
%
%
%       passed to 'compute_temporal_matrix' AND 'compute_spatial_matrix'
%
%           bins                        number of bins for spatial AND temporal matrix.
%                                           this value is overridden by 'time_bins' or 'space_bins'
%
%
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

% parameters to be passed on
%    space summary function
p.addParamValue('space_bins', 'default value');
%    time summary function
p.addParamValue('time_bins', 'default value');
p.addParamValue('time_color', 'default value');
p.addParamValue('frames', 'default value');
%    both
p.addParamValue('bins', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
time_params = make_struct_to_pass(p.Results,{'color','time_color','bins',{'time_bins','bins'}});
space_params = make_struct_to_pass(p.Results,{'frames','frames','bins',{'space_bins','bins'}});





% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% show output
if params.verbose
    fprintf('\nComputing something for %d cells...',length(cell_indices));
    start_time = clock; % note when it started
end


% loop through cells
for cc = 1:length(cell_indices)
    
    % get cell index and id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % ...
end


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

