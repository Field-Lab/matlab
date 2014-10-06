function [time cf] = get_correlation(datarun, cell_IDs, varargin)

% MY_FUNCTION     This function gets the ACF (if one cell ID given) or CCF (if 2 cell IDs
%                      given) and plots the correlogram
%                 correlation functions are calculated using compute_ccf
%
% usage:  datarun = get_correlation(datarun, cell_IDs, <params>)
%
% arguments:  datarun - datarun struct with fields
%                               datarun.cell_ids - list of cell ids
%                               datarun.cell_types - cell types listed in standard format
%                               datarun.spikes
%            cell_IDs - single cell (for autocorrelogram) or pair of cells (for cross-correlogram)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    plot
%
% optional parameters, their default values, and what they specify:
%
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty,
%                                       don't plot.  if -1, plot in current.
%
%
% parameters passed on to other functions. if not specified by the user, these parameters are not
% passed.
%       passed to 'compute_ccf'
%
%                         [ dt       -  bin size (s) [offset/128]
%                         [ offset   -  maximum offset for CCF (s) [100e-3]
%                         [ shuffle  -  string to determine type of shuffle
%                                        'none'  - standard CCF   [default]
%                                        'stim' - signal correlations
%                                        'noise' - noise correlations
%                                           Must include 'trial' if latter 
%                                           two options are used.
%                         [ trial    -  trial duration of shuffle
%
%
%
% September 2010, Lauren HRUBY Jepson
%%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

%p.keepUnmatched = true;

p.addRequired('datarun', @isstruct)
p.addRequired('cell_IDs', @isnumeric)
    
% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('options', struct('offset',100e-3,'scale','ms','shuffle','none'), @isstruct) %passed to compute_ccf

% resolve user input and default values
p.parse(datarun, cell_IDs, varargin{:});

% get params struct
params = p.Results;
options = params.options;



%% BODY OF THE FUNCTION

% show output
if params.verbose
    disp(['Computing correlation for cell(s) ', num2str(cell_IDs)]);
    start_time = clock; % note when it started
end


% get cell indices
cell_indices = get_cell_indices(datarun, cell_IDs);

spikes{1} = datarun.spikes{cell_indices(1)};

if length(cell_indices) == 1
    %cell_indices = [cell_indices cell_indices];
    spikes{2} = spikes{1};
else
    spikes{2} = datarun.spikes{cell_indices(2)};
end

[cf, time] = compute_ccf_fix(spikes{1}, spikes{2}, options);
if length(cell_indices) == 1 || cell_IDs(1) == cell_IDs(2) %autocorrelation, so keep only for time > 0
    time = time(time>0);
    cf = cf(time>0);
end


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);
if ~isempty(plot_axes)
%     if length(cell_indices) == 1 || cell_IDs(1) == cell_IDs(2) %autocorrelation, so plot only for time > 0
%         plot(1000*time(time>0), cf(time>0))
%     else
        plot(1000*time, cf)
%     end
    
    xlabel('ms')
%     if length(cell_IDs) == 1 || cell_IDs(1) == cell_IDs(2)
%         title(['autocorrelogram for cell ' num2str(cell_IDs)])
%     else
%         title(['cross-correlogram for cells ' num2str(cell_IDs(1)) ', ' num2str(cell_IDs(2))])
%     end
end


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end



