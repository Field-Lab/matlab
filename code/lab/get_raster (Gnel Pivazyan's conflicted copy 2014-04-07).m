function spikes_by_trials = get_raster(spike_times, trial_begin_times, varargin)
%
% raster_gdf(spike_times,cell_specification,trigger,varargin)
%  
% Plots spike rasters
%
% Inputs: 
%   spike_times         List of spike times
%   trial_begin_times   List of times at which each new trial begins
%   
% Outputs:
%   spikes_by_trials      cell array storing the spike times of each trial
%
% optional parameters
%   start           0           if you want to start plot at time other than
%                           trigger_times, start will be trigger_times +
%                           start
%   stop            []          if user wants to specify when to end trial
%   first_trial     1           allows user to start with a later trial
%   tic_color       'k'         color of raster ticks
%   tic_thickness   1         thickness of each tic
%   line_space      0.2       vertical space between tics
%   foa             0         calls to set_up_figure_or_axes function
%

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParamValue('start', 0, @isnumeric);
p.addParamValue('stop', [], @isnumeric);
p.addParamValue('tic_color', [0 0 0]);
p.addParamValue('line_space', 0.2, @isnumeric)
p.addParamValue('tic_thickness', 1);
p.addParamValue('first_trial', 1, @isnumeric)
p.addParamValue('foa', 0)
p.addParamValue('plot', true);
p.addParamValue('labels', true);

% parse inputs
p.parse(varargin{:});
params = p.Results;


%%% function begins here %%%

% figure out the duration of each trial
if isempty(params.stop)
    params.stop = mean(diff(trial_begin_times));
end

% setup the figure
if params.plot
    plot_axes = set_up_fig_or_axes(params.foa);
    xlabel('seconds')
    ylabel('trials')
end

% calculate the number of trials to plot
num_trials = length(trial_begin_times)-params.first_trial+1;

% initialize returned variable and counter
spikes_by_trials = cell(num_trials, 1);
trial_plot_line = 1; % counter

for i = params.first_trial:num_trials
    % shift spikes
    ith_shifted_times = spike_times - trial_begin_times(i);  
    % get spikes in trial frame
    ith_trial_times = ith_shifted_times(ith_shifted_times >= params.start & ith_shifted_times <= params.stop);
    % plot the spike times
    if params.plot
        if ~isempty(ith_trial_times)
            line([ith_trial_times, ith_trial_times], [trial_plot_line-1+params.line_space, trial_plot_line-params.line_space],...
                    'color', params.tic_color, 'LineWidth', params.tic_thickness)
        end
    end
    % keep track of trial plotted
    trial_plot_line = trial_plot_line + 1; 
    % store spike times
    spikes_by_trials{i} = ith_trial_times;
end

