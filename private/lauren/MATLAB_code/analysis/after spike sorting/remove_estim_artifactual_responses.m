function [spike_trains_out spikes_removed removed_bin] = remove_estim_artifactual_responses(spike_trains, stim_times, blank_window, varargin)
%REMOVE_ESTIM_ARTIFACTUAL_RESPONSES
%   removes artifactual responses to electrical stimulation (e.g. due to
%   merged axons) by removing any spikes within a specified time window
%   relative to the stimulus pulses
%
%   spike_train: cell array of vectors of spike times, in which each cell
%      element corresponds to one spike train (trial)
%   stim_times: vector of times when stimulus is applied
%   blank_window: window of time, relative to stimulus application, where
%      artifactual responses are expected (2x1 vector)
%   *** time units may be ms or samples, but must be consitent across
%   arguments!



p = inputParser;

p.addRequired('spike_trains', @iscell)
p.addRequired('stim_times', @isnumeric)
p.addRequired('blank_window', @isnumeric)


p.addParamValue('verbose', false, @islogical)
p.addParamValue('plot_results', false, @islogical)
p.addParamValue('figTitle', [], @ischar)
p.parse(spike_trains, stim_times, blank_window, varargin{:})

params = p.Results;

nTrials = length(spike_trains);

%keep track of which and how many spikes are removed
nRemoved = 0;
spikes_removed =   cell(size(spike_trains));
removed_bin =      cell(size(spike_trains));
spike_trains_out = cell(size(spike_trains));


for jj = 1:nTrials
    toRemove = false(size(spike_trains{jj}));
    for ii = 1:length(stim_times)
        blank_reg = stim_times(ii)+blank_window;
        toRemove(spike_trains{jj} >= blank_reg(1) & spike_trains{jj} <= blank_reg(2)) = true;
    end
    spikes_removed{jj} = spike_trains{jj}(toRemove);
    spike_trains_out{jj} = spike_trains{jj}(~toRemove);
    removed_bin{jj} = toRemove;
    nRemoved = nRemoved + sum(toRemove);
end



if params.verbose
    disp([num2str(nRemoved) ' artifactual spikes were removed from the collection of spike trains'])
end

if params.plot_results
    figure
    hold on
    for kk = 1:nTrials
        for jj = 1:length(spike_trains_out{kk})
            plot([spike_trains_out{kk}(jj) spike_trains_out{kk}(jj)], [kk-1 kk], 'k-', 'LineWidth', 1)
        end
        for jj = 1:length(spikes_removed{kk})
            plot([spikes_removed{kk}(jj) spikes_removed{kk}(jj)], [kk-1 kk], 'r-', 'LineWidth', 1)
        end
    end
    set(gca, 'yLim', [0 nTrials])
end

