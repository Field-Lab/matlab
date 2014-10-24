function [datarun, bins] = get_binned_spikes(datarun, stims, varargin)
% Doesn't yet do anything smart with triggers; just uses triggers(1).

opts = inputParser();
opts.addParamValue('triggers', datarun.triggers);
opts.addParamValue('cellspec', 'all');
opts.parse(varargin{:});
opts = opts.Results;

cell_indices = get_cell_indices(datarun, opts.cellspec);
 
% go through each cell
spike_rate = zeros(length(cell_indices), length(stims));
spike_times = cell(length(cell_indices), 1);
refresh_time = datarun.stimulus.refresh_time;
bins = opts.triggers(1) + (stims(1)-1:stims(end)).*refresh_time;
for cc = 1:length(cell_indices)
 
    cell_index = cell_indices(cc);
 
    % bin up spikes for entire duration
    spike_rate_ = histc(datarun.spikes{cell_index}, bins);
 
    % take spikes just from the relevant subset and time-shift to align peak frame with stimulus
    spike_rate_ = spike_rate_(stims);
 
    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate_)
        spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
    end
    
    % put into storage variables
    spike_rate(cc,:) = spike_rate_;
    spike_times{cc} = sort(spike_times_);
end

datarun.spike_rate(cell_indices,:) = spike_rate;