function datarun = getBinnedSpikes(datarun,start_stim,end_stim,stims,refresh_time,duration)

% get list of cells
cell_indices = get_cell_indices(datarun,'all');
 
% initialize storage variables
spike_rate = zeros(length(cell_indices),length(stims));
spike_times = cell(length(cell_indices),1);


% go through each cell
for cc = 1:length(cell_indices)
 
    cell_index = cell_indices(cc);
 
    % bin up spikes for entire duration
    spike_rate_ = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:datarun.duration);

    % take spikes just from the relevant subset and time-shift to align peak frame with stimulus
    %spike_rate_ = circshift(spike_rate_(start_stim:end_stim),time_offset);
    spike_rate_ = spike_rate_(start_stim:end_stim);
 
    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate_)
        spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
    end
    
    % put into storage variables
    spike_rate(cc,:) = spike_rate_;
    spike_times{cc} = sort(spike_times_);
    
end

datarun.spike_rate = spike_rate;