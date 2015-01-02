function shifted_trial_raster(start, stop, spikes, cell_idxs, trigger)
%SHIFTED_TRIAL_RASTER plots a raster plot of spikes started at different times
%   SHIFTED_TRIAL_RASTER(start, stop, spikes, cell_idxs)
%       start -> When to start displaying spikes where 0 is the trigger or
%                start of the trial
%       stop -> When to cease displaying spikes where 0 is the trigger or
%               start of the trial
%       spikes -> is a [num_neurons x 1] cell containing [num_spikes x 1] 
%                 vectors corresponding to all the spike times of that
%                 specific cell.  Contains all the neurons picked up by the
%                 array.
%       cell_idxs -> contains the indices (for spikes) to identify the
%                    neurons during each trial run.
%       trigger -> What the trigger is for each trial.  If trigger is a
%                  single value, SHIFTED_TRIAL_RASTER will use the same
%                  trigger for each neuron.  If trigger is an array, it
%                  must have the same number of elements as cell_idxs and
%                  SHIFTED_TRIAL_RASTER will use the corresponding trigger
%                  for each neuron's raster plot.
%
% Marvin Thielk 2013
% mthielk@salk.edu

if length(trigger) == 1
    trigger = repmat(trigger, size(cell_idxs));
end

psth_r = [];

for i = 1:length(cell_idxs)
    spk_times = spikes{cell_idxs(i)} - trigger(i);
    spk_idxs = find(spk_times >= start & spk_times <= stop);
    psth_r = [psth_r; 1000*spk_times(spk_idxs), repmat(length(cell_idxs)-i,[length(spk_idxs),1]);];
end

if ~isempty(psth_r)
    plot(psth_r(:,1),psth_r(:,2),'k.');%, 'MarkerSize',10
    axis([start*1000 stop*1000 0 length(cell_idxs)]);
end

end