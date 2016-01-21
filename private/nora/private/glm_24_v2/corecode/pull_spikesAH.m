% Function to pull spikes in a specified window along with negative spikes (happening just before window)
%AH update of Chaitu   10-4

%-- Arguments
% spikes_lng:  cell array of spike times, indexed by cell_ids
% cell_ids:    obvious   cells we are pulling from
% neighbors_id: cell ids of interest
% spikeT:      total number of bins in the time window, at spike train resolution
% duration:    duration of time window
% spike_dt:    timestep of spike train
% offset:      time (in seconds) offset from stimulus beginning

%-- Returns
% sp_times:    cell array of analog spike times of cells of interest, normalized to be in specified window
% D - logical: version of this
% negSpikes:   logical matrix of spikes occuring before window beginning

function [sp_times D negSpikes] = pull_spikesAH(spikes_lng,cell_ids,spikeT,duration,spike_dt,offset)

if (nargin < 7)
    offset = 0;
end

%if (isempty(cell_ids))
%    fprintf('pull_spikes: assuming spikes times are already ordered!\n');
%    cell_ids = 1:length(spikes_lng);
%end

if (nargout > 1)
    D = sparse(logical(false(spikeT,length(cell_ids))));
end

if (offset > 0 && nargout > 2)
    % NEG SPIKES HAS SOMETHING TO DO WITH SOME OFFSET
    negSpikes = sparse(logical(false(int32(offset/spike_dt),length(cell_ids)))); % Spikes occurring before time offset
else
    negSpikes = [];
end

sp_times = cell(length(cell_ids),1);

for k=1:length(cell_ids)
    spt = spikes_lng{k}; % pull spikes
    sp_times{k} = spt(spt > offset & spt < offset+duration) - offset; % subtract the time offset
    
    if (nargout > 1)
        sp_idx = ceil(sp_times{k}/spike_dt);
        D(sp_idx,k) = true;
    end
    
    if (offset > 0 && nargout > 2)
        negsp_idx = ceil(spt(spt > 0 & spt <= offset)/spike_dt);
        negSpikes(negsp_idx,k) = true;
    end
end
