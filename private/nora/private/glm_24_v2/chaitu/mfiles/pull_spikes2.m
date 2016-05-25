% Function to pull spikes in a specified window along with negative spikes (happening just before window)

%-- Arguments
% spikes_lng:  cell array of spike times, indexed by cell_ids
% cell_ids:    obvious
% cells_touse: cell ids of interest
% spikeT:      total number of bins in the time window, at spike train resolution
% duration:    duration of time window
% spike_dt:    timestep of spike train
% offset:      time (in seconds) offset from stimulus beginning

%-- Returns
% sp_times:    cell array of analog spike times of cells of interest, normalized to be in specified window
% D - logical: version of this
% negSpikes:   logical matrix of spikes occuring before window beginning

% version 2 -- modified by edoi, 2012-01-23.

function [sp_times D negSpikes] = pull_spikes2(spikes_lng,cell_ids,cells_touse,spikeT,duration,spike_dt,offset)

if (nargin < 7)
   offset = 0;
end

if (isempty(cell_ids))
   fprintf('pull_spikes: assuming spikes times are already ordered!\n');
   cell_ids = 1:length(spikes_lng);
end

if (nargout > 1)
   D = sparse(logical(false(spikeT,length(cells_touse))));
end

if (offset > 0 && nargout > 2)
   negSpikes = sparse(logical(false(int32(offset/spike_dt),length(cells_touse)))); % Spikes occurring before time offset
else
   negSpikes = [];
end

sp_times = cell(length(cells_touse),1);

if ~iscell(spikes_lng)
   tmp = spikes_lng; clear spikes_lng
   spikes_lng = cell(1);
   spikes_lng{1} = tmp;
end

for k=1:length(cells_touse)
   spt = spikes_lng{find(cell_ids ==  cells_touse(k),1)}; % pull spikes
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
