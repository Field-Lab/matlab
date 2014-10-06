function [rho, bin_centers] = binspikes(spikes, dt, duration)
% BINSPIKES:  Bin a spike train consisting of spike times
%
% Usage: [rho, time] = binspikes(spikes, dt)
%        [rho, time] = binspikes(spikes, dt, duration)
%
% Arguments:
%           spikes -- two forms:
%                        [ 1xN vector of spike times
%                        [
%                        [ Mx1 cell array where the i'th array has
%                        [   N(i) spike times.
%               dt -- bin size (same units as spikes)
%         duration -- duration of spike train (default is max spike time)
%
% Outputs:
%              rho -- two forms:
%                        MxN vector of binned spike times (M=1 if vector)
%                              Note: outputs sparse vector
%
%             time -- vector of time corresponding to each bin center
% 

% find duration
if nargin ~= 3
  if iscell(spikes)
    allspikes = cat(1,spikes{:});
    duration = max(allspikes);
  else
    duration = max(spikes);
  end
end


% calculate centers of each time resolution bin
bin_centers = [0 : dt : duration] + dt/2;

% allocate memory for binned spike times
if iscell(spikes)
  rho = sparse(length(spikes),length(bin_centers));
else
  rho = sparse(1,length(bin_centers));
end

% This step is slightly tricky. Each bin_center denotes an area
% of time (the time resolution) at which we are examining the 
% spike times. If a spike time falls with in that time bin, then
% add one to our count for that time resolution, and so on..

if ~iscell(spikes)
  rho = sparse(hist(spikes, bin_centers));
else
  for i=1:length(spikes)
    rho(i,:) = sparse(hist(spikes{i}, bin_centers));
  end
end

% display the cleanliness of the spike
cleanliness = (prod(size(rho)) - length(find(rho>1))) / prod(size(rho)) * 100;    
if (cleanliness < 99)
  disp(['---> ' num2str(cleanliness) '% of bins are clean']);
end

% clean the bins
rho = sign(rho);

% remove the trash on the last bin
rho = rho(:,1:end-1);
bin_centers = bin_centers(1:end-1);