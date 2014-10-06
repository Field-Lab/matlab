function arrayForCluster = generatePatternClusterStim1Secondary(centerElec, electrodes, relAmps)

% generates array of values representing relative amplitudes of pulses on each electrode,
% with first dimension corresponding to electrodes (first electrode = primary) and second dimension
% corresponding to different stimulus patterns
%
% electrodes: vector of electrode numbers (all electrodes in cluster), with any electrode
% being the primary electrode
%
% relAmps: relative amplitudes of secondary electrodes (do not include 0!!!)
%
% the first patterns are combinations of the primary and secondary
% electrodes
%
% next 2 patterns are primary electrode alone (each polarity)
% next 2*nSurroundingElectrodes are secondary electrodes alone (each polarity)
nNeighbors = length(electrodes)-1;
nRelAmps = length(relAmps);

iCenter = find(electrodes == centerElec);
iNeighbors = find(electrodes ~= centerElec);


arrayForCluster = zeros(length(electrodes), nNeighbors*(nRelAmps+2)+2);
for j = 1:nRelAmps
    for k = 1:nNeighbors
      arrayForCluster(iNeighbors(k), k+(j-1)*nNeighbors) = relAmps(j);
      arrayForCluster(iCenter,   k+(j-1)*nNeighbors) = 1;
    end
end

% primary electrodes alone
arrayForCluster(iCenter, nNeighbors*nRelAmps+1) = 1;
arrayForCluster(iCenter, nNeighbors*nRelAmps+2) = -1;

% secondary electrodes alone
for i = 1:nNeighbors
    arrayForCluster(iNeighbors(i), nNeighbors*nRelAmps + 2*(i-1) + 3) = 1;
    arrayForCluster(iNeighbors(i), nNeighbors*nRelAmps + 2*(i-1) + 4) = -1;
end