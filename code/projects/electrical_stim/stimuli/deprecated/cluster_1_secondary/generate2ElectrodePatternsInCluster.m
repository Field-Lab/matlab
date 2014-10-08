function arrayForCluster = generate2ElectrodePatternsInCluster(clusterElectrodes)

% generates array of values representing relative amplitudes of pulses on each electrode,
% with first dimension corresponding to electrodes (first electrode = primary) and second dimension
% corresponding to different dimensions
%
% clusterElectrodes: vector of electrode numbers (all electrodes in cluster), with first electrode
% being the primary electrode (should be changed to nClusterElectrodes, because actual electrode
% identities are not used!!!)
%
% the first 16*nSurroundingElectrodes patterns are combinations of the primary and secondary
% electrodes, in the following order:
%
% primary = 1, secondary #1 = 0.25
% primary = 1, secondary #2 = 0.25
% primary = 1, secondary #3 = 0.25
% ... (rest of secondaries)
% primary = -1, secondary #1 = 0.25
% primary = -1, secondary #2 = 0.25
% primary = -1, secondary #3 = 0.25
% ... (rest of secondaries)
% primary = 1, secondary #1 = -0.25
% ... (rest of secondaries)
% primary = -1, secondary #1 = -0.25
% ... (rest of secondaries)
% primary = 1, secondary #1 = 0.5
% ... (rest of secondaries)
% primary = -1, secondary #1 = 0.5
% ... (rest of secondaries, rest of polarity combinations)
% primary = 1, secondary #1 = 0.75
% ... (rest of secondaries, rest of polarity combinations)
% primary = 1, secondary #1 = 1

% next 2 patterns are primary electrode alone (each polarity)
% next 2*nSurroundingElectrodes are secondary electrodes alone (each polarity)


nSurroundingElectrodes = length(clusterElectrodes)-1;
arrayForCluster = zeros(length(clusterElectrodes), nSurroundingElectrodes*16+2);
for j = 1:4 %4 different secondary amplitudes
    for k = 1:nSurroundingElectrodes
      %both +
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4) = 0.25*j;
      arrayForCluster(1,   k+(j-1)*nSurroundingElectrodes*4) = 1;
      %primary -, secondary +
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes) = 0.25*j;
      arrayForCluster(1,   k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes) = -1;
      %primary +, secondary -
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*2) = -0.25*j;
      arrayForCluster(1,   k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*2) = 1;
      %primary -, secondary -
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*3) = -0.25*j;
      arrayForCluster(1,   k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*3) = -1;
    end
end

% primary electrodes alone
arrayForCluster(1, nSurroundingElectrodes*16+1) = 1;
arrayForCluster(1, nSurroundingElectrodes*16+2) = -1;

% secondary electrodes alone
for i = 1:nSurroundingElectrodes
    arrayForCluster(1+i, nSurroundingElectrodes*16 + 2*(i-1) + 3) = 1;
    arrayForCluster(1+i, nSurroundingElectrodes*16 + 2*(i-1) + 4) = -1;
end