function arrayForCluster = generate2ElectrodePatternsInCluster(clusterElectrodes)

nSurroundingElectrodes = length(clusterElectrodes)-1;
arrayForCluster = zeros(length(clusterElectrodes), nSurroundingElectrodes*16+2);
for j = 1:4
    for k = 1:nSurroundingElectrodes
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4) = 0.25*j;
      arrayForCluster(1, k+(j-1)*nSurroundingElectrodes*4) = 1;
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes) = 0.25*j;
      arrayForCluster(1, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes) = -1;
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*2) = -0.25*j;
      arrayForCluster(1, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*2) = 1;
      arrayForCluster(1+k, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*3) = -0.25*j;
      arrayForCluster(1, k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*3) = -1;
    end
end
arrayForCluster(1, nSurroundingElectrodes*16+1) = 1;
arrayForCluster(1, nSurroundingElectrodes*16+2) = -1;
for i = 1:nSurroundingElectrodes
    arrayForCluster(1+i, nSurroundingElectrodes*16 + 2*(i-1) + 3) = 1;
    arrayForCluster(1+i, nSurroundingElectrodes*16 + 2*(i-1) + 4) = -1;
end