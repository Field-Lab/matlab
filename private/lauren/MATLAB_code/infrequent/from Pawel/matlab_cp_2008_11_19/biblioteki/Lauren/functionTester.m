% code for generating arguments necessary to test various matlab functions



%% manualPCACluster.m tester %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

% loads some data to test with (make sure path is correct!)
load '/netapp/snle/home/lhruby/Desktop/testData.mat'

% calls function
clusterIndex = manualPCACluster(PCAData);

% plots representation of clusterIndex
clusterColor = hsv(max(clusterIndex));
nTraces = size(PCAData, 1);
nElectrodes = size(PCAData, 2);
nSamples = size(PCAData, 3);

prinCompArray = zeros(nTraces, nElectrodes*nSamples);
for i = 1:nElectrodes % concatenates traces on different electrodes electrode for each pulse
    prinCompArray(:, (i-1)*nSamples + 1 : i*nSamples) = PCAData(:,i,:);
end

[PCACoef, PCAScore] = princomp(prinCompArray);


figure
hold on

for i = 1:max(clusterIndex) % plots traces as different colors representing the cluster they belong to.
    for j = 1:length(clusterIndex)
        if clusterIndex(j) == i
            plot(PCAScore(j,1), PCAScore(j,2), '.', 'MarkerFaceColor', clusterColor(i,:), 'MarkerEdgeColor', clusterColor(i,:), 'MarkerSize', 20)
        end
    end
end

hold off








%% extractDataSubset3D tester %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

data = zeros(3, 10, 5);
for i = 1:3
    for j = 1:10
        for k = 1:5
            data(i,j,k) = i + j*10 + k*100;
        end
    end
end

d1Extract = [0 1 1];
d2Extract = [1 0 0 1 1 0 0 0 0 0];
d3Extract = [0 1 1 1 0];

extractedData = extractDataSubset3D(data, d1Extract, d2Extract, d3Extract);






%% PCChooser.m tester %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

% loads some data to test with (make sure path is correct!)
load '/netapp/snle/home/lhruby/Desktop/testData.mat'

nTraces = size(PCAData, 1);
nElectrodes = size(PCAData, 2);
nSamples = size(PCAData, 3);

% concatenates traces on different electrodes for each pulse
prinCompArray = zeros(nTraces, nElectrodes*nSamples);
for i = 1:nElectrodes
    prinCompArray(:, (i-1)*nSamples + 1 : i*nSamples) = PCAData(:,i,:);
end

% runs PCA
[PCACoef, PCAScore] = princomp(prinCompArray);

% calls function
[PCx PCy] = PCChooser(PCAScore);

