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

% plots traces as different colors representing the cluster they belong to.
for j = 1:length(clusterIndex)
    plot(PCAScore(j,1), PCAScore(j,2), '.', 'MarkerFaceColor', clusterColor(clusterIndex(j),:), 'MarkerEdgeColor', clusterColor(clusterIndex(j),:), 'MarkerSize', 20)
end

hold off


%% autoClusterImm.m tester %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

% loads some data to test with (make sure path is correct!)
load '/netapp/snle/home/lhruby/Desktop/testData.mat'

nTraces = size(PCAData, 1);
nElectrodes = size(PCAData, 2);
nSamples = size(PCAData, 3);

prinCompArray = zeros(nTraces, nElectrodes*nSamples);
for i = 1:nElectrodes % concatenates traces on different electrodes electrode for each pulse
    prinCompArray(:, (i-1)*nSamples + 1 : i*nSamples) = PCAData(:,i,:);
end

[PCACoef, PCAScore] = princomp(prinCompArray);

clusterIndex = autoClusterImm(PCAScore);

figure
hold on

% plots traces as different colors representing the cluster they belong to.
clusterColor = hsv(max(clusterIndex));
for j = 1:length(clusterIndex)
    plot(PCAScore(j,1), PCAScore(j,2), '.', 'MarkerFaceColor', clusterColor(clusterIndex(j),:), 'MarkerEdgeColor', clusterColor(clusterIndex(j),:), 'MarkerSize', 20)
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

%% plotStimulusTraces.m tester %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-1 1],'FontSize',20,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1);

ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

%cd /Data.noindex/Lauren/2008-08-26-0/
cd D:\2008-08-26-0;
filename_movie='movie008';
l=length(filename_movie);
filename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)];

patternNumberToPlot = 5;

number_of_PD_chunk = 1;

%channels = [1];

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
channels=electrodeMap.getAdjacentsTo(3,1)';

figHandle = plotStimulusTraces(filename, number_of_PD_chunk, patternNumberToPlot, channels, ArrayID, FigureProperties, NS_GlobalConstants);

