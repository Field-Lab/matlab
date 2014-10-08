% a script to test a few filters for removing/minimizing artifacts

clear all
close all

ArtifactSubtraction = 0;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArrayID = 1;

DataPath= '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007/';
ArtifactPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data015/';
PatternNumber = 1;
MovieNumber = 24;

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
centerChannel = 13;
Radius = 1;
ChannelsToUse=electrodeMap.getAdjacentsTo(centerChannel,Radius)';

%returns data traces: nTraces x nElectrodes x nSamples
[dataTraces, ArtifactDataTraces, Channels]=NS_ReadPreprocessedData(DataPath, ArtifactPath, ArtifactSubtraction, PatternNumber,MovieNumber);
artifactTraces = NS_ReadPreprocessedData(ArtifactPath, ArtifactPath, ArtifactSubtraction, PatternNumber,MovieNumber);

nTraces  = size(dataTraces, 1);
nSamples = size(dataTraces, 3);

%% gets "template" ei from .ei file

pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 227;

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples

%% difference function
deriv1 = squeeze(diff(dataTraces(:,centerChannel,:),1,3));
deriv2 = diff(deriv1,1,2);

deriv1Artifact = squeeze(diff(artifactTraces(:,centerChannel,:),1,3));
deriv2Artifact = squeeze(diff(deriv1Artifact,1,2));

figure
subplot(3,2,1)
plot(1:nSamples,squeeze(dataTraces(:,13,:)), 'k-')
subplot(3,2,3)
plot(1:nSamples-1, deriv1, 'm-')
subplot(3,2,5)
plot(1:nSamples-2, deriv2, 'b-')

subplot(3,2,2)
plot(1:nSamples,squeeze(artifactTraces(:,13,:)), 'k-')
subplot(3,2,4)
plot(1:nSamples-1, deriv1Artifact, 'm-')
subplot(3,2,6)
plot(1:nSamples-2, deriv2Artifact, 'b-')


%% difference function with smoothing of first derivative
movAvg = zeros(nTraces, nSamples+1);
for i = 1:nTraces
    movAvg(i,:) = conv(ones(1,3),squeeze(deriv1(i,:)));
end
deriv2smooth = diff(movAvg,1,2);

figure
title('data')
subplot(3,1,1)
plot(1:nSamples,squeeze(dataTraces(:,13,:)),'k-')
subplot(3,1,2)
plot(1:nSamples+1, movAvg)
subplot(3,1,3)
plot(1:nSamples, deriv2smooth)

%% difference function of template
eiTrace = squeeze(ei(1,centerChannel+1,:));
eiDeriv1 = diff(eiTrace);
eiDeriv2 = diff(eiDeriv1);

figure
title('ei signal')
subplot(3,1,1)
plot(eiTrace)
subplot(3,1,2)
plot(eiDeriv1)
subplot(3,1,3)
plot(eiDeriv2)


%deriv1template = squeeze(diff


%% matched filter

% time-reversing templates

revEiTrace = flipdim(eiTrace,1);
revEiDeriv1 = flipdim(eiDeriv1,1);
revEiDeriv2 = flipdim(eiDeriv2,1);

% filtering

filteredData = zeros(nTraces, length(revEiTrace) + nSamples - 1);
filteredDeriv1 = zeros(nTraces, length(revEiDeriv1) + size(deriv1, 2) - 1);
filteredDeriv2 = zeros(nTraces, length(revEiDeriv2) + size(deriv2, 2) - 1);

for i = 1:nTraces
    filteredData(i,:) = conv(revEiTrace, squeeze(dataTraces(i,centerChannel,:)));
    filteredDeriv1(i,:) = conv(revEiDeriv1, deriv1(i,:));
    filteredDeriv2(i,:) = conv(revEiDeriv2, deriv2(i,:));
end

figure
title('after filtering with ei signal (matched filter)')
subplot(3,1,1)
plot(1:size(filteredData, 2), filteredData)
subplot(3,1,2)
plot(1:size(filteredDeriv1, 2), filteredDeriv1)
subplot(3,1,3)
plot(1:size(filteredDeriv2, 2), filteredDeriv2)

%% fabricated error function

% truncating ei traces
eiShort = eiTrace(13:30);
eiDeriv1Short = eiDeriv1(13:29);
eiDeriv2Short = eiDeriv2(13:28);

error = zeros(nTraces, nSamples - length(eiShort));
errorDeriv1 = zeros(nTraces, nSamples - length(eiShort));
errorDeriv2 = zeros(nTraces, nSamples - length(eiShort));
for i = 1:nTraces
    for j = 1:nSamples - length(eiShort)
        error(i,j) = norm(eiShort - squeeze(dataTraces(i,centerChannel,j:j+length(eiShort)-1)));
        errorDeriv1(i,j) = norm(eiDeriv1Short - deriv1(i,j:j+length(eiDeriv1Short)-1)');
        errorDeriv2(i,j) = norm(eiDeriv2Short - deriv2(i,j:j+length(eiDeriv2Short)-1)');
    end
end

figure
subplot(3,1,1)
plot(1:size(error,2),error)
subplot(3,1,2)
plot(1:size(errorDeriv1,2),errorDeriv1)
subplot(3,1,3)
plot(1:size(errorDeriv2,2),errorDeriv2)






