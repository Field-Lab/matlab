function [waveformAvg nCells] = spontWaveform(spontSample, primaryIndex, samplerate, spontVectorLength, neighborsInSaveElec)

% Locates spikes within a sample of spontaneous activity and runs spike-sorting on it.
%
% inputs: 
%
% spontSample = array of spontaneous activity on a set of electrodes (primary and nearest neighbor)
% with dimensions (electrode) x (time sample)
% 
% primaryIndex = index of the primary electrode in spontSample
% spontVectorLength = length (in samples) of each inter-stimulus segment of spontaneous activity
% samplerate = sample rate, in samples / ms
% nSaveElec

% ***threshold-crossing events that begin within 2 ms of the end of a segment of recording or 1 ms
% of the beginning of a segment of recording are thrown out


nElec = length(squeeze(spontSample(:,1)));
spontVectorLength = length(squeeze(spontSample(1,:)))/totalPulses;

%% spike-finding on spontaneous sample

sigma = robust_std(spontSample(primaryIndex, :));
threshold = sigma*3;

threshCross = spontSample(primaryIndex,:) <= -threshold;
threshCrossSeg = bwlabel(threshCross);
nRawSpikes = max(threshCrossSeg);

rawSpikeStartIndex = zeros(nRawSpikes,1);
for i = 1:nRawSpikes;
    spikeSegment = threshCrossSeg == i;
    rawSpikeStartIndex(i) = find(spikeSegment,1 ,'first');
end

clear spikeSegment threshCrossSeg threshCross

%% removing spikes near "edges" of spontaneous vectors

badSpikes = zeros(nRawSpikes, 1);

%removing spikes near "edges" of spontaneous vectors
for i = 0:totalPulses
    tooLate = rawSpikeStartIndex >= i*spontVectorLength - 2*samplerate;
    tooEarly = rawSpikeStartIndex <= i*spontVectorLength + samplerate;
    badSpikes = tooLate.*tooEarly + badSpikes;
end

%remove spikes too close to other spikes (avoid counting a spike twice)
%(if two threshold-crossings are within 0.5 ms of eachother, the first is deleted)
for i = 1:nRawSpikes-1
    if abs(rawSpikeStartIndex(i) - rawSpikeStartIndex(i+1)) < samplerate/2
        badSpikes(i) = 1;
    end
end

%goodSpikes = array of ones and zeros representing which spike times in rawSpikeStartIndex correspond
%with spikes that meet above criteria
goodSpikes = 1 - badSpikes;


nCleanSpikes = sum(goodSpikes);
cleanSpikeStartIndexSparse = goodSpikes.*rawSpikeStartIndex;
cleanSpikeStartIndex = zeros(nCleanSpikes, 1);

clear rawSpikeStartIndex

k = 0;
for i = 1:nRawSpikes
    if cleanSpikeStartIndexSparse(i) ~= 0;
       k = k + 1;
       cleanSpikeStartIndex(k) = cleanSpikeStartIndexSparse(i);
    end
end

clear cleanSpikeStartIndexSparse tooEarly tooLate badSpikes goodSpikes nRawSpikes

%% extracting spontaneous spike waveforms
% takes waveform starting 10 samples before negative peak and extending 14 samples beyond negative
% peak (25 total)

cleanSpikeMinIndex = zeros(nCleanSpikes);
spikeVectors = zeros(nCleanSpikes, 25*nSpontPCAElec);

figure
hold on
for i = 1:nCleanSpikes
    spikeChunk = spontSample(primaryIndex, cleanSpikeStartIndex(i)-samplerate : cleanSpikeStartIndex(i)+samplerate);
    [ans,tempMin] = min(spikeChunk); %#ok<NOANS>
    cleanSpikeMinIndex(i) = tempMin + cleanSpikeStartIndex(i) - samplerate - 1;
    k = 0;
    for j = 1:nSaveElec
        if (neighborsInSaveElec(j)||j==primaryIndex)
            k = k + 1;
            spikeVectors(i,(k-1)*25+1:k*25) = spontSample(j, cleanSpikeMinIndex(i)-10 : cleanSpikeMinIndex(i)+14);
        end
    end
    plot(spikeVectors(i,:))
end
hold off

%% script to run jeff's PCA functions

clear dataset

dataset.raw_spikes = [(1:size(spikeVectors, 1))' spikeVectors];

% num samples
dataset.samples_before = 10;
dataset.samples_after = 14;
dataset.electrodes = 1:nSpontPCAElec;

% fill in other parameters
dataset.window_length = dataset.samples_before + dataset.samples_after + 1;
dataset.triggers = [];
dataset.mdf_file = [];

[dataset] = spike_prepare_spikes(dataset,'pca',1:nSpontPCAElec,99,1);
%99 refers to figure number, 1 refers to projection index

nCells = input('Number of cells to pull out:\n');
close

proj_struct = dataset.pca{1};
proj_struct_new = spike_select_several_clusters(proj_struct, 98, nCells); %  3 refers to number of clusters to select, 98 refers to figure number

%% finds average spike waveform for spontaneously firing cells

waveformAvgTemp = zeros(nCells, 25*nSpontPCAElec);

for i = 1:nCells
    waveformSum = zeros(1, 25*nSpontPCAElec);
    for j = 1:length(proj_struct_new.selected_indices{i});
        waveformSum = squeeze(spikeVectors(proj_struct_new.selected_indices{i}(j), :)) + waveformSum;
    end
    waveformAvgTemp(i,:) = waveformSum/length(proj_struct_new.selected_indices{i});
end

waveformAvg = zeros(nCells, nSpontPCAElec, 25);
for i = 1:nSpontPCAElec
    for j = 1:nCells
        waveformAvg(j,i,:) = waveformAvgTemp(j, (i-1)*25 + 1 : i*25);
    end
end
