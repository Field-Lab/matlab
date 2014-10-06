function [waveformAvg nCells] = spontWaveformReworked(spontSample, primaryIndex, samplerateMs, spontSegmentLength)

% Locates spikes within a sample of spontaneous activity and runs spike-sorting on it.
%
% inputs: 
%
% spontSample = array of spontaneous activity on a set of electrodes (primary and nearest neighbor)
% with dimensions (sample) x (electrode)
% 
% primaryIndex = index of the primary electrode in spontSample dimension 2
% spontSegmentLength = length (in samples) of each segment of spontaneous activity
% samplerateMs = sample rate, in samples / ms

% ***threshold-crossing events that begin within 2 ms of the end of a segment of recording or 1 ms
% of the beginning of a segment of recording are thrown out


nElec = size(spontSample, 2);
nSegments = size(spontSample, 1)/spontSegmentLength;

if rem(size(spontSample, 1), spontSegmentLength) ~= 0
    error('the first dimension of spontSample isn''t a multiple of spontSegmentLength.  Aborting.')
end

%% spike-finding on spontaneous sample

sigma = robust_std(spontSample(:, primaryIndex));
threshold = sigma*3;

%corrects for any mean offset
meanValue = mean(spontSample, 1);
for i = 1:nElec
    spontSample(:,i) = spontSample(:,i) - meanValue(i);
end

threshCross = spontSample(:, primaryIndex) <= -threshold;
threshCrossSeg = bwlabel(threshCross);
nRawSpikes = max(threshCrossSeg);

rawSpikeStartIndex = zeros(nRawSpikes, 1);
for i = 1:nRawSpikes;
    spikeSegment = threshCrossSeg == i;
    rawSpikeStartIndex(i) = find(spikeSegment, 1,'first');
end

clear spikeSegment threshCrossSeg threshCross

%% removing spikes near "edges" of spontaneous Segments

badSpikes = zeros(nRawSpikes, 1);

% removing spikes near "edges" of spontaneous Segments (beginning of threshold-crossing within 
% 2 ms of end of segment or within 1 ms of beginning of segment)

for i = 0:nSegments
    tooLate = rawSpikeStartIndex >= i*spontSegmentLength - 2*samplerateMs;
    tooEarly = rawSpikeStartIndex <= i*spontSegmentLength + samplerateMs;
    badSpikes = tooLate.*tooEarly + badSpikes;
end

%remove spikes too close to other spikes (avoid counting a spike twice)
%(if two threshold-crossings are within 0.5 ms of eachother, the first is deleted)
for i = 1:nRawSpikes-1
    if abs(rawSpikeStartIndex(i) - rawSpikeStartIndex(i+1)) < samplerateMs/2
        badSpikes(i) = 1;
    end
end

% goodSpikes = array of ones and zeros representing which spike times in rawSpikeStartIndex correspond
% with spikes that meet above criteria
%goodSpikes = ~badSpikes;

cleanSpikeStartIndex = rawSpikeStartIndex(~badSpikes);

clear tooEarly tooLate badSpikes nRawSpikes rawSpikeStartIndex

%% extracting spontaneous spike waveforms
% takes waveform starting 10 samples before negative peak and extending 14 samples beyond negative
% peak (25 total)

nCleanSpikes = length(cleanSpikeStartIndex);
cleanSpikeMinIndex = zeros(nCleanSpikes);
spikeVectors = zeros(nCleanSpikes, 25*nElec);

figure
hold on
for i = 1:nCleanSpikes
    spikeChunk = spontSample(cleanSpikeStartIndex(i)-samplerateMs : cleanSpikeStartIndex(i)+samplerateMs, primaryIndex);
    [ans,tempMin] = min(spikeChunk); %#ok<NOANS>
    cleanSpikeMinIndex(i) = tempMin + cleanSpikeStartIndex(i) - samplerateMs - 1;
    for j = 1:nElec
        spikeVectors(i,(j-1)*25+1:j*25) = spontSample(cleanSpikeMinIndex(i)-10 : cleanSpikeMinIndex(i)+14, j);
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
dataset.electrodes = 1:nElec;

% fill in other parameters
dataset.window_length = dataset.samples_before + dataset.samples_after + 1;
dataset.triggers = [];
dataset.mdf_file = [];

[dataset] = spike_prepare_spikes(dataset,'pca',1:nElec,99,1);
%99 refers to figure number, 1 refers to projection index

nCells = input('Number of cells to pull out:\n');
close

proj_struct = dataset.pca{1};
proj_struct_new = spike_select_several_clusters(proj_struct, 98, nCells); %  3 refers to number of clusters to select, 98 refers to figure number

%% finds average spike waveform for spontaneously firing cells

waveformAvgTemp = zeros(nCells, 25*nElec);

for i = 1:nCells
    waveformSum = zeros(1, 25*nElec);
    for j = 1:length(proj_struct_new.selected_indices{i});
        waveformSum = squeeze(spikeVectors(proj_struct_new.selected_indices{i}(j), :)) + waveformSum;
    end
    waveformAvgTemp(i,:) = waveformSum/length(proj_struct_new.selected_indices{i});
end

waveformAvg = zeros(nCells, nElec, 25);
for i = 1:nElec
    for j = 1:nCells
        waveformAvg(j,i,:) = waveformAvgTemp(j, (i-1)*25 + 1 : i*25);
    end
end





