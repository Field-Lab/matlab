%determines SNR for a raw data file

%%
clear all
close all

samplerate = 20000;

pathToData = '/Volumes/Lee/Data/Lauren/spike_size_analysis/2008-03-25-5/data001/data001000.bin';
primaryElec = 3;

sampleStart = 10;       %start of where sample should be taken, in seconds
sampleLength = 10;      %amount of data to be extracted, in seconds

% if sample is made up of smaller chunks of time, this is the length of each chunk (in samples), 
% otherwise set to samplerate*sampleLength
spontSegmentLength = samplerate*sampleLength;

%opens data file
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(pathToData);

%gets data from sample
data = rawFile.getData(samplerate*sampleStart, samplerate*(sampleStart + sampleLength));

% may not be necessary
addpath /Applications/MATLAB74/work/Lauren/analysis/jeffSpikeSorting

[xCoords yCoords] = getElectrodeCoords61();

electrodeCluster = [];
for i = 1:64
    if norm([xCoords(primaryElec) - xCoords(i), yCoords(primaryElec) - yCoords(i)]) < 2.1
        electrodeCluster = [electrodeCluster i]; %#ok<AGROW>
    end
end
primaryIndex = find(electrodeCluster == primaryElec);

%offset of one because of difference in indexing
spontSample = data(:, electrodeCluster + 1);

[waveformAvg, nCells] = spontWaveformReworked(spontSample, primaryIndex, samplerate/1000, spontSegmentLength);

nElec = length(electrodeCluster);

for i = 1:nCells
    maxValue = max(max(waveformAvg(i,:,:)));
    minValue = min(min(waveformAvg(i,:,:)));
    figure
    set(gcf, 'Position', [50 50 1000 200])
    for j = 1:nElec
        subplot(1, nElec, j)
        plot(squeeze(waveformAvg(i,j,:)))
        axis([0 25 minValue maxValue])
    end
end


