cd '/Data.noindex/Lauren/electrical_stimuli/2008-11-09/'; 

%pathToNeuronsFile = '/Volumes/Lee/Analysis/Lauren/2008-04-22-2/data001-map/data001-map.neurons';
pathToNeuronsFile = '/Volumes/Lee/Analysis/Lauren/2008-11-09/data001-map/data001-map.neurons';
cellIDs = [932 752 796 947 816 710];
cellElectrodes = [63 48 41 8 55 51]; %must correspond with cellIDs
stimAmps = [0.91 1.2 0.53 0.43 1.1 1.5]; % in uA, peak cathodic phase, corresponing with cellIDs and cellElectrodes

movingBarStartTime = 0; %in samples (beginning of one visual stimulus repetition)
movingBarEndTime = 60000; %last sample won't be included in retrieved spike times


nSamples = (movingBarEndTime - movingBarStartTime);
nCells = length(cellIDs);

electrodes=1:64;
Array=zeros(64);
for i = 1:length(stimAmps)
    Array(cellElectrodes(i),cellElectrodes(i)) = stimAmps(i);
end


spikeTimesInWindow = getNeuronSpikeTimes(pathToNeuronsFile, cellIDs, movingBarStartTime, movingBarEndTime);

totalSpikeTimes = [];
spikeTimeCellIds = [];
for i = 1:length(spikeTimesInWindow)
    totalSpikeTimes = [totalSpikeTimes; spikeTimesInWindow{i}];
    spikeTimeCellIds = [spikeTimeCellIds; cellIDs(i)*ones(length(spikeTimesInWindow{i}),1)];
end

totalSpikeTimesWithIds = [totalSpikeTimes, spikeTimeCellIds];
totalSpikeTimesWithIdsSorted = sortrows(totalSpikeTimesWithIds,1);


movieChunkArray = zeros(1, 6 + size(totalSpikeTimesWithIds,1)*3);
movieChunkArray(1,6) = nSamples; % 7th value in movieChunk gives length of movie in samples

for i = 1:size(totalSpikeTimesWithIds,1)
    j = 7+(i-1)*3; %every third value after the first 5
    movieChunkArray(1,j)   = totalSpikeTimesWithIdsSorted(i,1); %time at which to play pattern
    movieChunkArray(1,j+1) = cellElectrodes(cellIDs==totalSpikeTimesWithIdsSorted(i,2)); %pattern number
    movieChunkArray(1,j+2) = 1; %scaling factor: for now always = 1
end

movieChunkArray = [length(movieChunkArray) movieChunkArray]; %first value in movie chunks denotes length of vector
MovieChunksFile = [1 movieChunkArray]; %first value of movie chunk file denotes number of chunks in file

fid = fopen('moving_bar_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('moving_bar_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('moving_bar_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 