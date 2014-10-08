
cd '/Volumes/War/Analysis/Lauren/2008-11-10-3';
pathToNeuronsFile = '/Volumes/War/Analysis/Lauren/2008-11-10-3/data002-map/data002-map.neurons';
%cellIDs = [466 543 542 571 886 721 708 692 646 31 226];


% cells{1} = [407 422 466 467 481 497];
% cells{2} = [318 349 377 378 409 603 916];
% cells{3} = [228 261 274 289 335];
% cells{4} = [6 38 78 98 157 170 172 202];
% cells{5} = [574 589 618 619 634 663 664];
% cells{6} = [814 873 903 919];
% cells{7} = [754 767];
% cells{8} = [917 932 948];
% cells{9} = [635 648 678 695 723 737];
% 
% 
% cellElectrodes = [29 28 19 6 43 59 52 61 46]; %must correspond with cellIDs
% stimAmps = [1 1 2.5 1.81 1.20 1 1.30 1.81 2.0]; % in uA, corresponing with cellIDs and cellElectrodes

cellElectrodes = [19 6 52];
stimAmps = [2.4 2.2 1.3];

movingBarStartTime = 60000; %in samples (beginning of one visual stimulus repetition)
movingBarEndTime = 120000; %last sample won't be included in retrieved spike times


nSamples = (movingBarEndTime - movingBarStartTime);
nCells = length(cells);

electrodes=1:64;
Array=zeros(64);
for i = 1:length(stimAmps)
    Array(cellElectrodes(i),cellElectrodes(i)) = stimAmps(i);
end

% spikeTimesInWindow = catAndCleanDuplicateCellSpikeTimes(pathToNeuronsFile, cells, [movingBarStartTime movingBarEndTime]);
% %spikeTimesInWindow = getNeuronSpikeTimes(pathToNeuronsFile, cellIDs, movingBarStartTime, movingBarEndTime);
% 
% totalSpikeTimes = [];
% spikeTimeCellIds = [];
% for i = 1:length(spikeTimesInWindow)
%     totalSpikeTimes = [totalSpikeTimes; spikeTimesInWindow{i}];
%     spikeTimeCellIds = [spikeTimeCellIds; i*ones(length(spikeTimesInWindow{i}),1)];
% end
% 
% totalSpikeTimesWithIds = [totalSpikeTimes, spikeTimeCellIds];
% totalSpikeTimesWithIdsSorted = sortrows(totalSpikeTimesWithIds,1);
% 
% 
% movieChunkArray = zeros(1, 6 + size(totalSpikeTimesWithIds,1)*3);
% movieChunkArray(1,6) = nSamples; % 7th value in movieChunk gives length of movie in samples
% 
% for i = 1:size(totalSpikeTimesWithIds,1)
%     j = 7+(i-1)*3; %every third value after the first 5
%     movieChunkArray(1,j)   = totalSpikeTimesWithIdsSorted(i,1); %time at which to play pattern
%     movieChunkArray(1,j+1) = cellElectrodes(totalSpikeTimesWithIdsSorted(i,2)); %pattern number
%     movieChunkArray(1,j+2) = 1; %scaling factor: for now always = 1
% end
% 
% movieChunkArray = [length(movieChunkArray) movieChunkArray]; %first value in movie chunks denotes length of vector
% MovieChunksFile = [1 movieChunkArray]; %first value of movie chunk file denotes number of chunks in file
% 
% plotStimulusMovie(MovieChunksFile, Array, electrodes);

%cd '/Data.noindex/Lauren/electrical_stimuli'; 
fid = fopen('moving_bar_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('moving_bar_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('moving_bar_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 