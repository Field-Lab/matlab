
%pathToData = '/Volumes/Palace/Data/Lauren/2000-12-14-1/data053/data053.bin';
pathToNeuronsFile = '/Volumes/Brokedown/Analysis/Lauren/2010-10-18-3/data002-from-data000/data002-from-data000.neurons';

nSamples = 20000; %make sure this fits with length of moving bar

cellIDs =        [36   48   247  393  531  694];
cellElectrodes = [2    4    17   21   36   47];
stimAmps =       [1.1  0.68 .35  1.1  .55  0.9];

%cellIDs = [154 334 466 588 767 813];
%cellElectrodes = [29 14 41 7 53 39]; %must correspond with cellIDs
%stimAmps = [0.2 0.3 0.4 0.2 0.3 0.4]; %in uA, corresponing with cellIDs and cellElectrodes


%% get triggers to determine start/end times of moving bar repitions

datarun.names.rrs_neurons_path = pathToNeuronsFile;
datarun = load_neurons(datarun);

triggers = datarun.triggers;
triggerInts = diff(triggers);
chosenInts = [];

triggersStartInd = 1:2:199;
triggersEndInd = 2:2:200;
triggersStart = triggers(triggersStartInd);
triggersEnd = triggers(triggersEndInd);


% prevTrigInt = false;
% prevPrevTrigInt = false;
% for i = 1:length(triggerInts)
%     if abs(triggerInts(i) - 0.4) < 0.02 %close to 0.4 seconds
%         if ~(prevTrigInt || prevPrevTrigInt) %first moving bar in set of 12
%             chosenInts = [chosenInts i]; %#ok<AGROW>
%         end
%         prevPrevTrigInt = prevTrigInt;
%         prevTrigInt = true;
%     else
%         prevPrevTrigInt = prevTrigInt;
%         prevTrigInt = false;
%     end
% end

%chosenInts = chosenInts(1:2:end); %take only white bars (every other)

%movingBarStartTimes = [0 40000 80000 120000 16000];
%movingBarLength = 40000;

nCells = length(cellIDs);
%nBarReps = length(chosenInts);
nBarReps = length(triggersStart);

%load all spike times from moving bar repetitions
spikeTimesAll = cell(nBarReps, 1);
interval = cell(nBarReps, 1);

%spikeTimesAll = datarun.spikes{get_cell_indices(datarun, cellIDs)};
for i = 1:nBarReps
    spikeTimesAll{i} = cell(nCells, 1);
    %interval{i} = [triggers(chosenInts(i)) triggers(chosenInts(i)+1)];
    interval{i} = [triggersStart(i)+0.45 triggersEnd(i)-0.4155];
    for j = 1:nCells
        spikeTimesTemp = datarun.spikes{get_cell_indices(datarun, cellIDs(j))};
        spikeTimesAll{i}{j} = spikeTimesTemp(spikeTimesTemp > interval{i}(1) & spikeTimesTemp < interval{i}(2)) - interval{i}(1); %sets start of interval to time = 0
        
    %spikeTimesAll{i} = getNeuronSpikeTimes(pathToNeuronsFile, cellIDs,
    %triggers(chosenInts(i))*20000, triggers(chosenInts(i)+1)*20000);
    end
end


%plot rasters
figure
maxNSpikes = 0;
for i = 1:nCells
    axes('position', [0.1 (nCells - i + 1)/(nCells+1) 0.8 1/(nCells+2)])
    hold on
    for k = 1:nBarReps
    %k=1
        for j = 1:length(spikeTimesAll{k}{i})
            plot([spikeTimesAll{k}{i}(j) spikeTimesAll{k}{i}(j)], [k-0.5 k+0.5], 'k-', 'LineWidth', 1)
        end
        maxNSpikes = max([maxNSpikes length(spikeTimesAll{k}{i})]);
    end
    hold off
    %set(gca, 'xLim', [0 nSamples/20000], 'yLim', [0 nBarReps])
    set(gca, 'xLim', [0 nSamples/20000])
    xlabel('time (s)')
    ylabel(['cell' num2str(cellIDs(i))])
end

%% save out as tab-delimited ascii file



spikeTimesCat = zeros(nCells*nBarReps, maxNSpikes);
for ii = 1:nBarReps
    for jj = 1:nCells
        nSpikesCurrent = length(spikeTimesAll{ii}{jj});
        (ii-1)*nCells + jj
        spikeTimesCat((ii-1)*nCells + jj, 1:nSpikesCurrent) = spikeTimesAll{ii}{jj};
    end
end

save('moving_bar_data_off', 'spikeTimesCat', '-ascii', '-tabs')

%% finds total costs between pairs of moving bar repetitions
% costs1 = zeros(nBarReps, nBarReps);
% costs10 = zeros(nBarReps, nBarReps);
% 
% dists1 = [];
% dists10 = [];
% 
% for i = 1:nBarReps
%     for j = i+1:nBarReps
%         sumCost1 = 0;
%         sumCost10 = 0;
%         for k = 1:nCells
%             d1 = spkd(spikeTimesAll{i}{k}, spikeTimesAll{j}{k}, 1/1000);
%             sumCost1 = sumCost1 + d1;
%             d10 = spkd(spikeTimesAll{i}{k}, spikeTimesAll{j}{k}, 10/1000);
%             sumCost10 = sumCost10 + d10;
%         end
%         dists1 = [dists1 sumCost1];
%         dists10 = [dists10 sumCost10];
%         costs1(i,j) = sumCost1;
%         costs10(i,j) = sumCost10;
%     end
% end
% 
% linkmap1 = linkage(dists1);
% linkmap10 = linkage(dists10);
% 
% figure
% dendrogram(linkmap1);
% title('cost parameter = 1 ms')
% 
% figure
% dendrogram(linkmap10);
% title('cost parameter = 10 ms')
% 
% 
%%


electrodes = cellElectrodes;
Array = zeros(length(electrodes), length(electrodes));
for i = 1:length(stimAmps)
    Array(i,i) = stimAmps(i);
end

%%
spikeTimesToUseIndeces = [17 40 25 34 20];

%for now, just use first response to moving bar
spikeTimesToUse = spikeTimesAll{spikeTimesToUseIndeces(1)};

%convert to samples
for i = 1:length(spikeTimesToUse)
    spikeTimesToUse{i} = spikeTimesToUse{i}*20000;
end

% 

% concatenate spike times for all cells and round down to nearest ms
spikeTimesOrigCat = [];
spikeTimesCat = [];
spikeTimes2Cat = [];
patterns = [];
%patterns2 = [];
for i = 1:length(spikeTimesToUse)
    spikeTimesOrig = floor(spikeTimesToUse{i}); %rounded to nearest sample
    spikeTimes = 20*floor(spikeTimesToUse{i}/20); %rounded to nearest 1 ms
    spikeTimes2 = 40*floor(spikeTimesToUse{i}/40); %rounded to nearest 2 ms
    
    spikeTimesOrigCat = [spikeTimesOrigCat; spikeTimesOrig];
    spikeTimesCat = [spikeTimesCat; spikeTimes]; %#ok<AGROW>
    spikeTimes2Cat = [spikeTimes2Cat; spikeTimes2]; %#ok<AGROW>
    patterns = [patterns; i*ones(length(spikeTimesToUse{i}),1)]; %#ok<AGROW>
    %patterns2 = [patterns2; i*ones(length(spikeTimesToUse{i}),1)];

end

%keyboard

% sorts spikeTimesCat and patterns by time
[spikeTimesCat iSort] = sort(spikeTimesCat);
patterns = patterns(iSort);



MovieChunks = 3; %2 chunks
Chunk = NS_MovieChunkGenerationForExperiment(spikeTimesOrigCat, nSamples, patterns);
MovieChunks = [MovieChunks Chunk];
Chunk = NS_MovieChunkGenerationForExperiment(spikeTimesCat, nSamples, patterns);
MovieChunks = [MovieChunks Chunk];
Chunk2 = NS_MovieChunkGenerationForExperiment(spikeTimes2Cat, nSamples, patterns);
MovieChunks = [MovieChunks Chunk2];    


MovieChunksFile=MovieChunks;

%%

fid = fopen('moving_bar1_off_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('moving_bar1_off_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('moving_bar1_off_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 