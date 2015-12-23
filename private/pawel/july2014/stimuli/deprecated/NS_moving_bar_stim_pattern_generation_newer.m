
clear all

pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data002-from-data000/data002-from-data000.neurons';

nSamples = 20000; %make sure this is enough time for moving bar responses

nBarReps = 40; % number of moving bar repetitations for each color in visual stimulus
barColor = 'all'; % 'white', 'black', or 'all'

% cellIDs =        [35   181  407  511  683  814];
% cellElectrodes = [64   13   22   33   49   59];
% stimAmps =       [.4   1.6  1.1  1.8  0.7  0.9];



cellIDs =        [36   48   247  393  531  694];
cellElectrodes = [2    4    17   21   36   47];
stimAmps =       [1.1  0.68 .35  1.1  .55  0.9];


startPadding = 0.45; %how long after start trigger response is visible (in seconds)
endPadding = 0.4155; %how long before end trigger reponses seems to end (in seconds)

% %% get triggers to determine start/end times of moving bar repetitions
% 
% datarun.names.rrs_neurons_path = pathToNeuronsFile;
% datarun = load_neurons(datarun);
% 
% triggers = datarun.triggers;
% 
% nTriggers = length(triggers);
% 
% if (strcmpi(barColor, 'all') && nTriggers ~= 2*nBarReps) || (~strcmpi(barColor, 'all') && nTriggers ~= 4*nBarReps)
%     warndlg('unexpected number of triggers--check that nReps value is correct')
% end
% 
% if strcmpi(barColor, 'white')
%     triggersStart = triggers(1:4:nBarReps*4-3); %start of moving bar
%     triggersEnd = triggers(2:4:nBarReps*4-2); %end of moving bar
% elseif strcmpi(barColor, 'black')
%     triggersStart = triggers(3:4:nBarReps*4-1); %start of moving bar
%     triggersEnd = triggers(4:4:nBarReps*4); %end of moving bar
% elseif strcmpi(barColor, 'all')
%     triggersStart = triggers(1:2:nBarReps*2-1); %start of moving bar
%     triggersEnd = triggers(2:2:nBarReps*2); %end of moving bar
% else
%     errordlg('invalid bar color specified')
% end
% 
% 
% if min(triggersEnd - triggersStart) < 1
%     warndlg('moving bar takes less than 1 second--check that triggers are identified correctly')
% end
% 
% nCells = length(cellIDs);
% 
% 
% %load all spike times from moving bar repetitions
% spikeTimesAll = cell(nBarReps, 1);
% interval = cell(nBarReps, 1);
% 
% for i = 1:nBarReps
%     spikeTimesAll{i} = cell(nCells, 1);
%     interval{i} = [triggersStart(i)+startPadding triggersEnd(i)-endPadding];
%     for j = 1:nCells
%         spikeTimesTemp = datarun.spikes{get_cell_indices(datarun, cellIDs(j))};
%         % retrieves spike times within window of time
%         spikeTimesAll{i}{j} = spikeTimesTemp(spikeTimesTemp > interval{i}(1) & spikeTimesTemp < interval{i}(2));
%         spikeTimesAll{i}{j} = spikeTimesAll{i}{j} - interval{i}(1); % sets start of time window to time = 0
%     end
% end
% 
% %%
% 
% %plot rasters
% figure
% 
% for i = 1:nCells
%     axes('position', [0.1 (nCells - i + 1)/(nCells+1) 0.8 1/(nCells+2)])
%     hold on
%     for k = 1:nBarReps
%         for j = 1:length(spikeTimesAll{k}{i})
%             plot([spikeTimesAll{k}{i}(j) spikeTimesAll{k}{i}(j)], [k-1 k], 'k-', 'LineWidth', 1)
%         end
%     end
%     hold off
%     
%     set(gca, 'yLim', [0 nBarReps], 'xlim', [0 interval{1}(2) - interval{1}(1)])
%     if i == nCells
%         xlabel('time (s)')
%     else
%         set(gca, 'xtick', [])
%     end
%     ylabel(['cell' num2str(cellIDs(i))])
% end
% 
% %% sniff tests
% 
% 
% %add spikes at same time to all repetitions of originals
% for i = 1:nCells
%     randSpikes = 0.8*rand(5+i,1);
%     for k = 1:nBarReps
%         spikeTimesAll{k}{i} = [spikeTimesAll{k}{i}; randSpikes];
%     end
% end
% 
% %% compute pair-wise Victor spike train distances (costs)
% 
% costParam = 10/1000;
% %costParam = 10/1000; %in seconds; half of the distance that a spike moves that is equivalent (in cost) to deleting and inserting spike
% 
% costParam = 1/costParam; %cost per second to move spike
% 
% costs = zeros(nBarReps, nBarReps);
% 
% for i = 1:nBarReps
%     for j = i+1:nBarReps
%         sumCost = 0;
%         for k = 1:nCells %calculate pairwise cost for each cell and sum
%             d = spkd(spikeTimesAll{i}{k}, spikeTimesAll{j}{k}, costParam);
%             sumCost = sumCost + d;
%         end
%         costs(i,j) = sumCost;
%     end
% end
% 
% %%
% 
% if 1 %plots linkage map
%     %converts matrices to "pdist" format
%     dists = [];
%     for i = 1:nBarReps
%         dists = [dists costs(i,i+1:end)];
%     end
% 
%     linkmap = linkage(dists);
% 
%     figure
%     dendrogram(linkmap, 0);
% end
% 
% %% for each moving bar response, calculates mean distance to all other moving bar responses
% 
% costs = costs + costs'; %converst from URT to symmetric matrix
% meanCosts = mean(costs,1);
% 
% if 1 %plots histogram of mean distances
%     figure
%     hist(meanCosts)
% end
% 
% %%
% 
% [y order] = sort(meanCosts);
% 
% %% plot raster of subset of reps
% 
% repsToPlot = 1:40;
% %repsToPlot = order(1:10);
% 
% figure
% for ii = 1:nCells
%     axes('position', [0.1 (nCells - ii + 0.8)/(nCells+1) 0.8 1/(nCells+3)])
%     hold on
%     for kk = 1:length(repsToPlot)
%         for jj = 1:length(spikeTimesAll{repsToPlot(kk)}{ii})
%             plot([spikeTimesAll{repsToPlot(kk)}{ii}(jj) spikeTimesAll{repsToPlot(kk)}{ii}(jj)],...
%                 [kk-1 kk], 'k-', 'LineWidth', 1)
%         end
%     end
%     hold off
%     set(gca, 'yLim', [0 length(repsToPlot)], 'xlim', [0 interval{1}(2) - interval{1}(1)])
%     
%     if ii == nCells
%         xlabel('time (s)')
%     else
%         set(gca, 'xtick', [])
%     end
%     
%     ylabel(['cell' num2str(cellIDs(ii))])
%     if ii == 1
%         title('10 ''best'' trials')
%     end
% end
% 
% %%
% 
% repsToPlot = order(end-9:end);
% 
% figure
% for ii = 1:nCells
%     axes('position', [0.1 (nCells - ii + 0.8)/(nCells+1) 0.8 1/(nCells+3)])
%     hold on
%     for kk = 1:length(repsToPlot)
%         for jj = 1:length(spikeTimesAll{repsToPlot(kk)}{ii})
%             plot([spikeTimesAll{repsToPlot(kk)}{ii}(jj) spikeTimesAll{repsToPlot(kk)}{ii}(jj)],...
%                 [kk-1 kk], 'k-', 'LineWidth', 1)
%         end
%     end
%     hold off
%     set(gca, 'yLim', [0 length(repsToPlot)], 'xlim', [0 interval{1}(2) - interval{1}(1)])
%     xlabel('time (s)')
%     ylabel(['cell' num2str(cellIDs(ii))])
%     
%     if ii == nCells
%         xlabel('time (s)')
%     else
%         set(gca, 'xtick', [])
%     end
%     if ii == 1
%         title('10 ''worst'' trials')
%     end
% end


%% generate stimulus files


electrodes = cellElectrodes;
Array = zeros(length(electrodes), length(electrodes));
for i = 1:length(stimAmps)
    Array(i,i) = stimAmps(i);
end

%%
%spikeTimesToUseIndeces = [38 18 15 30 3];

%for now, just use first response to moving bar
%spikeTimesToUse = spikeTimesAll{spikeTimesToUseIndeces(5)};

%convert to samples
for i = 1:length(spikeTimesToUse)
    for j = 1:length(spikeTimesToUse{i}
        spikeTimesToUse{i}{j} = spikeTimesToUse{i}*20000;
    end
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

fid = fopen('moving_bar5_update_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('moving_bar5_update_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('moving_bar5_update_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 