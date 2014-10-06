function MovieChunksFile = generateMovieSpatiotempProbeNew(nCenters, timeShiftInMs, periodInMs, offsets, details)

% generates movieChunksFile for stimulus file generation
%
% offsets: a vector of offset values specifying delays between subsequent pulses within a sequence
% (in ms)
% specifies order in which patterns are played during stimulus (time each pattern is played)
%
% movie chunk length is limited to <=2 seconds (nSamples <= 40000)
%
% 2010-03-19: added second repetition of center-alone stimulus in final chunk so that total reps
% match other stimuli (since other stimuli have 2x random orders)
%
% note: in order for preprocessing to work correctly (correct parsing of
% raw data), each center electrode pulse amplitude must be in a separate
% movie chunk



prePulseBin = ~~details.prePulse;
stimPulseAmps = details.stimPulseAmps;
pPairsAmps = details.pairsCenterAmps;
clusterID = details.clusterID;


nOffsets = length(offsets);

nPrepulseAmpCombos = zeros(nCenters,1);
nPairedPulseCombosPerPrimAmp = zeros(nCenters,1);
for i = 1:nCenters
    nPrepulseAmpCombos(i) = sum(clusterID == i & prePulseBin); %number of prepulse-amplitude combinations in ith cluster
    nPairedPulseCombosPerPrimAmp(i) = sum(clusterID == i & pPairsAmps == 1); %assumes that there is
    %an equal number of paired pulse combinations for each primary electrode amplitude within cluster
end

maxNPrepulseAmpCombos = max(nPrepulseAmpCombos);
maxNPairedPulseCombosPerPrimAmp = max(nPairedPulseCombosPerPrimAmp);

%nPulsePairsPerPrimaryAmp = maxNPrepulseAmpCombos*nOffsets; % number of prepulse-primary pulse pairs for each cluster/stim pulse amplitude
%nStimuliPerPrimaryAmp = nPulsePairsPerPrimaryAmp + 1 + maxNPairedPulseCombosPerPrimAmp;


%timePerAmp = 20*periodInMs*nStimuliPerPrimaryAmp; %max per movie chunk = 40,000 (2 s)

offsetsPerChunk = floor(40000/(maxNPrepulseAmpCombos*2*periodInMs*20)); %number of delays that will fit into a 2s movie chunk

chunksPerAmp = ceil(nOffsets/offsetsPerChunk)+1; %+1 because first chunk has simultaneous and center-alone stim


if chunksPerAmp*max(stimPulseAmps) > 12
    warndlg([num2str(chunksPerAmp*max(stimPulseAmps)) ' movie chunks required'])
end

%put all the paired-pulse stimuli and the center-alone stimuli in the first
%movie chunk

if 20*periodInMs*2*(maxNPairedPulseCombosPerPrimAmp + 1) > 40000 %factor of 2 for 2 random orders
    error('current code isn''t flexible enough to handle this many surrounding electrode/prepulse amplitude combinations!')
end


period = periodInMs*20; %converts to samples

iChunk = 0;
for i = 1:max(stimPulseAmps)
   iChunk = iChunk+1;
   
   %first chunk for ith stimPulseAmp = paired-pulse stimuli and the center-alone stimuli
   pattern{iChunk} = []; %#ok<AGROW>
   times{iChunk} = []; %#ok<AGROW>
   for j = 1:nCenters
       centerAloneInd = find(stimPulseAmps == i & clusterID == j);
       pairsInClusterInd = find(clusterID == j & pPairsAmps == i);
       if ~isempty(pairsInClusterInd)
           nPairs = length(pairsInClusterInd); %should be the same as nPairedPulseCombosPerPrimAmp(j)
           if nPairs ~= nPairedPulseCombosPerPrimAmp(j)
               keyboard
           end
           pattern{iChunk} = [pattern{iChunk} pairsInClusterInd(randperm(nPairs)) centerAloneInd]; %#ok<AGROW> %first random order
           pattern{iChunk} = [pattern{iChunk} pairsInClusterInd(randperm(nPairs)) centerAloneInd]; %#ok<AGROW> %second random order
           timesTemp = period:period:period*2*(nPairs+1);
           times{iChunk} = [times{iChunk} timesTemp + (j-1)*period/nCenters]; %#ok<AGROW,*AGROW>
       end
   end
       
   iDelay = 0;
   for j = 2:chunksPerAmp
       iChunk = iChunk + 1;
       currentTime = period;
       pattern{iChunk} = []; %#ok<AGROW>
       times{iChunk} = []; %#ok<AGROW>
       for k = 1:offsetsPerChunk
           iDelay = iDelay + 1;
           
           if iDelay <= length(offsets)
               for m = 1:nCenters
                   if any(stimPulseAmps == i & clusterID == m)
                       prepulseInd = find(prePulseBin & clusterID == m);
                       nPulsePairs = length(prepulseInd);
                       centerAloneInd = find(stimPulseAmps == i & clusterID == m);
                       
                       pattern{iChunk} = [pattern{iChunk} prepulseInd(randperm(nPulsePairs))]; %#ok<AGROW> %first random order of prepulses
                       pattern{iChunk} = [pattern{iChunk} prepulseInd(randperm(nPulsePairs))]; %#ok<AGROW> %second random order of prepulses
                       timesTemp = currentTime:period:currentTime + period*2*nPulsePairs - 1;
                       times{iChunk} = [times{iChunk} timesTemp - offsets(iDelay)*20 + (m-1)*period/nCenters]; %#ok<AGROW>
                       pattern{iChunk} = [pattern{iChunk} centerAloneInd*ones(1,2*nPulsePairs)]; %#ok<AGROW> %stim pulses
                       times{iChunk} = [times{iChunk} timesTemp + (m-1)*period/nCenters]; %#ok<AGROW>
                   end
               end
               currentTime = period*floor(max(times{iChunk})/period) + period;
           end
       end
   end
end


%orders movies and times, and determines trigger interval
maxMovieLength = 0;
for i = 1:length(pattern) %loops through movie chunks
    [times{i} patternOrder] = sort(times{i}, 2, 'ascend'); %#ok<AGROW>
    pattern{i} = pattern{i}(patternOrder); %#ok<AGROW>
    
    maxMovieLength = max([maxMovieLength times{i}(end)]);
end
maxMovieLength = maxMovieLength + period;

nSamples = 2000*ceil(maxMovieLength/2000); %rounds up to the nearest 0.1 seconds
warndlg(['set trigger interval to ' num2str(nSamples/20000) ' seconds'])


% makes a movie chunk for each row of Pattern and concatenates MovieChunks
MovieChunks = length(pattern);
for i = 1:length(pattern)
    patternChunk = pattern{i};
    timesChunk = times{i} + timeShiftInMs*20;
    
    Chunk = NS_MovieChunkGenerationForExperiment(timesChunk, nSamples, patternChunk);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end
    
MovieChunksFile=MovieChunks;