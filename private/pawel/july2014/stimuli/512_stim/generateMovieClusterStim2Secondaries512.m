function MovieChunksFile = generateMovieClusterStim2Secondaries512(clusterIDs, timeShiftInMs, delayInMs, smallAmpPatterns)


% usage:  MovieChunksFile = generateMovieClusterStim2Secondaries512(clusterIDs, timeShiftInMs, delayInMs)
%
% generates movieChunksFile for stimulus file generation
%       - patterns associated with different electrode clusters are interleaved
%       - maximum number of patterns associated with a single cluster is used to determine number of
%         movie chunks to break movie into and the number of samples required in each movie chunk
%       - patterns associated with a given cluster are randomly grouped into sets of patterns that
%         will fit into single movie chunks, and within each movie chunk the set of patterns is
%         applied in two random orders
%         * each pattern is applied twice (actual number of pattern repetitions will be DOUBLE the number specified
%         in labview)
%
% arguments:     clusterIDs - a vector of integers specifying which patterns correspond with
%                             which clusters; each cluster is given a unique integer value, and the vector
%                             indeces correspond with the pattern numbers
%
%                 delayInMs - defines length of time between pattern application on SAME cluster
%             timeShiftInMs - defines time offset before starting stimulus application (typically 0)
%                             *code untested for nonzero values (may cause problems with choosing
%                             movie chunk lengths)
%          smallAmpPatterns - patterns that should be separated out into their own movie chunk to 
%                             avoid current range conflicts (see
%                             generatePatternClusterStim2Secondaries.m)
%
%
% outputs:  MovieChunksFile - array of values that specifies when each stimulus pattern is played, in the
%                             format read by the Labview Software
%
% author: Lauren (SNL-E), based on function written by Pawel Hottowy
% 2009-08-24
% checked/commented 2010-10-13

% 2014-07-07 edited for STIM512, Lauren Grosberg (line 67, max nSamples per
% trigger interval is 10000 for 512 system, divide by 9000 due to rounding errors later)


%convert inputs to samples
delay = delayInMs*20;
timeShift = timeShiftInMs*20; %in sampling periods (50 microseconds)

nClusters = max(clusterIDs); %number of electrode "clusters" (primary + set of secondaries)

%% determine how many chunks (nChunksTotal) are needed and how long each chunk will be (nSamples)

nSmallAmpPatternsEachCluster = zeros(nClusters, 1);
for ii = 1:nClusters
    nSmallAmpPatternsEachCluster(ii) = sum(clusterIDs == ii & smallAmpPatterns);
end

maxPatternsInCluster = 0; %maximum number of patterns associated with a single electrode cluster
maxPatternsInClusterSmallAmps = 0;
randGroupings = cell(nClusters, 1); %random ordering of patterns within each cluster
randGroupingsSmallAmps = cell(nClusters, 1);
for i = 1:nClusters
    maxPatternsInCluster = max(maxPatternsInCluster, sum(clusterIDs==i & ~smallAmpPatterns)); %only count the normal patterns (not small amps patterns)
    maxPatternsInClusterSmallAmps = max(maxPatternsInClusterSmallAmps, sum(clusterIDs==i & smallAmpPatterns));
    patternsInCluster = find(clusterIDs==i & ~smallAmpPatterns);
    patternsInClusterSmallAmps = find(clusterIDs==i & smallAmpPatterns);
    randGroupings{i} = patternsInCluster(randperm(length(patternsInCluster))); %so that pattern numbers associated with each cluster are randomly divided between movie chunks
    randGroupingsSmallAmps{i} = patternsInClusterSmallAmps(randperm(length(patternsInClusterSmallAmps)));
end

%patterns on different electrode clusters are interleaved, so length of time required is set by
%cluster with largest number of associated patterns
nChunks = ceil((timeShift + maxPatternsInCluster*2*delay)/9000); %factor of 2 is for 2 random orders
nPatternsPerChunk = ceil(maxPatternsInCluster/nChunks);

%how many extra chunks are required for small-amp patterns
nChunksSmallAmps = ceil(maxPatternsInClusterSmallAmps/nPatternsPerChunk);

nChunksTotal = nChunks + nChunksSmallAmps;

%patterns are divided roughly evenly between movie chunks (for cluster with most patterns)
nSamples = ceil((timeShift + nPatternsPerChunk*2*delay)/2000)*2000; %division/multiplication used to round to nearest tenth of a second

h = msgbox(['Number of movie chunks required = ' num2str(nChunksTotal), 10, 10, 'Set trigger interval to ' num2str(nSamples/20000) ' seconds']);
boxPos = get(h, 'position');
set(h, 'position', [400 500 boxPos(3) boxPos(4)])

%% generates lists of pattern ids and application times for each movie chunk

pattern = cell(nChunksTotal, 1);
times = cell(nChunksTotal, 1);

patternInd = 0; %keeps track of what pattern number index to start with for a given movie chunk
for i = 1:nChunks
    for j = 1:nClusters
        %determines pattern numbers to be used in this movie chunk for this cluster
        if length(randGroupings{j}) >= patternInd + nPatternsPerChunk %not the last movie chunk that includes patterns from this cluster
            pToUse = randGroupings{j}(patternInd + 1:patternInd + nPatternsPerChunk);
        else %last movie chunk that includes patterns from this cluster
            pToUse = randGroupings{j}(patternInd + 1:end);
        end
        %adds pattern numbers in 2 random orders to pattern{i} and associated pattern application
        %times to times{i}
        pattern{i} = [pattern{i} pToUse(randperm(length(pToUse)))...
            pToUse(randperm(length(pToUse)))]; %#ok<AGROW>
        timesTemp = 0:delay:delay*(2*length(pToUse)-1);
        times{i} = [times{i} timesTemp + (j-1)*delay/nClusters];
    end
    patternInd = patternInd+nPatternsPerChunk;
end

patternInd = 0;
for i = (1:nChunksSmallAmps) + nChunks
    for j = 1:nClusters
        %determines pattern numbers to be used in this movie chunk for this cluster
        if length(randGroupingsSmallAmps{j}) >= patternInd + nPatternsPerChunk %not the last movie chunk that includes patterns from this cluster
            pToUse = randGroupingsSmallAmps{j}(patternInd + 1:patternInd + nPatternsPerChunk);
        else %last movie chunk that includes patterns from this cluster
            pToUse = randGroupingsSmallAmps{j}(patternInd + 1:end);
        end
        %adds pattern numbers in 2 random orders to pattern{i} and associated pattern application
        %times to times{i}
        pattern{i} = [pattern{i} pToUse(randperm(length(pToUse)))...
            pToUse(randperm(length(pToUse)))]; %#ok<AGROW>
        timesTemp = 0:delay:delay*(2*length(pToUse)-1);
        times{i} = [times{i} timesTemp + (j-1)*delay/nClusters];
    end
    patternInd = patternInd+nPatternsPerChunk;
end

%sorts patterns and times by time and adds timeShift
for i = 1:nChunksTotal
    [times{i} patternOrder] = sort(times{i}, 2, 'ascend'); %#ok<AGROW>
    times{i} = times{i} + timeShift;
    pattern{i} = pattern{i}(patternOrder); %#ok<AGROW>
end

%checks that each pattern is included exactedly twice in one of the movie chunks
for i = 1:length(clusterIDs)
    nReps = zeros(length(pattern),1);
    for j = 1:length(pattern)
        nReps(j) = sum(pattern{j}==i);
    end
    if ~sum(nReps==2) == 1 || ~sum(nReps==0) == length(pattern)-1
        warndlg(['pattern' num2str(i) 'doesn''t appear twice in one movie chunk and 0 times in all other chunks'])
    end
end

%% generates movie file in labview-readable format using specified patterns (pattern) and applications times (times)

MovieChunks = nChunksTotal; %first value indicates number of chunks
for i = 1:nChunksTotal
    chunkPatterns = pattern{i};
    chunkTimes = times{i};
    Chunk = NS_MovieChunkGenerationForExperiment(chunkTimes, nSamples, chunkPatterns);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> 
end

MovieChunksFile=MovieChunks;


