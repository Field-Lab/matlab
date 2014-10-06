function [PatternOrder, MovieChunksFile]=generateMovieClusterStim1Secondary(nClusters, clusterIDs, TimeShiftInMs, DelayInMs)

% generates movieChunksFile for stimulus file generation
%
% specifies order in which patterns are played during stimulus (time each pattern is played)

delay = DelayInMs*20/nClusters;

patternsInClusters = cell(nClusters, 1);
nPatternsInClusters = zeros(nClusters, 1);
patternOrders = cell(nClusters, 1);
for i = 1:nClusters;
    patternsInClusters{i} = find(clusterIDs == i);
    nPatternsInClusters(i) = sum(clusterIDs == i);
    
    patternOrders{i} = [randPerm(nPatternsInClusters(i)) randPerm(nPatternsInClusters(i))] + patternsInClusters{i}(1)-1;
end

nEvents = nClusters*max(nPatternsInClusters)*2; %factor of 2 is for 2 random orders

%interleaves the pattern orders
eventOrder = zeros(1,nEvents);
for i = 1:nClusters
    eventOrder(i:nClusters:(2*nPatternsInClusters(i)-1)*nClusters+i) = patternOrders{i};
end


nSamplesTotal = nEvents*delay;

nMovieChunks = ceil(nSamplesTotal/40000); %can only have 2 s in each movie chunk

if nMovieChunks > 1
    patternsPerMovieChunk = ceil(nEvents/nMovieChunks)*ones(nMovieChunks,1);
    patternsPerMovieChunk(end) = nEvents - patternsPerMovieChunk(1)*(nMovieChunks-1); %whatever's left over
else
    patternsPerMovieChunk = nEvents;
end

Pattern = zeros(nMovieChunks, max(patternsPerMovieChunk));
iEvent = 0;
for i = 1:nMovieChunks
    Pattern(i, 1:patternsPerMovieChunk(i)) = eventOrder(iEvent+1 : iEvent + patternsPerMovieChunk(i));
    iEvent = iEvent + patternsPerMovieChunk(i);
end


TimeShift = TimeShiftInMs*20; %in sampling periods (50 microseconds)
%Delay = round(DelayInMs*20/nClusters); %number of samples between patterns on different clusters


if nMovieChunks > 1
    nSamples = max(patternsPerMovieChunk)*delay+100;
    nSamples = 2000*ceil(nSamples/2000); %rounds up to nearest 0.1 seconds
else
    nSamples = 2000*ceil(nSamplesTotal/2000); %rounds up to nearest 0.1 seconds
end

warndlg(['set trigger interval to ' num2str(nSamples/20000) 's']);

% makes a movie chunk for each row of Pattern and concatenates MovieChunks
MovieChunks = nMovieChunks;
for i=1:nMovieChunks;
    patternsChunk = Pattern(i,:);
    times = 0:delay:delay*(length(patternsChunk)-1) + TimeShift;
    
    %removes patterns that are zeros in Pattern matrix
    for j = length(patternsChunk):-1:1
        if patternsChunk(j) == 0
            patternsChunk(j) = [];
            times(j) = [];
        end
    end
    Chunk = NS_MovieChunkGenerationForExperiment(times, nSamples, patternsChunk);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end

PatternOrder=Pattern;
MovieChunksFile=MovieChunks;

