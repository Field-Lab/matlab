function MovieChunksFile = generateMovieAxonStim(NumberOfClusters, Array, TimeShiftInMs, DelayInMs, NumberOfSamples)

AS=size(Array);
NumberOfClusters=round(AS(1)/7);
ClusterLength=round(AS(2)/NumberOfClusters);

nChunks = 2;

Patterns=zeros(nChunks,ClusterLength*NumberOfClusters);

for i=1:NumberOfClusters
    StartIndex=(i-1)*ClusterLength;
    Patterns(1,i:NumberOfClusters:i+NumberOfClusters*(ClusterLength-1))=randperm(ClusterLength)+(i-1)*ClusterLength;
    Patterns(2,i:NumberOfClusters:i+NumberOfClusters*(ClusterLength-1))=randperm(ClusterLength)+(i-1)*ClusterLength;
end


TimeShift=TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay=round(DelayInMs*20/NumberOfClusters);

MovieChunks = [];
for j = 1:nChunks
    Times=[];
    PatternsOut=[];
    for i=1:size(Patterns,2)
        Time=TimeShift+(i-1)*Delay;
        Multipliers=Array(:,Patterns(j,i));
        if max(max(abs(Multipliers)))>0
            Times = [Times Time];
            PatternsOut = [PatternsOut Patterns(j,i)];
        else
            i
        end
    end

    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsOut);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW>
end


MovieChunksFile = [nChunks MovieChunks];%first value indicates number of chunks