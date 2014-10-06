function [PatternOrder,MovieChunksFile]=NS_MovieChunksForExperiment2elModified(NumberOfClusters,TimeShiftInMs,DelayInMs,NumberOfSamples);

Pattern=zeros(5,50*NumberOfClusters);
PatternsPerCluster=110;

for i=1:NumberOfClusters
    offset=PatternsPerCluster*(i-1);
    Pattern(1,i:NumberOfClusters:i+NumberOfClusters*47)=[randperm(24) randperm(24)]+offset; % center with 0.25 for perifery el.
    Pattern(2,i:NumberOfClusters:i+NumberOfClusters*47)=[randperm(24) randperm(24)]+offset+24; % center with 0.5 for perifery el.
    Pattern(3,i:NumberOfClusters:i+NumberOfClusters*49)=[97+offset randperm(24)+offset+48 97+offset randperm(24)+offset+48]; % center with 0.75 for perifery el., plus center alone positive
    Pattern(4,i:NumberOfClusters:i+NumberOfClusters*49)=[98+offset randperm(24)+offset+72 98+offset randperm(24)+offset+72]; % center with 1 for perifery el., plus center alone negative
    Pattern(5,i:NumberOfClusters:i+NumberOfClusters*47)=[randperm(12) randperm(12) randperm(12) randperm(12)]+offset+98; % only perifery electrodes (both polarities)
end
TimeShift=TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay=round(DelayInMs*20/NumberOfClusters);

MovieChunks=[5]
for i=1:5
    if i<3
        l=48*NumberOfClusters;
    else
        l=50*NumberOfClusters;
    end
    if i==5
        l=48*NumberOfClusters;
    end
    Patterns=Pattern(i,1:l);
    Times=[TimeShift:Delay:TimeShift+Delay*(l-1)];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
    MovieChunks=[MovieChunks Chunk];
end

PatternOrder=Pattern;
MovieChunksFile=MovieChunks;