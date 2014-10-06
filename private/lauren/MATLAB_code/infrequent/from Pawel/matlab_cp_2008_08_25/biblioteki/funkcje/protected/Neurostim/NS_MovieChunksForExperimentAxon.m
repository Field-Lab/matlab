function [PatternsOut,Times,MovieChunk]=NS_MovieChunksForExperimentAxon(NumberOfClusters,Array,TimeShiftInMs,DelayInMs,NumberOfSamples);

AS=size(Array);
NumberOfClusters=round(AS(1)/7)
ClusterLength=round(AS(2)/NumberOfClusters)

Patterns=zeros(1,ClusterLength*NumberOfClusters);

for i=1:NumberOfClusters
    StartIndex=(i-1)*ClusterLength;
    Patterns(1,i:NumberOfClusters:i+NumberOfClusters*(ClusterLength-1))=randperm(ClusterLength)+(i-1)*ClusterLength;
end

TimeShift=TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay=round(DelayInMs*20/NumberOfClusters);

Times=[];
PatternsOut=[];
for i=1:length(Patterns)
    Time=TimeShift+(i-1)*Delay;        
    Multipliers=Array(:,Patterns(i));
    if max(max(abs(Multipliers)))>0
        Times=[Times Time];
        PatternsOut=[PatternsOut Patterns(i)];
    else
        i
    end
end
size(Times)

MovieChunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsOut);