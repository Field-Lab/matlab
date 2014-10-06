function [PatternOrder,MovieChunksFile]=NS_MovieChunksForExperiment2el(NumberOfClusters,TimeShiftInMs,DelayInMs,NumberOfSamples);

Pattern=zeros(NumberOfClusters,50*NumberOfClusters);

for i=1:NumberOfClusters
    Pattern(1,i:NumberOfClusters:i+NumberOfClusters*47)=[randperm(24) randperm(24)]+98*(i-1);
    Pattern(2,i:NumberOfClusters:i+NumberOfClusters*47)=[randperm(24) randperm(24)]+98*(i-1)+24;
    %Pattern(3,i:NumberOfClusters:i+NumberOfClusters*47)=[randperm(24) randperm(24)]+98*(i-1)+48;
    Pattern(3,i:NumberOfClusters:i+NumberOfClusters*49)=[97+98*(i-1) randperm(24)+98*(i-1)+48 97+98*(i-1) randperm(24)+98*(i-1)+48];
    %Pattern(4,i:NumberOfClusters:i+NumberOfClusters*47)=[randperm(24) randperm(24)]+98*(i-1)+72;
    Pattern(4,i:NumberOfClusters:i+NumberOfClusters*49)=[98+98*(i-1) randperm(24)+98*(i-1)+72 98+98*(i-1) randperm(24)+98*(i-1)+72];
end
TimeShift=TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay=round(DelayInMs*20/NumberOfClusters);

MovieChunks=[4];
for i=1:4
    l1=48*NumberOfClusters;
    l2=58*NumberOfClusters;
    if i<3
        l=48*NumberOfClusters;
    else
        l=50*NumberOfClusters;;
    end
    Patterns=Pattern(i,1:l);
    Times=[TimeShift:Delay:TimeShift+Delay*(l-1)];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
    MovieChunks=[MovieChunks Chunk];
end

PatternOrder=Pattern;
MovieChunksFile=MovieChunks;