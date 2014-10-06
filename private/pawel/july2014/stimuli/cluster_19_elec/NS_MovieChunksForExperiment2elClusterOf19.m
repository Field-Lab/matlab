function [PatternOrder,MovieChunksFile]=NS_MovieChunksForExperiment2elClusterOf19(NumberOfClusters,TimeShiftInMs,DelayInMs,NumberOfSamples);

Pattern=zeros(5,42);
N=41;

Pattern(1,1:2*N)=[randperm(N) randperm(N)];
Pattern(2,1:2*N)=[randperm(N) randperm(N)]+N;
Pattern(3,1:2*N)=[randperm(N) randperm(N)]+2*N;
Pattern(4,1:2*N)=[randperm(N) randperm(N)]+3*N;
Pattern(5,1:2*(N+1))=[randperm(N+1) randperm(N+1)]+4*N;

TimeShift=TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay=DelayInMs*20;

MovieChunks=[5]
for i=1:5
    if i<5
        l=2*N;
    else 
        l=2*(N+1);
    end
    Patterns=Pattern(i,1:l);
    Times=[TimeShift:Delay:TimeShift+Delay*(l-1)];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
    MovieChunks=[MovieChunks Chunk];
end

PatternOrder=Pattern;
MovieChunksFile=MovieChunks;