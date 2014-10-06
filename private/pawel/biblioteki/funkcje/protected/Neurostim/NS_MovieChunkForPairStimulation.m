function MovieChunks=NS_MovieChunkForPairStimulation(StimElectrodes,TimeShift,DelayBetweenPairs,DelaysInPair,NumberOfSamples);

NumberOfElectrodes=length(StimElectrodes);
TotalNumberOfPairs=NumberOfElectrodes*(NumberOfElectrodes-1)/2

[A,B]=NS_GeneratePairs(length(StimElectrodes));
%length(DelaysInPair)*DelayBetweenPairs
NumberOfPairsPerMovie=NumberOfSamples/(length(DelaysInPair)*DelayBetweenPairs)
NumberOfPairsPerMovie=floor(NumberOfPairsPerMovie);
NumberOfMovies=TotalNumberOfPairs/NumberOfPairsPerMovie

MovieChunks=[NumberOfMovies];
for i=1:NumberOfMovies
    Pairs=B(NumberOfPairsPerMovie*(i-1)+1:min(NumberOfPairsPerMovie*i,TotalNumberOfPairs),:);
    Patterns1=[];
    Times=[];
    for i=1:length(DelaysInPair) % for each value of relative delay
        Patterns1=[Patterns1 Pairs(:,1)' Pairs(:,2)'];
        Times0=TimeShift+DelayBetweenPairs*length(Pairs(:,1))*(i-1);
        TimesToAdd=[Times0:DelayBetweenPairs:Times0+DelayBetweenPairs*(NumberOfPairsPerMovie-1)];
        Times=[Times TimesToAdd TimesToAdd+DelaysInPair(i)];
    end
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,StimElectrodes(Patterns1));
    MovieChunks=[MovieChunks Chunk];
    
    %size(Times)
    %size(Patterns1)
end