function MovieChunks=NS_MovieChunksForExperimentFrequencyScan(PatternNumber,Periods,DelayInMs,NumberOfSamples);

MovieChunks=[length(Periods)]; % starting value
for i=1:length(Periods)
    Times=[DelayInMs:Periods(i):NumberOfSamples];
    Patterns=ones(1,length(Times))*PatternNumber;
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
    MovieChunks=[MovieChunks Chunk];
end