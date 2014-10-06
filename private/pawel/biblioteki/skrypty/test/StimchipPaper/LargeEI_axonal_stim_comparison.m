FileName='E:\pawel\analysis\retina\2010-09-21-0\data003\ClusterFile_003_id227_2';
MovieNumber=22;
PatternNumber=56;
WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber,100);

DataPath='E:\pawel\analysis\retina\2010-09-21-0\data003';
ArtifactDataPath='E:\pawel\analysis\retina\2010-09-21-0\data003';
[DataTraces2,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,PatternNumber,MovieNumber,100,0);



TraceisWithArtifacts=find(WaveformTypes==1);
TracesWithSpikes=find(WaveformTypes==2);

Artifact=mean(DataTraces2(TraceisWithArtifacts,:,:));

TracesWithoutArtifact=NS512_SubtractArtifact(DataTraces2(TracesWithSpikes,:,:),Artifact);

ChannelsPlot=[8 10 13 16 5 7 11 14 18 2 3 6 12 17 64 1 4 58 60 61 54 56 53 55 19 49 47 41 15 44];
[CorrectedTraces,EI,TimingsForSpikesMinimumsCorrected]=NS512_TimingsForDetectedNeuron3(TracesWithoutArtifact,ChannelsPlot);