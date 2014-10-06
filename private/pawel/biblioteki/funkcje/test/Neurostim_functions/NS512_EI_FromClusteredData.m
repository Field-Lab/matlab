function [EI,Channels]=NS512_EI_FromClusteredData(DataPath,ArtifactDataPath,ArtifactSubtraction,ClusterFilePath,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);

ClusterIndexes=NS_ReadClusterFileAll(ClusterFilePath);
SCI=size(ClusterIndexes);

NumberOfTraces=min(SCI(3),TracesNumberLimit);
WaveformTypes=ClusterIndexes(MovieNumber,PatternNumber,1:NumberOfTraces);
WaveformTypes=reshape(WaveformTypes,1,NumberOfTraces);
Falses=find(WaveformTypes==1);
Successes=find(WaveformTypes==2);

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
Artifact=mean(DataTraces(Falses,:,:));
Spike=mean(DataTraces(Successes,:,:));
EI=Spike-Artifact;