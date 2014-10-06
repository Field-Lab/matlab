clear;

NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ArrayID=500;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

PatternNumber=54;
MovieNumber=39;

PatternNumber=58;
MovieNumber=45;

DataPath='E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data003e';
ArtifactDataPath=DataPath;
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,PatternNumber,MovieNumber,120,0);

ClusterFileName='E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data003e\ClusterFile_003_n55';
WaveformTypes=NS_ReadClusterFile(ClusterFileName,MovieNumber,PatternNumber,200);
Artifacts=find(WaveformTypes==1);
Spikes=find(WaveformTypes==2);

ArtifactTraces=DataTraces(Artifacts,:,:);
Artifact=mean(ArtifactTraces);
TracesWithoutArtifact=NS512_SubtractArtifact(DataTraces,Artifact);

ChannelsPlot=[132 1:131 133:192];
[CorrectedTraces,EI,UniSpikesIndicCorrected]=NS512_TimingsForDetectedNeuron2(TracesWithoutArtifact(Spikes,:,:),ChannelsPlot);

a=CorrectedTraces(:,109,:);
a1=reshape(a,62,50);
figure(15)
plot(a1')
%break;
a=max(EI')-min(EI'); %peak to peak for each channel
BadChannels=find(a<4);
EI(BadChannels,:)=0;
WritePathFigs='C:\home\pawel\nauka\MEA2010\prezentacja\animacje\culture';
EI2=EI([2:132 1 133:192],:);
M=NS_SaveMovieFromSignature2(EI2,Channels,[452],ArrayID,WritePathFigs,NS_GlobalConstants);