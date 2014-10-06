clear;
Electrodes=[62:8:118 66:8:122 61:8:117 65:8:121 268:8:324 264:8:320 267:8:323 263:8:319 266:8:322 262:8:318 265:8:321 261:8:317];

NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ArrayID=500;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

PatternNumber=16;
MovieNumber=125;

DataPath='E:\pawel\analysis\retina\2009-11-27-0\data001\April2010\files';
ArtifactDataPath=DataPath;
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,PatternNumber,MovieNumber,100,0);

ClusterFileName='E:\pawel\analysis\retina\2009-11-27-0\data001\April2010\files\ClusterFile_001_ID_1264';
WaveformTypes=NS_ReadClusterFile(ClusterFileName,MovieNumber,PatternNumber,100);
Artifacts=find(WaveformTypes==1);
Spikes=find(WaveformTypes==2);

ArtifactTraces=DataTraces(Artifacts,:,:);
Artifact=mean(ArtifactTraces);
TracesWithoutArtifact=NS512_SubtractArtifact(DataTraces,Artifact);

ChannelsPlot=[85 1:84 86:512];
[CorrectedTraces,EI,UniSpikesIndicCorrected]=NS512_TimingsForDetectedNeuron3(TracesWithoutArtifact(Spikes,:,:),ChannelsPlot);

a=CorrectedTraces(:,1,:);
b=reshape(a,73,50);
figure(1)
plot(b')

EI2=EI([2:85 1 86:512],:);
%EI3=zeros(512,50);
%EI3(Electrodes,:)=EI2(Electrodes,:);
EI3=EI2(Electrodes,:);
WritePathFigs='C:\home\pawel\nauka\MEA2010\prezentacja\animacje\proba1';
M=NS_SaveMovieFromSignature2(EI3,Electrodes,[85 287],ArrayID,WritePathFigs,NS_GlobalConstants);