function Amp=NS_2008_10_01(DataPath,ArtifactDataPath,ClusterFilePath,NS_GlobalConstants,ArtifactSubtraction,Movies,PatternNumber,Channels,electrode);

for i=1:length(Movies)
    MovieNumber=Movies(i)-3;
    [DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
    DataTraces=DataTraces(:,Channels,:);
    WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
    c1=find(WaveformTypes==1);
    Traces=DataTraces(c1,:,:);
    EI1=NS_CalculateEI(Traces);
    
    MovieNumber=Movies(i);
    [DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
    DataTraces=DataTraces(:,Channels,:);
    WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
    c1=find(WaveformTypes==1);
    Traces=DataTraces(c1,:,:);
    EI2=NS_CalculateEI(Traces);
    
    EIsDiff(i,:,:)=reshape(EI2,1,length(Channels),60);
end

eldata=EIsDiff(:,electrode,:);
seldata=size(eldata);
for i=1:seldata(1)
    d=eldata(i,1,:);
    Amp(i)=max(d)-min(d);
end