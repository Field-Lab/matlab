DataPath='D:\analysis\2008-08-26-0\data006_proba3';
ArtifactDataPath='D:\analysis\2008-08-26-0\data011_proba3';
ClusterFilePath='D:\analysis\2008-08-26-0\data006_proba3\clusters006_3';
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArtifactSubtraction=1;

MovieNumber=40;
PatternNumber=3;
Channels=[1:64];
[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
DataTraces=DataTraces(:,Channels,:);
WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
c1=find(WaveformTypes==1);
Traces=DataTraces(c1,:,:);
EI1=NS_CalculateEI(Traces);

MovieNumber=52;
PatternNumber=3;
Channels=[1:64];
[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
DataTraces=DataTraces(:,Channels,:);
WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
c1=find(WaveformTypes==2);
Traces=DataTraces(c1,:,:);
EI1=NS_CalculateEI(Traces)*1.5;

MovieNumber=61;
PatternNumber=1;
Channels=[1:64];
[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
DataTraces=DataTraces(:,Channels,:);
WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
c1=find(WaveformTypes==1);
Traces=DataTraces(c1,:,:);
EI2=NS_CalculateEI(Traces)/2;

MovieNumber=73;
PatternNumber=21;
Channels=[1:64];
[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
DataTraces=DataTraces(:,Channels,:);
WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
c1=find(WaveformTypes==1);
Traces=DataTraces(c1,:,:);
EI3=NS_CalculateEI(Traces)/1.5;

MovieNumber=70;
PatternNumber=60;
Channels=[1:64];
[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
DataTraces=DataTraces(:,Channels,:);
WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
c1=find(WaveformTypes==1);
Traces=DataTraces(c1,:,:);
EI4=NS_CalculateEI(Traces)/1.5;

MovieNumber=58;
PatternNumber=14;
Channels=[1:64];
[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
DataTraces=DataTraces(:,Channels,:);
WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
c1=find(WaveformTypes==1);
Traces=DataTraces(c1,:,:);
EI5=NS_CalculateEI(Traces);


EIsDiff(1,:,:)=reshape(EI1,1,64,60);
EIsDiff(2,:,:)=reshape(EI2,1,64,60);
EIsDiff(3,:,:)=reshape(EI3,1,64,60);
EIsDiff(4,:,:)=reshape(EI4,1,64,60);
EIsDiff(5,:,:)=reshape(EI5,1,64,60);

FigureProperties=struct('FigureNumber',102,'Subplot',[2 3 3],'TimeRange',[5 25],'AmplitudeRange',[-200 100],'FontSize',13,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'y'],'LineWidth',2,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(EIsDiff/0.44,Channels,[1 2 3 4 5],1,FigureProperties,NS_GlobalConstants);
%M=NS_SaveMovieFromSignature(EI4,Channels,1,FigureProperties,NS_GlobalConstants);


MovieNumber=46;
PatternNumber=3;
ArtifactSubtraction=1;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Channels=electrodeMap.getAdjacentsTo(PatternNumber,1)';
Movies=[40:3:61];
%Movies=40;
EIsDiff=zeros(length(Movies)+1,length(Channels),60);

Cluster2=[1 1 1 1 1 2 2 1];
Cluster1=[2 2 2 2 2 1 1 3];

for i=1:length(Movies)
    i
    MovieNumber=Movies(i);
    [DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
    DataTraces=DataTraces(:,Channels,:);
    WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
    %ArtifactDataTraces=ArtifactDataTraces(:,Channels,:);    
    c1=find(WaveformTypes==Cluster1(i));
    Traces=DataTraces(c1,:,:);
    EI1=NS_CalculateEI(Traces);
    
    c2=find(WaveformTypes==Cluster2(i));
    Traces=DataTraces(c2,:,:);
    EI2=NS_CalculateEI(Traces);
    EIsDiff(i,:,:)=EI2-EI1;
    Efficacy(i)=length(c2)
end

FileName='005000.bin';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-08-26-0\data005\data005.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-08-26-0\data005\data005.neurons');
idList = neuronFile.getIDList();
NeuronID=31;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-7;
cd D:\2008-08-26-0\data005;
FileName=['D:\2008-08-26-0\data005\data' FileName]
[RAWtraces,signal]=NS_AverageTraces(FileName,Timings1-1,Channels,[-2 57],NS_GlobalConstants);
signal=signal';
ss=size(signal);
for i=1:ss(1)
    signal(i,:)=signal(i,:)-mean([signal(i,ss(2)) signal(i,1)]);
end
EIsDiff(9,:,:)=signal;

FigureProperties=struct('FigureNumber',13,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-200 100],'FontSize',16,'Colors',['k' 'g' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2);

y=NS_PlotClustersOfSignaturesOnArrayLayout(EIsDiff/0.44,Channels,[1 1 1 1 1 1 1 1 0],1,FigureProperties,NS_GlobalConstants);
%y=NS_PlotClustersOfSignaturesOnArrayLayout(EIsDiff,Channels,ones(length(Movies),1),1,FigureProperties,NS_GlobalConstants);




PatternNumber=3;
ArtifactSubtraction=1;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Channels=electrodeMap.getAdjacentsTo(PatternNumber,1)';
Movies=[35:3:50];
%Movies=40;
clear EIsDiff;
EIsDiff=zeros(length(Movies)+1,length(Channels),60);

Cluster1=[1 1 1 1 1 1];
Cluster2=[2 2 2 2 2 5];

for i=1:length(Movies)
    i
    MovieNumber=Movies(i);
    [DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
    DataTraces=DataTraces(:,Channels,:);
    WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
    %ArtifactDataTraces=ArtifactDataTraces(:,Channels,:);    
    c1=find(WaveformTypes==Cluster1(i));
    Traces=DataTraces(c1,:,:);
    EI1=NS_CalculateEI(Traces);
    
    c2=find(WaveformTypes==Cluster2(i));
    Traces=DataTraces(c2,:,:);
    EI2=NS_CalculateEI(Traces);
    EIsDiff(i,:,:)=EI2-EI1;
    Efficacy(i)=length(c2)
end

FileName='005000.bin';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-08-26-0\data005\data005.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-08-26-0\data005\data005.neurons');
idList = neuronFile.getIDList();
NeuronID=31;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-9;
cd D:\2008-08-26-0\data005;
[RAWtraces,signal]=NS_AverageTraces(FileName,Timings1-1,Channels,[-2 57],NS_GlobalConstants);
signal=signal';
ss=size(signal);
for i=1:ss(1)
    signal(i,:)=signal(i,:)-mean([signal(i,ss(2)) signal(i,1)]);
end
EIsDiff(7,:,:)=signal;

FigureProperties=struct('FigureNumber',33,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-200 100],'FontSize',16,'Colors',['k' 'g' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2);

y=NS_PlotClustersOfSignaturesOnArrayLayout(EIsDiff/0.44,Channels,[1 1 1 1 1 1 0],1,FigureProperties,NS_GlobalConstants);




break;


MovieNumber=46;
PatternNumber=15;
ArtifactSubtraction=1;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Channels=electrodeMap.getAdjacentsTo(PatternNumber,1)';
Movies=[46:3:58];
%Movies=40;
clear EIsDiff;
EIsDiff=zeros(length(Movies)+1,length(Channels),60)
size(EIsDiff)

Cluster1=[2 1 1 1 1 2 2 1];
Cluster2=[1 2 2 2 2 1 1 3];

for i=1:length(Movies)
    i
    MovieNumber=Movies(i);
    [DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
    DataTraces=DataTraces(:,Channels,:);
    WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
    %ArtifactDataTraces=ArtifactDataTraces(:,Channels,:);    
    c1=find(WaveformTypes==Cluster1(i));
    Traces=DataTraces(c1,:,:);
    EI1=NS_CalculateEI(Traces);
    
    c2=find(WaveformTypes==Cluster2(i));
    Traces=DataTraces(c2,:,:);
    EI2=NS_CalculateEI(Traces);
    EIsDiff(i,:,:)=EI2-EI1;
    Efficacy(i)=length(c2)
end

FileName='000000.bin';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-08-26-0\data000\data000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-08-26-0\data000\data000.neurons');
idList = neuronFile.getIDList();
NeuronID=211;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-6;
cd D:\2008-08-26-0\data000;
[RAWtraces,signal]=NS_AverageTraces(FileName,Timings1-1,Channels,[-2 57],NS_GlobalConstants);
signal=signal';
ss=size(signal);
for i=1:ss(1)
    signal(i,:)=signal(i,:)-mean([signal(i,ss(2)) signal(i,1)]);
end

EIsDiff(5,:,:)=0;
EIsDiff(6,:,:)=signal;
FigureProperties=struct('FigureNumber',23,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',20,'Colors',['k' 'g' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2);

y=NS_PlotClustersOfSignaturesOnArrayLayout(EIsDiff,Channels,[1 1 1 1 1 0],1,FigureProperties,NS_GlobalConstants);
%y=NS_PlotClustersOfSignaturesOnArrayLayout(EIsDiff,Channels,ones(length(Mo
%vies),1),1,FigureProperties,NS_GlobalConstants);