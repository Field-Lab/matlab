DataPath='D:\analysis\2008-08-26-0\data006_proba3';
ArtifactDataPath='D:\analysis\2008-08-26-0\data011_proba3';
ClusterFilePath='D:\analysis\2008-08-26-0\data006_proba3\clusters006_3';
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArtifactSubtraction=1;

MovieNumber=76;
PatternNumber=22;
Channels=[1:64];
[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,55);

DataTraces=DataTraces(:,Channels,:);
WaveformTypes=NS_ReadClusterFile(ClusterFilePath,MovieNumber,PatternNumber);
c1=find(WaveformTypes==1);
Traces=DataTraces(c1,:,:);
EI=NS_CalculateEI(Traces);
figure(1)
clf
el=[28 34];
for i=1:length(el)
    figure(1)
    s=EI(el(i),:);
    plot(s)
    hold on;
end
grid on
%break;
FigureProperties=struct('FigureNumber',102,'Subplot',[2 3 3],'TimeRange',[5 25],'AmplitudeRange',[-200 100],'FontSize',13,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'y'],'LineWidth',2,'YLabel','input signal [\muV]');
FigureProperties.FigureNumber=105;
%M=NS_SaveMovieFromSignature(EI,Channels,1,FigureProperties,NS_GlobalConstants);
cd D:\analysis\2008-08-26-0\movies;
%movie2avi(M,['Stim_el' num2str(PatternNumber) 'ManyAmps'],'fps',15,'quality',100);

cd D:\2008-08-26-0\;
%[DataTraces,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber);
FileName='000';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-08-26-0\data000\data000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-08-26-0\data000\data000.neurons');
idList = neuronFile.getIDList();
NeuronID=346;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings=spikeTimes(1,1:min(50,length(spikeTimes)))-6;
channels=[11:24 26:29];
Offsets=ones(1,length(channels))*(-360);
Traces=NS_ReadManyTracesFromRaw(FileName,channels,Timings,0,80,Offsets,NS_GlobalConstants);

types=ones(1,length(Timings));

FigureProperties=struct('FigureNumber',103,'Subplot',[2 3 3],'TimeRange',[0 65],'AmplitudeRange',[-120 60],'FontSize',18,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'y'],'LineWidth',1,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(Traces/0.44,channels,types,1,FigureProperties,NS_GlobalConstants);

[DataTraces1,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,55);
[DataTraces2,ArtifactDataTraces,Channels1]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,70);
channels=[16:24 27:29];
DT=zeros(200,length(channels),60);
DT(1:100,:,:)=DataTraces1(:,channels,:);
DT(101:200,:,:)=DataTraces2(:,channels,:);

types(1:100)=1;
types(101:200)=2;
FigureProperties=struct('FigureNumber',113,'Subplot',[2 3 3],'TimeRange',[0 60],'AmplitudeRange',[-250 100],'FontSize',18,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'y'],'LineWidth',1,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(DT/0.44,channels,types,1,FigureProperties,NS_GlobalConstants);