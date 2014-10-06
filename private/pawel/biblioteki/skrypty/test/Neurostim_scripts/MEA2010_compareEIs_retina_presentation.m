NS_GlobalConstants=NS_GenerateGlobalConstants(512);

% * * * * * Spontaneous EI * * * * 
FileName='E:\pawel\data\retina\2009-11-27-0\data000';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\retina\2009-11-27-0\data000\data000000\data000000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\retina\2009-11-27-0\data000\data000000\data000000.neurons');
idList = neuronFile.getIDList();
NeuronID=3631;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
N=500;
L=100;
Timings1=spikeTimes(1,1:min(N,length(spikeTimes)))-14;
Data=zeros(513,L);
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(FileName); 
for i=1:length(Timings1)
    RawData=double(rawFile.getData(Timings1(i),L)');    
    Data=Data+RawData;
end
m=mean(Data')';
for i=1:513
    Data(i,:)=Data(i,:)-m(i);
end

T(2,:,:)=Data(2:513,:)/N;

%  * * * * * Elicited EI  * * * * 
FilePath='F:\analysis\retina\2009-11-27-0\data001\ClusterFile_001';
ClusterIndex0=NS_ReadClusterFile(FilePath,58,42,200);
ClusterIndex=ClusterIndex0(1:100);

DataPath='F:\analysis\retina\2009-11-27-0\data001';
[DataTraces0,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,42,58,200,0);
DataTraces=DataTraces0(1:100,:,:);

artifacts=find(ClusterIndex==1);
%artifacts=find(artifacts0<103);
spikes=find(ClusterIndex==2);
clear Artifact;
Artifact=mean(DataTraces(artifacts,:,:));
S=mean(DataTraces(spikes,:,:));
Spike=S-Artifact;
m=mean(Spike,3);
for i=1:192
    Spike(1,i)=Spike(1,i)-m(i);
end
T(1,:,:)=Spike;
%Chns=[444 445 451 452 453 460 461];
Chns=[233 234 241 242 243 249 250];

% * * * * * Plot  * * * * 
T1=T(:,Chns,:);
T1=T1./0.27;
Chns=[233 234 241 242 243 249 250];
FigureProperties=struct('FigureNumber',5,'TimeRange',[0 59],'AmplitudeRange',[-600 200],'FontSize',24,'Colors',['b' 'r' 'k' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',2,'YLabel','Signal [\muV]');
% both EIs:
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3Rotated(T1,[Chns],[1 2],500,FigureProperties,NS_GlobalConstants,[144],[]);
% only stimulated EI:
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3Rotated(T1(1,:,:),[Chns],[1],500,FigureProperties,NS_GlobalConstants,[144],[]);

h=gcf;
FullName=['C:\home\pawel\nauka\MEA2010\prezentacja\figures_pok107\' 'EI_retina_stim_only_rotated'];            
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 7]);
set(h,'PaperPosition',[0 0 10 7]);  
print(h, '-dtiff', '-r120', FullName);