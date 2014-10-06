NS_GlobalConstants=NS_GenerateGlobalConstants(512);

% * * * * * Spontaneous EI * * * * 
FileName='E:\pawel\data\cultures\2009-11-20-0\data002';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data002\data002.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data002\data002.neurons');
idList = neuronFile.getIDList();
NeuronID=6647;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
N=500;
L=100;
Timings1=spikeTimes(1,1:min(N,length(spikeTimes)))-17;
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
FilePath='E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data003e\ClusterFile_003_n65';
ClusterIndex0=NS_ReadClusterFile(FilePath,54,63,200);
ClusterIndex=ClusterIndex0(1:100);

DataPath='E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data003e\';
[DataTraces0,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,63,54,200,0);
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
fff=T(1,321:512,1:100);
T(1,321:512,:)=Spike(1,:,1:100)*2.6;

Chns=[444 445 451 452 453 460 461];

% * * * * * Plot  * * * * 
T1=T(:,Chns,:);
T1=T1./0.84;
Chns=[233 234 241 242 243 249 250];
FigureProperties=struct('FigureNumber',6,'TimeRange',[0 56],'AmplitudeRange',[-80 40],'FontSize',24,'Colors',['b' 'r' 'k' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',2,'YLabel','Signal [\muV]');
%y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks4(T1,[Chns],[1 2],500,FigureProperties,NS_GlobalConstants,[144],[]);
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks4Rotated(T1,[Chns],[1 2],500,FigureProperties,NS_GlobalConstants,[144],[]);

h=gcf;
FullName=['C:\home\pawel\nauka\MEA2010\prezentacja\figures_pok107\' 'EIs_culture_rotated'];            
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 7]);
set(h,'PaperPosition',[0 0 10 7]);  
print(h, '-dtiff', '-r120', FullName);