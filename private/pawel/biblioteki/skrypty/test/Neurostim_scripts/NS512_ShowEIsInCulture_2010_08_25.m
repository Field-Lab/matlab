NS_GlobalConstants=NS_GenerateGlobalConstants(512);

% * * * * * Spontaneous EI * * * * 
%FileName='E:\pawel\data\cultures\2009-11-20-0\data002'
%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data002\data002.params');
%neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data002\data002.neurons');
%idList = neuronFile.getIDList();
%NeuronID=6647;
%spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
%N=500;
%L=100;
%Timings1=spikeTimes(1,1:min(N,length(spikeTimes)))-17;
%Data=zeros(513,L);
%rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(FileName); 
%for i=1:length(Timings1)
%    RawData=double(rawFile.getData(Timings1(i),L)');    
%    Data=Data+RawData;
%end
%m=mean(Data')';
%for i=1:513
%    Data(i,:)=Data(i,:)-m(i);
%end

%T(2,:,:)=Data(2:513,:)/N;

%  * * * * * Elicited EI  * * * * 
Pattern=197;
Movie=153;
RecElectrode=193;

ClusterFileName=['G:\analysis\2010-07-29-0\data\ClusterFile_002_el' num2str(RecElectrode)];
ClusterIndex0=NS_ReadClusterFile(ClusterFileName,Movie,Pattern,200);
ClusterIndex=ClusterIndex0(1:100);
ClusterIndex(1:34)=1;
ClusterIndex(35:100)=2;

DataPath='G:\analysis\2010-07-29-0\data';
ArtifactDataPath='G:\analysis\2010-07-29-0\artifacts';
%(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
[DataTraces0,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,1,Pattern,Movie,200,0);
DataTraces=DataTraces0(1:100,:,:);

artifacts=find(ClusterIndex==1);
%artifacts=find(artifacts0<103);
spikes=find(ClusterIndex==2);
clear Artifact;
Artifact=mean(DataTraces(artifacts,:,:));
S=mean(DataTraces(spikes,:,:));
Spike=S-Artifact;

Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);
Chns=find(Indexes(:,1)<17 & Indexes(:,1)>0);
%Chns=[65:320];

%m=mean(Spike,3);
%for i=Chns
%    Spike(1,i)=Spike(1,i)-m(i);    
%end

DataTracesFinal=DataTraces;
for i=1:100
    DataTracesFinal(i,:,:)=DataTraces(i,:,:); %-Artifact(1,:,:);
end

WritePathFigs='G:\analysis\2010-07-29-0\figures2';
FigureProperties=struct('FigureNumber',6,'TimeRange',[0 140],'AmplitudeRange',[-80 40],'FontSize',14,'Colors',['b' 'r' 'k' 'g'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',2,'YLabel','Signal [\muV]');
%y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3(Spike(1,Chns,:),Chns,[1],500,FigureProperties,NS_GlobalConstants,[RecElectrode],[Pattern]);
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3(Spike(1,Chns,:),Chns,[1],500,FigureProperties,NS_GlobalConstants,[81 82 102 122 166 193 209 242 296],[Pattern]);
FullName=[WritePathFigs '\' 'p' num2str(Pattern) '_m' num2str(Movie) 'Switch_EI'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 13]);
set(h,'PaperPosition',[0 0 10 13]); 
print(h, '-dtiff', '-r240', FullName);

FigureProperties=struct('FigureNumber',7,'TimeRange',[0 140],'AmplitudeRange',[-80 40],'FontSize',14,'Colors',['b' 'k' 'r' 'w'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',0.5,'YLabel','Signal [\muV]');
%y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3(DataTracesFinal(:,Chns,:),Chns,ClusterIndex,500,FigureProperties,NS_GlobalConstants,[RecElectrode],[Pattern]);
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3(DataTracesFinal(:,Chns,:),Chns,ClusterIndex,500,FigureProperties,NS_GlobalConstants,[81 82 102 122 166 193 209 242 296],[Pattern]);
%FullName=[WritePathFigs '\' 'p' num2str(Pattern) '_m' num2str(Movie) '_el' num2str(RecElectrode) '_Traces'];            
FullName=[WritePathFigs '\' 'p' num2str(Pattern) '_m' num2str(Movie) '_Switch_Traces'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 13]);
set(h,'PaperPosition',[0 0 10 13]); 
print(h, '-dtiff', '-r240', FullName);