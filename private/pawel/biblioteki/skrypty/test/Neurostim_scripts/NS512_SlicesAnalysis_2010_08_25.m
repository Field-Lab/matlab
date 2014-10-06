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
DataPath='G:\analysis\2010-07-29-0\data';
ArtifactDataPath='G:\analysis\2010-07-29-0\artifacts';
%FilePath='I:\analysis\2010-07-29-0\data\ClusterFile_002a';
Pattern=197;
%Pattern=126;
Movie=153;
%Movie=113;
RecElectrodes=[81 82 102 122 166 193 209 242 296];
%RecElectrodes=[72 81 102 163 189 211];
figure(Movie);
clf;
hold on;
for i=1:length(RecElectrodes)
    ClusterFileName=['G:\analysis\2010-07-29-0\data\ClusterFile_002_el' num2str(RecElectrodes(i))];
    ClusterIndex0=NS_ReadClusterFile(ClusterFileName,Movie,Pattern,200); 
    ClusterIndex=(sign(ClusterIndex0-1.5)+1)./2;
    %subplot(2,3,i);
    figure(Movie);
    %h=plot(-ClusterIndex(1:100)*0.5+i+0.25,'bd-');    
    h=plot(ClusterIndex(1:100)*0.5-i+10-0.25,'bd-');    
    set(h,'MarkerSize',5);
    h=gca;
    set(h,'YTick',[1:length(RecElectrodes)]);    
    set(h,'YTickLabel',RecElectrodes([length(RecElectrodes):-1:1]));
    %set(h,'YDir','reverse')    
    set(h,'FontSize',14);
    h=xlabel('Trial number');
    set(h,'FontSize',14);
    h=ylabel('Electrode number & stimulation efficacy');    
    
    s=find(ClusterIndex==1);
    %plot(s,ones(1,length(s))*(-2*i)+1,'rd');   
    [DataTraces0,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,1,Pattern,Movie,200,0);
    figure(10);
    subplot(3,3,i);
    s0=DataTraces0(1:100,RecElectrodes(i),:);
    s=reshape(s0,100,140);
    h=plot([0:139]/20,reshape(s0,100,140)');
    FaultSpikes=find(ClusterIndex0(1:100)==3);
    spikes=find(ClusterIndex(1:100)==1);
    arts=find(ClusterIndex(1:100)==0);
    set(h(arts),'Color','k');
    set(h(spikes),'Color','r');
    set(h(FaultSpikes),'Color','g');
    axis([0 7 -80 40]);
    grid on;
    xlabel('Time [ms]');
    ylabel('Amplitude [mV]');
    h=gca;
    set(h,'FontSize',14);
    set(h,'XTick',0:1:7);
    h=text(0.5,-70,['el. ' num2str((RecElectrodes(i)))]);
    set(h,'FontSize',14);
    
    figure(11);
    subplot(3,3,i);
    s0=DataTraces0(1:100,RecElectrodes(i),:);
    s=reshape(s0,100,140);
    h=plot(reshape(s0,100,140)');    
    set(h(33:100),'Color','r');
    set(h(1:32),'Color','k');
    axis([0 140 -80 40]);
    grid on;
end
%break;
figure(Movie);
FullName=['G:\analysis\2010-07-29-0\figures\' 'p' num2str(Pattern) '_m' num2str(Movie) '9_neurons_successes'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 4]);
set(h,'PaperPosition',[0 0 10 4]); 
print(h, '-dtiff', '-r240', FullName);

figure(10);
FullName=['G:\analysis\2010-07-29-0\figures\' 'p' num2str(Pattern) '_m' num2str(Movie) '9_neurons_Traces'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 7]);
set(h,'PaperPosition',[0 0 10 7]); 
print(h, '-dtiff', '-r240', FullName);

break;

ClusterIndex0=NS_ReadClusterFile(FilePath,Movie,Pattern,200);
ClusterIndex=ClusterIndex0(1:100);

DataPath='I:\analysis\2010-07-29-0\data';
ArtifactDataPath='I:\analysis\2010-07-29-0\artifacts';
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
Chns=find(Indexes(:,1)<27 & Indexes(:,1)>9);
%Chns=[65:320];

%m=mean(Spike,3);
%for i=Chns
%    Spike(1,i)=Spike(1,i)-m(i);    
%end

DataTracesFinal=DataTraces;
for i=1:100
    DataTracesFinal(i,:,:)=DataTraces(i,:,:); %-Artifact(1,:,:);
end

WritePathFigs='I:\analysis\2010-07-29-0\figures';
FigureProperties=struct('FigureNumber',6,'TimeRange',[0 140],'AmplitudeRange',[-80 40],'FontSize',14,'Colors',['b' 'r' 'k' 'g'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',2,'YLabel','Signal [\muV]');
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3(Spike(1,Chns,:),Chns,[1],500,FigureProperties,NS_GlobalConstants,[RecElectrode],[Pattern]);
FullName=[WritePathFigs '\' 'p' num2str(Pattern) '_m' num2str(Movie) '_el' num2str(RecElectrode) '_EI'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 13]);
set(h,'PaperPosition',[0 0 10 13]); 
print(h, '-dtiff', '-r240', FullName);

FigureProperties=struct('FigureNumber',7,'TimeRange',[0 140],'AmplitudeRange',[-80 40],'FontSize',14,'Colors',['b' 'k' 'r' 'w'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',0.5,'YLabel','Signal [\muV]');
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks3(DataTracesFinal(:,Chns,:),Chns,ClusterIndex,500,FigureProperties,NS_GlobalConstants,[RecElectrode],[Pattern]);
FullName=[WritePathFigs '\' 'p' num2str(Pattern) '_m' num2str(Movie) '_el' num2str(RecElectrode) '_Traces'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 13]);
set(h,'PaperPosition',[0 0 10 13]); 
print(h, '-dtiff', '-r240', FullName);