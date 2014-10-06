NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ArrayID=500;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
Radius=1;

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-200 100],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
DataPath='D:\analysis\retina\2009-11-27-0\data001';
WritePathFigs='D:\analysis\retina\2009-11-27-0\data001\figures';

PatternNumber=24;
MovieNumber=47;
N=10;

[DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);    
DataTraces=DataTraces0(1:100,Channels,1:50);
[Traces,Artifact,c]=NS512_SubtractLocalArtifact(DataTraces,N);
y=NS512_PlotClustersOfSignaturesOnArrayLayout(Artifact,Channels,1,ArrayID,FigureProperties,NS_GlobalConstants);        
%EventsPerChannel=sum(Events)
%a=max(EventsPerChannel)
%[p,index]=find(EventsPerChannel==a)
%WaveformTypes=Events(:,index(1));
WaveformTypes=ones(1,100);
FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-200 100],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
y=NS512_PlotClustersOfSignaturesOnArrayLayout(Traces,Channels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);        
Events=NS512_DetectSpikes(Traces,25,3);