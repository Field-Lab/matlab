NS_GlobalConstants=NS_GenerateGlobalConstants(512);
FigureProperties=struct('FigureNumber',4,'Subplot',[2 3 3],'TimeRange',[0 60],'AmplitudeRange',[-20 10],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','signal [mV]');

DataPath='D:\analysis\2009-11-20-0\data003c';
ClusterFilePath='D:\analysis\2009-11-20-0\data003b\ClusterFile_003';
[EI1,Channels1]=NS512_EI_FromClusteredData(DataPath,DataPath,0,ClusterFilePath,48,37,100,0);
DataPath='D:\analysis\2009-11-20-0\data003d';
ClusterFilePath='D:\analysis\2009-11-20-0\data003b\ClusterFile_003';
[EI2,Channels2]=NS512_EI_FromClusteredData(DataPath,DataPath,0,ClusterFilePath,48,37,100,0);
%GoodChannels=[1:32];
EI=[EI1 EI2];
EI=EI(1,10:89,:);
Channels=[Channels1' Channels2']';
Channels=Channels(10:89);
[StChannels,Amplitudes]=NS_StimulatedChannels(DataPath,43,54,[1:512],NS_GlobalConstants)
%[EI,Channels]=NS512_EI_FromClusteredData(DataPath,DataPath,0,ClusterFilePath,48,37,100,0);

y=NS_PlotClustersOfSignaturesOnArrayLayout(EI,Channels,1,500,FigureProperties,NS_GlobalConstants);