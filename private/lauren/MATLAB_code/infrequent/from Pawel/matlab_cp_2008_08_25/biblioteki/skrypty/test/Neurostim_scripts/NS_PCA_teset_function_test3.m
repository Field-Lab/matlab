BadElectrodes=[9 25 28 31 33 37 41 57 64];

StimElectrodes=16;
RecElectrodes=[16 13];

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
%RecElectrodes = electrodeMap.getAdjacentsTo(RecPrimaryElectrode,1);
RecElectrodes=NS_RemoveBadChannels(RecElectrodes,BadElectrodes);

Movies=[22:29];

ReadPath='C:\praca\data\2008-03-21-0';
WritePath='C:\praca\analiza\MEA2008\prezentacja\obrazki2';
%ReadPath='E:\2008-06-02-0';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 60],'AmplitudeRange',[-250 250],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
FileName='003';

Movie=22;
NumberOfClusters=1;
ClustersForDiff(1,:)=[1,1];
FigureProperties.Colors=['g' 'k' 'b' 'm' 'k'];
%Types=NS_ClusterAndPrint_meeting(StimElectrodes,RecElectrodes,Movie,BadElectrodes,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);

%break;

FigureProperties.Colors=['g' 'k' 'r' 'm' 'k']
%FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-250 250],'FontSize',14,'Colors',['g' 'r' 'b' 'm' 'k']);
%FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 60],'AmplitudeRange',[-250 250],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
Movie=23;
NumberOfClusters=2;
ClustersForDiff=[2 1];
TracesNumbers=[1:3 5:63 65:98];
%Types=NS_ClusterAndPrint_meeting(StimElectrodes,RecElectrodes,Movie,BadElectrodes,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);

Movie=27;
NumberOfClusters=2;
%ClustersForDiff(1,:)=[1,1];
FigureProperties.Colors=['r' 'b' 'r' 'm' 'k'];
%Types=NS_ClusterAndPrint_meeting(StimElectrodes,RecElectrodes,Movie,BadElectrodes,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);

%break;

Movie=27;
NumberOfClusters=3;
ClustersForDiff(1,:)=[2,1];
%TracesNumbers=[1:3 5:63 65:98];
TracesNumbers=[1:73 75:80 82:98];
TracesNumbers=[1:68 70:82 84:98];
%Types=NS_ClusterAndPrint_meeting(StimElectrodes,RecElectrodes,Movie,BadElectrodes,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);

%break;

Movie=28;
NumberOfClusters=4;
ClustersForDiff(1,:)=[1,3];
ClustersForDiff(2,:)=[4,3];
ClustersForDiff(3,:)=[2,4];
ClustersForDiff(4,:)=[2,1];
%ClustersForDiff=[1 2];
%TracesNumbers=[1:3 5:63 65:98];
TracesNumbers=[1:73 75:80 82:98];
FigureProperties.Colors=['b' 'r' 'b' 'm' 'k'];
Types=NS_ClusterAndPrint_meeting(StimElectrodes,RecElectrodes,Movie,BadElectrodes,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);

