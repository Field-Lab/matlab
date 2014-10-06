BadElectrodes=[9 25 28 31 33 37 41 57 64];

StimElectrodes=[1:64];
StimElectrodes=NS_RemoveBadChannels(StimElectrodes,BadElectrodes)
RecElectrodes=[0];

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
%RecElectrodes = electrodeMap.getAdjacentsTo(RecPrimaryElectrode,1);
RecElectrodes=NS_RemoveBadChannels(RecElectrodes,BadElectrodes)

Movies=[1:25];

ReadPath='E:\2008-06-02-0'; %'C:\praca\data\2008-03-21-0';
WritePath='E:\analysis\2008-06-02-0\2008-07-21'; %'C:\praca\analiza\MEA2008\prezentacja\obrazki2';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 60],'AmplitudeRange',[-250 250],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
FileName='005';

Types=NS_PrintOverlaidResponsesManyAmplitudes(StimElectrodes,RecElectrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath,FigureProperties);
FigureProperties.Colors=['g' 'k' 'b' 'm' 'k'];
%Types=NS_ClusterAndPrint_meeting(StimElectrodes,RecElectrodes,Movie,BadElectrodes,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);
