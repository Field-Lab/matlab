%ElectrodeNumbers=[20 22 24];

%ElectrodeNumbers=[42 43 45];
%ElectrodeNumbers=[14]
BadElectrodes=[9 25 28 31 33 37 41 57 64];
%Electrodes=[];
%for i=ElectrodeNumbers
%    active=1;
%    for j=BadElectrodes
%        if i==j
%            active=0;
%        end
%    end
%    if active==1
%        Electrodes=[Electrodes i];
%    end
%end
%Electrodes;
%Movies=[27:33];

StimElectrodes=[21 17 12 6 3 2 22 15 18 14 4 11];
RecPrimaryElectrode=17;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
RecElectrodes = electrodeMap.getAdjacentsTo(RecPrimaryElectrode,1);
RecElectrodes=NS_RemoveBadChannels(RecElectrodes,BadElectrodes);
RecElectrodes=[22 15 4 1 21 17 12 6 3 18 14 11 7];
StimElectrodes=21;
%StimElectrodes=[6 7];
Movies=[18:25];

%ReadPath='E:\2008-03-21-0';
ReadPath='E:\2008-06-02-0';

%WritePath='E:\analysis\2008-03-21-0\forMEA';
WritePath='E:\analysis\2008-06-02-0\2008-06-09\data005';
FileName='005';
%Types=NS_ClusterAndPrint(StimElectrodes,RecElectrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);

%WritePath='E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_bi';
WritePath='E:\analysis\2008-06-02-0\data006';
FileName='006';
%Types=NS_ClusterAndPrint(Electrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);

%WritePath='E:\analysis\2008-03-21-0\2008_05_05\report\subtr\50us_tri';
WritePath='E:\analysis\2008-06-02-0\2008-06-11\data010CleanClusters';
FileName='010';
Types=NS_ClusterAndPrint(StimElectrodes,RecElectrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);

%WritePath='E:\analysis\2008-03-21-0\2008_05_05\report\subtr\50us_bi';
FileName='006';
%Types=NS_ClusterAndPrint(Electrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);

%WritePath='E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri_2nd_run';
FileName='007';
%Types=NS_ClusterAndPrint(Electrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);