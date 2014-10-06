ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges,'ArrayID',1);

BadElectrodes=[9 25 28 31 33 37 41 57 64];
StimElectrodes=21;
RecElectrodes=[22 15 4 1 21 17 12 6 3 18 14 11 7 5 8 10 13 16 19];
RecElectrodes=[1:8 10:24 27];
RecElectrodes=NS_RemoveBadChannels(RecElectrodes,BadElectrodes);
%Movies=[18:25];
MovieNumber1=24;
MovieNumber2=25;

ReadPath='D:\2008-06-02-0';

WritePath='E:\analysis\2008-06-02-0\2008-06-09\data005';
FileName='005';
%Types=NS_ClusterAndPrint(StimElectrodes,RecElectrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);

WritePath='E:\analysis\2008-06-02-0\data006';
FileName='006';
%Types=NS_ClusterAndPrint(Electrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);

%WritePath='E:\analysis\2008-06-02-0\2008-06-09\data010';
WritePath='D:\analysis\2008-06-02-0\2008-06-13\data010CleanClusters';
FileName='010';
DiffEI=NS_PrintDiffEIs(StimElectrodes,RecElectrodes,MovieNumber1,MovieNumber2,BadElectrodes,ReadPath,FileName,WritePath);
%T(2,:,:)=DiffEI';
sdiff=size(DiffEI)

FileName='007';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-06-02-0\VisionNeurons\data007\data007\data007000\data007000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-06-02-0\VisionNeurons\data007\data007\data007000\data007000.neurons');
idList = neuronFile.getIDList();
NeuronID=153;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-14;
TimeRange=[-0 50];
signal=NS_AverageTraces(FileName,Timings1,RecElectrodes,TimeRange,NS_GlobalConstants);
signal=(signal+370);
%T(1,:,:)=signal';
data=zeros(2,24,51);
data(1,:,:)=DiffEI(1,:,:);
data(2,:,:)=signal';
FigureProperties=struct('FigureNumber',53,'TimeRange',TimeRange,'AmplitudeRange',[-400 400],'FontSize',12,'Colors',['g' 'b' 'r' 'y']);
y=NS_PlotDifferentTracesOnArrayLayout(data,RecElectrodes,0,1,FigureProperties,NS_GlobalConstants);

FigureProperties=struct('FigureNumber',63,'TimeRange',TimeRange,'AmplitudeRange',[-400 400],'FontSize',12,'Colors',['g' 'b' 'r' 'y']);
M=NS_SaveMovieFromSignature(signal',RecElectrodes,1,FigureProperties,NS_GlobalConstants);
%movie(M);