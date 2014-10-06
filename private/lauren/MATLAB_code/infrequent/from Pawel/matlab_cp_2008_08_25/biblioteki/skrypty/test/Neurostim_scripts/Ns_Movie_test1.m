ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges,'ArrayID',1);

RecElectrodes=[1:64];
BadElectrodes=[9 25 28 31 33 37 41 57 64];
RecElectrodes=NS_RemoveBadChannels(RecElectrodes,BadElectrodes);

FileName='007';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-06-02-0\VisionNeurons\data007\data007\data007000\data007000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-06-02-0\VisionNeurons\data007\data007\data007000\data007000.neurons');
%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-06-02-0\VisionNeurons\data000\data000\data000000\data000000.params');
%neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-06-02-0\VisionNeurons\data000\data000\data000000\data000000.neurons');
idList = neuronFile.getIDList();
NeuronID=153;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-14;
TimeRange=[-0 50];
signal=NS_AverageTraces(FileName,Timings1,RecElectrodes,TimeRange,NS_GlobalConstants);
signal=(signal+370);
%T(1,:,:)=signal';

FigureProperties=struct('FigureNumber',63,'TimeRange',TimeRange,'AmplitudeRange',[-400 400],'FontSize',12,'Colors',['g' 'b' 'r' 'y']);
M=NS_SaveMovieFromSignature(signal',RecElectrodes,1,FigureProperties,NS_GlobalConstants);
%movie(M);