ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;

cd C:\praca\data\2008-02-04-0;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

FileName='005';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.neurons');
idList = neuronFile.getIDList();
NeuronID=886;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings=spikeTimes(1,1:400);
Channels=[53:56 59:61];
TimeStart=-10;
NumberOfSamples=50;
Offsets=ones(length(Channels))*(-370);

Traces=NS_ReadManyTracesFromRaw(FileName,Channels,Timings,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);

figure(3)
EI=NS_CalculateEI(Traces);
WaveformTypes=zeros(length(Channels));
ArrayID=1;
FigureNumber=1;
AmplitudeRange=[-200 200];
FigureProperties=struct('FigureNumber',4,'TimeRange',[-10 40],'AmplitudeRange',AmplitudeRange,'FontSize',13,'Colors',['k' 'r' 'b' 'y']);

y=NS_PlotSignatureOnArrayLayout(EI,Channels,WaveformTypes,ArrayID,FigureNumber,FigureProperties,NS_GlobalConstants);