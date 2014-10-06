ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;

cd C:\praca\data\2008-02-04-0;
FileName='005';
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.neurons');

idList = neuronFile.getIDList();
NeuronID=886;
spikeTimes = neuronFile.getSpikeTimes(886)';

Channel=56;

TimeRange=[-20 60];

Limit=100;
AmplitudeRange=[-100 100];
FigureProperties=struct('TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',11);

Channels=[53:56 58:62];
%Channels=[1:64];

Timings=spikeTimes(1,1:min(1000,length(spikeTimes)));
%NS_PlotManyTracesOnArrayLayout(FileName,Channels,Timings,1,24,FigureProperties,NS_GlobalConstants);

signal=NS_AverageTraces(FileName,Timings,Channels,TimeRange,NS_GlobalConstants);
%FigureProperties=struct('TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',11);
FigureProperties=struct('FigureNumber',3,'TimeRange',TimeRange,'AmplitudeRange',[-100 100],'FontSize',11,'Colors',['g' 'b' 'r' 'y']);
OffsetCancellation=2;
OffsetSamples=[1:4];
WaveformTypes=ones(1,length(Channels));
y=NS_PlotTracesOnArrayLayout(signal',Channels,WaveformTypes,OffsetCancellation,OffsetSamples,1,32,FigureProperties,NS_GlobalConstants);
