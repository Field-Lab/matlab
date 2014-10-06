ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges,'ArrayID',1);

cd C:\praca\data\2008-02-04-0;
FileName='004';
Channel=56;
MovieNumber=87;
TimeRange=[-10 70];
Limit=100;
AmplitudeRange=[-700 -500];
AdditionalChannels=[Channel+1];
Channels=[53:56 59:61];
T=zeros(2,length(Channels),TimeRange(2)-TimeRange(1)+1);

%[CI,MI]=NS_Report(FileName,1,NS_GlobalConstants);

%1) Plotting the stimulation waveform
figure(1);
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);
[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);

X=1000;
z=zeros(5,2);
Amp=30;
Pulse(1,:)=[0 0 0 -Amp -Amp Amp Amp Amp 0];
Pulse3=[Pulse z];
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse3,Status,Channel,0,X,'b-',NS_GlobalConstants);
set(PlotPointer,'LineWidth',2);
hold on
Pulse(1,:)=[Amp/2 Amp/2 Amp/2 -Amp -Amp Amp/2 Amp/2 Amp/2 0];
Pulse1=[z Pulse z];
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse1,Status,Channel,13,X,'b-',NS_GlobalConstants);
set(PlotPointer,'LineWidth',2);
Pulse(1,:)=[Amp*2/3 Amp*2/3 Amp*2/3 -Amp -Amp Amp/3 Amp/3 Amp/3 0];
Pulse2=[z Pulse z];
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse2,Status,Channel,29,X,'b-',NS_GlobalConstants);
%set(PlotPointer,'LineWidth',2);

grid off;
axis([-0.25 2.25 -1 0.5]);

%2) Plotting many traces on figure - only artifacts vs artifacts+spikes
FigureNumber=25;
FigureProperties=struct('FigureNumber',FigureNumber,'TimeRange',TimeRange,'AmplitudeRange',[-100 350],'FontSize',14,'Colors',['k' 'r' 'b' 'y']);
Delay=0;
WaveformTypes=ones(1,length(Timings));
WaveformTypes([5 8 14 26 45 62 63 69 72 73 79 82 95 97 103 106 116 124 130 142])=0;
%(FileName,Timings,WaveformTypes,Channel,OffsetCancellation,OffsetSamples,FigureProperties,Delay,NS_GlobalConstants)
y=NS_PlotManyTracesOnFigure(FileName,Timings,WaveformTypes,Channel,3,[-365],FigureProperties,Delay,NS_GlobalConstants);

%3) Finding the EIs based on elicited spikes
a=find(WaveformTypes==1);
signal1=NS_AverageTraces(FileName,Timings(a),Channels,TimeRange,NS_GlobalConstants);
a=find(WaveformTypes==0);
signal2=NS_AverageTraces(FileName,Timings(a),Channels,TimeRange,NS_GlobalConstants);
T(2,:,:)=signal1'-signal2';

%4) Finding EI for spontaneous activity
FileName='005';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.neurons');
idList = neuronFile.getIDList();
NeuronID=886;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-14;
signal=NS_AverageTraces(FileName,Timings1,Channels,TimeRange,NS_GlobalConstants);
T(1,:,:)=signal';

%6) Plotting both EIs
FigureProperties=struct('FigureNumber',3,'TimeRange',TimeRange,'AmplitudeRange',[-100 100],'FontSize',13,'Colors',['g' 'b' 'r' 'y']);
OffsetCancellation=2;
OffsetSamples=[1:-TimeRange(1)-1];
OffsetCancellation=2;
OffsetSamples=[1:-TimeRange(1)-4];
y=NS_PlotDifferentTracesOnArrayLayout(T,Channels,OffsetCancellation,OffsetSamples,FigureProperties,NS_GlobalConstants);

% new
data=rawFile.getData(Timings(j)+TimeRange(1),TimeRange(2)-TimeRange(1)+1);
signal=double(data(:,Channels(i)+1)); %first channel is a TTL channel  
h=NS_PlotManyTracesOnArrayLayoutNew(Waveforms,Channels,WaveformTypes,ArrayID,FigureNumber,FigureProperties,NS_GlobalConstants);