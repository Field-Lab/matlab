ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges,'ArrayID',1);

cd H:\2008-02-19\2008-02-19-1;
path=pwd;
FileName='013';
Channel=11;
MovieNumber=4;
TimeRange=[-15 125];
Limit=100;
AmplitudeRange=[-200 200];
%AdditionalChannels=[Channel+1];
Channels=[13 14 16 18 19];
Channels=[6 7 10 11 12 13 14];
T=zeros(1,length(Channels),TimeRange(2)-TimeRange(1)+1);

%[CI,MI]=NS_Report(FileName,1,NS_GlobalConstants);

%1) Plotting the stimulation waveform
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);
[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);
figure(1);
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,Channel,20,'b-',NS_GlobalConstants);
set(PlotPointer,'LineWidth',2);

%2) Plotting many traces on figure - only artifacts vs artifacts+spikes
FigureNumber=23;
FigureProperties=struct('FigureNumber',FigureNumber,'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',13,'Colors',['k' 'r' 'b' 'y']);
Delay=0;
WaveformTypes=ones(1,length(Timings));
%WaveformTypes([47 48 61 78 79 80 88 89 90 91 92 93])=0; %channel 16, movie 17
%WaveformTypes([13 14 17 25 26 27 32 33 34 39 40 44 47 49 52 53 54 57 59 60 65 66 67 72 73 77 78 86 87 88 91 92 93 94 95 98 99])=0; %channel 11, movie 17
WaveformTypes([13 19 27 31 32 47 49 52 60 64 69 70 74 75 80 90 91 94 95 96])=0; %channel 11, movie 4
%WaveformTypes([40 53:57 76:78 92:95])=0; %channel 16, movie 15
%WaveformTypes(a)=0;
%clear a;
Channel=13;
l1=1;
l2=100;
Timings=Timings(1,l1:l2);
WaveformTypes=WaveformTypes(1,l1:l2);
%WaveformTypes=WaveformTypes(1,50:100);
y=NS_PlotManyTracesOnFigure(FileName,Timings,WaveformTypes,Channel,2,[1:10],FigureProperties,Delay,NS_GlobalConstants);
%artifact16=NS_AverageTraces(FileName,Timings(a),Channels,TimeRange,NS_GlobalConstants)
%3) Finding the EIs based on elicited spikes
a=find(WaveformTypes==1);
[RAWtraces,signal1]=NS_AverageTraces([path '\data' FileName],Timings(a),Channels,TimeRange,NS_GlobalConstants);
a=find(WaveformTypes==0);
[RAWtraces,signal2]=NS_AverageTraces([path '\data' FileName],Timings(a),Channels,TimeRange,NS_GlobalConstants);
T(1,:,:)=signal1'-signal2';

%4) Finding EI for spontaneous activity
FileName='001';
%cd H:\2008-02-19\2008-02-19-0;
full_path='H:\2008-02-19\2008-02-19-0\data001\data001000.bin';
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\retina\2008-02-19-0\2008-02-19-0\data001\data001.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\retina\2008-02-19-0\2008-02-19-0\data001\data001.neurons');
idList = neuronFile.getIDList();
%NeuronID=226;
NeuronID=152;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-10;
L=141;
N=200;
spikes2=zeros(N,numel(Channels),L); %spike,channel,sample
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-20,L)';
    d1=d0(Channels+1,:);
    spikes2(i,:,:)=d1;
end




WaveformTypes=ones(1,length(spikeTimes));
%FigureProperties=struct('FigureNumber',15,'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',13,'Colors',['k' 'r' 'b' 'y']);
%Channel=11;
%y=NS_PlotManyTracesOnFigure(FileName,Timings1,WaveformTypes,Channel,2,[1:10],FigureProperties,0,NS_GlobalConstants);
FileName='H:\2008-02-19\2008-02-19-0\data001';
signal=NS_AverageTraces(FileName,Timings1,Channels,TimeRange,NS_GlobalConstants);
T(2,:,:)=signal';

%cd E:\2008-02-19-1;
%6) Plotting both EIs
FigureProperties=struct('FigureNumber',3,'TimeRange',TimeRange,'AmplitudeRange',[-150 50],'FontSize',11,'Colors',['g' 'b' 'r' 'y']);
OffsetCancellation=2;
OffsetSamples=[1:-TimeRange(1)-1];
OffsetCancellation=2;
OffsetSamples=[1:-TimeRange(1)-4];
y=NS_PlotDifferentTracesOnArrayLayout(T,Channels,OffsetCancellation,OffsetSamples,FigureProperties,NS_GlobalConstants);
