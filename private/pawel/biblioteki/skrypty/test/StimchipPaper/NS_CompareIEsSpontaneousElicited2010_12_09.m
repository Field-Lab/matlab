ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges,'ArrayID',1);

FontName='Arial';

cd F:\2008-02-04-0;
FileName='004';
Channel=56;
MovieNumber=87;
TimeRange=[-10 50];
Limit=100;
AmplitudeRange=[-700 -500];
AdditionalChannels=[Channel+1];
Channels=[53:56 59:61];
T=zeros(2,length(Channels),TimeRange(2)-TimeRange(1)+1);

FontSize=20;
BoxLineWidth=2;

%1) Plotting the stimulation waveform
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);
Timings=Timings+8;

%2) Plotting many traces on figure - only artifacts vs artifacts+spikes
FigureNumber=225;
FigureProperties=struct('FigureNumber',FigureNumber,'TimeRange',TimeRange,'AmplitudeRange',[-100 400],'FontSize',FontSize,'FontName',FontName,'Colors',['k' 'r' 'b' 'y']);
Delay=0;
WaveformTypes=ones(1,length(Timings));
WaveformTypes([5 8 14 26 45 62 63 69 72 73 79 82 95 97 103 106 116 124 130 142])=0;
%(FileName,Timings,WaveformTypes,Channel,OffsetCancellation,OffsetSamples,FigureProperties,Delay,NS_GlobalConstants)
y=NS_PlotManyTracesOnFigure(FileName,Timings-1,WaveformTypes,Channel,3,[-365],FigureProperties,Delay,NS_GlobalConstants);
hold on;

%2b) stimulation waveform
x=[-1 1 1 3 3 5 5 7 7 9];
x=x*0.05-0.35;
y=[0 0 2 2 -3 -3 1 1 0 0];
Skalowanie=9;
h=plot(x,y*Skalowanie+335);
set(h,'LineWidth',2);h=text(0.15,335,'Stim');
set(h,'FontSize',FontSize);
set(h,'FontName',FontName);

%2c) Record switch timing
x=[0 1 1 9 9 10];
x=x*0.05-0.4;
y=[1 1 0 0 1 1];
h=plot(x,y*16+370);
set(h,'LineWidth',2);
h=text(0.15,378,'Rec');
set(h,'FontSize',FontSize);
set(h,'FontName',FontName);

%2d) set the scale and plot to file
h=gca;
set(h,'LineWidth',BoxLineWidth);
set(h,'XLim',[-0.5 2.5]);
hc = figure(FigureNumber);
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[7 5]);
set(hc,'PaperPosition',[0 0 7 5]);
print(hc, '-dtiff', '-r400', 'C:\home\pawel\nauka\Stimchip_paper\obrazki\Traces');

%3) Plot magnified version of the previous plot
FigureNumberZoom=226;
FigureProperties=struct('FigureNumber',FigureNumberZoom,'TimeRange',[-8 6],'AmplitudeRange',[-100 400],'FontSize',FontSize,'FontName',FontName,'Colors',['k' 'r' 'b' 'y']);
y=NS_PlotManyTracesOnFigure(FileName,Timings-1,WaveformTypes,Channel,3,[-365],FigureProperties,Delay,NS_GlobalConstants);
h=gca;
set(h,'LineWidth',BoxLineWidth);
set(h,'XLim',[-0.4 0.3]);
hold on;

%3a) stim waveform
x=[-1 1 1 3 3 5 5 7 7 9];
x=x*0.05-0.35;
y=[0 0 2 2 -3 -3 1 1 0 0];
Skalowanie=9;
h=plot(x,y*Skalowanie+335);
set(h,'LineWidth',2);
h=text(0.115,335,'Stim');
set(h,'FontSize',FontSize);
set(h,'FontName',FontName);

%3b) Record switch timing
x=[0 1 1 9 9 10];
x=x*0.05-0.4;
y=[1 1 0 0 1 1];
h=plot(x,y*18+370);
set(h,'LineWidth',2);
h=text(0.115,378,'Rec');
set(h,'FontSize',FontSize);
set(h,'FontName',FontName);

%3c) Set scale and print
h=gca;
set(h,'LineWidth',BoxLineWidth);
set(h,'XTick',[-0.4:0.1:0.3]);
hc = figure(FigureNumberZoom);
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[7 5]);
set(hc,'PaperPosition',[0 0 7 5]);
print(hc, '-dtiff', '-r400', 'C:\home\pawel\nauka\Stimchip_paper\obrazki\TracesZoom');

%4) Finding the EIs based on elicited spikes and plotting it in fig. 133
full_path='F:\2008-02-04-0\data004';
a=find(WaveformTypes==1);
size(a)
[RAWtraces,signal1]=NS_AverageTraces(full_path,Timings(a)-1,Channels,TimeRange,NS_GlobalConstants);
a=find(WaveformTypes==0);
size(a)
[RAWtraces,signal2]=NS_AverageTraces(full_path,Timings(a)-1,Channels,TimeRange,NS_GlobalConstants);
clear T;
T(1,:,:)=signal1'-signal2';
T0(1,:,:)=signal1'-signal2';
FigureProperties=struct('FigureNumber',103,'TimeRange',TimeRange,'AmplitudeRange',[-150 100],'FontSize',FontSize,'FontName',FontName,'Colors',['y' 'r' 'k' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-100 -50 0 50 100]);
OffsetCancellation=2;
OffsetSamples=[1:-TimeRange(1)-4];
y=NS_PlotDifferentTracesOnArrayLayout(T/0.84,Channels,OffsetCancellation,OffsetSamples,FigureProperties,NS_GlobalConstants);

T3=T(1,1,:);
y=NS_PlotDifferentTracesOnArrayLayout(T3/0/84,[56],OffsetCancellation,OffsetSamples,FigureProperties,NS_GlobalConstants);

%5) Finding EI for spontaneous activity
FileName='005';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\retina\2008-02-04-0\data005\data005000\data005000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\retina\2008-02-04-0\data005\data005000\data005000.neurons');
idList = neuronFile.getIDList();
NeuronID=886;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-6;
full_path='F:\2008-02-04-0\data005';
[RAWtraces,signal]=NS_AverageTraces(full_path,Timings1-1,Channels,TimeRange,NS_GlobalConstants);
T(2,:,:)=signal';

%6) Plot stimulation EI for stimulated channel

spike=T(1,4,:);
Spike=reshape(spike,1,61);
FigureNumberSpike=301;
figure(FigureNumberSpike);
t=[-0.5:0.05:2.5];
h=plot(t,Spike/0.84+2);
grid on;
set(h,'Color','r')
set(h,'LineWidth',2.5);
h=gca;
set(h,'LineWidth',BoxLineWidth);
set(h,'XLim',[-0.5 2.5]);
set(h,'YLim',[-125 50]);
set(h,'FontSize',FontSize);
set(h,'FontName',FontName);
set(h,'XTick',[-0.5:0.5:2.5]);
set(h,'YTick',[-125:25:50]);
xlabel('Time [ms]');
ylabel('Signal [\muV]');
hold on;

x=[-1 1 1 3 3 5 5 7 7 9];
x=x*0.05-0.35;
y=[0 0 2 2 -3 -3 1 1 0 0];
Skalowanie=3.5;
h=plot(x,y*Skalowanie+15);
set(h,'LineWidth',2);
h=text(0.15,15,'Stim');
set(h,'FontSize',FontSize);
set(h,'FontName',FontName);

x=[0 1 1 9 9 10];
x=x*0.05-0.4;
y=[1 1 0 0 1 1];
h=plot(x,y*5+32);
set(h,'LineWidth',2);
h=text(0.15,34.5,'Rec');
set(h,'FontSize',FontSize);
set(h,'FontName',FontName);

h=gca;
set(h,'LineWidth',BoxLineWidth);
set(h,'XLim',[-0.5 2.5]);
hc = figure(FigureNumberSpike);
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[7 5]);
set(hc,'PaperPosition',[0 0 7 5]);
print(hc, '-dtiff', '-r400', 'C:\home\pawel\nauka\Stimchip_paper\obrazki\Spike');

%6) Plotting both EIs
FigureProperties=struct('FigureNumber',3,'TimeRange',TimeRange,'AmplitudeRange',[-150 100],'FontSize',FontSize,'FontName',FontName,'Colors',['k' 'r' 'k' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-100 -50 0 50 100]);
OffsetSamples=[1:-TimeRange(1)-1];
OffsetCancellation=2;
OffsetSamples=[1:-TimeRange(1)-4];
h1=NS_PlotDifferentTracesOnArrayLayout(T/0.84,Channels,OffsetCancellation,OffsetSamples,FigureProperties,NS_GlobalConstants);

hold on;
%break;
%h=gca;
axes(h1(4));
x=[0 1 1 3 3 5 5 7 7 8];
x=x*0.05-0.35;
y=[0 0 2 2 -3 -3 1 1 0 0];
h=plot(x,y*8+50);
set(h,'LineWidth',2);

h=gca;
set(h,'XLim',[-0.5 2.5]);
hc = figure(3);
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[7.3 5.3]);
set(hc,'PaperPosition',[0 0 7.3 5.3]);
print(hc, '-dtiff', '-r400', 'C:\home\pawel\nauka\Stimchip_paper\obrazki\Spike_EIs');