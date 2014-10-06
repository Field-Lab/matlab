ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;

cd F:\2008-02-19\2008-02-19-1;
FileName='013';
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

%[CI,MI]=NS_Report(FileName,1,NS_GlobalConstants);

MovieNumber=17;
TimeRange=[-10 170];
Limit=100;
RecChannel=13;
Channels=[7 10 11 13 14 16 18 19];

[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,11,MovieNumber,NS_GlobalConstants);
figure(51);
subplot(4,1,1);
Amplitude=NS_PlotStimulationPulsesTrain(Pulse,Status,11,100,400,4,20,'r-',NS_GlobalConstants);
gca;
axis([-5 75 -0.6 0.6]);
%grid off;
[Pulse2,Status]=NS_FindPulseShapeForMovie(FileName,16,MovieNumber,NS_GlobalConstants);
subplot(4,1,2);
Amplitude=NS_PlotStimulationPulsesTrain(Pulse2,Status,16,100,400,4,20,'b-',NS_GlobalConstants);
gca;
axis([-15 65 -1 1]);
%grid off;

Channel=11;
[Timings11,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);
Timings=Timings11(1,1:Limit);
WaveformTypes([13 14 17 25 26 27 32 33 34 39 40 44 47 49 52 53 54 57 59 60 65 66 67 72 73 77 78 86 87 88 91 92 93 94 95 98 99])=0; %channel 11
a11=[13 14 17 25 26 27 32 33 34 39 40 44 47 49 52 53 54 57 59 60 65 66 67 72 73 77 78 86 87 88 91 92 93 94 95 98 99]; %channel 11
[Rawtraces,artifact11]=NS_AverageTraces([pwd '\data' FileName],Timings(a11),Channels,TimeRange,NS_GlobalConstants); % artifact related tu pulse on channel 11
artifact11=artifact11-mean(artifact11(1:-TimeRange(1)-2));
figure(11);
plot(artifact11);

Channel=16;
[Timings16,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);
Timings=Timings16(1,1:Limit);
a16=[47 48 61 78 79 80 88 89 90 91 92 93];
[RAWtraces,artifact16]=NS_AverageTraces([pwd '\data' FileName],Timings(a16),Channels,TimeRange,NS_GlobalConstants); % artifact related tu pulse on channel 11
artifact16=artifact16-mean(artifact16(1:-TimeRange(1)-2));
figure(16);
plot(artifact16);

Channel=60;
[Timings60,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,Channel,MovieNumber,NS_GlobalConstants);
Timings=Timings60(1,1:Limit);
a=[1:100];
[RAWtraces,artifact60]=NS_AverageTraces([pwd '\data' FileName],Timings(a),Channels,TimeRange,NS_GlobalConstants); % artifact related tu pulse on channel 11
artifact60=artifact60-mean(artifact60(1:-TimeRange(1)-2));
figure(60);
plot(artifact60);

%RecChannel=11;
TimeStart=14210100;
TimeLength=100000;
t=[1:TimeLength]/20;
%rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile('D:\2008-02-19\2008-02-19-1\data013');
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile([pwd '\data' FileName]);
data = rawFile.getData(TimeStart+1,TimeLength);
SignalTrace=zeros(1,TimeLength);
SignalTrace=double(data(:,RecChannel+1)');
clear data;
figure(1)
plot(t,SignalTrace);

Timings11a=Timings11+TimeRange(1)-TimeStart;
%Timings11=Timings11(1,1:400);
ble=100
SignalTrace=NS_SubtractArtifact_StimchipPaper(SignalTrace,artifact11(1:ble,5)'+17,Timings11a(1,2:24));
Timings16a=Timings16+TimeRange(1)-TimeStart;
SignalTrace=NS_SubtractArtifact_StimchipPaper(SignalTrace,artifact16(1:ble,5)'+13,Timings16a(1,2:24));
Timings60a=Timings60+TimeRange(1)-TimeStart;
SignalTrace=NS_SubtractArtifact_StimchipPaper(SignalTrace,artifact60(1:ble,5)'+13,Timings60a(1,2:24));

S=SignalTrace(1,9801:12000);
figure(53);
clf;
hold on;

x=[-1 1 1 3 3 5 5 7 7 9];
y=[0 0 2 2 -3 -3 1 1 0 0];
Skalowanie=4;

SamplesToShow=max(x)+T1*20:min(x)+3+T2*20;

for i=1:4
    T1=5+(i-1)*20;
    T2=15+(i-1)*20;
    
    SamplesToShow=max(x)+T1*20:min(x)+2+T2*20;
    S1=S(1,SamplesToShow)+380;
    h=plot(SamplesToShow/20-0.5,S1);
    set(h,'Color','k');
    
    SamplesToShow=max(x)+(T1+10)*20:min(x)+2+(T2+10)*20;
    S1=S(1,SamplesToShow)+380;
    h=plot(SamplesToShow/20-0.5,S1);
    set(h,'Color','k');
                   
    x1=x*0.05-0.35+T1;        
    h=plot(x1,y*Skalowanie+40);
    set(h,'LineWidth',2);
    set(h,'Color','r');  
    
    x1=x*0.05-0.35+T2;        
    h=plot(x1,y*Skalowanie+40);
    set(h,'LineWidth',2);
    set(h,'Color','b');  
end
grid on;
h=gca;
%set(h,'XTick',[0:5:80])
set(h,'FontSize',20);
xlabel('Time [ms]');
ylabel('Amplitude [\muV]');
axis([0 80 -70 60]);


FigureNumber=22;
FigureProperties=struct('FigureNumber',FigureNumber,'TimeRange',TimeRange,'AmplitudeRange',[-800 -400],'FontSize',13,'Colors',['k' 'r' 'b' 'y']);

size(Timings);

Delay=0;
WaveformTypes=ones(1,length(Timings));

%WaveformTypes([13 14 17 25 26 27 32 33 34 39 40 44 47 49 52 53 54 57 59 60 65 66 67 72 73 77 78 86 87 88 91 92 93 94 95 98 99])=0; %channel 11
%WaveformTypes([47 48 61 78 79 80 88 89 90 91 92 93])=0; %channel 16;
WaveformTypes([])=0;

y=NS_PlotManyTracesOnFigure_StimchipPaper(FileName,Timings,WaveformTypes,Channel,2,[1:10],FigureProperties,Delay,NS_GlobalConstants);

break;
%y=NS_PlotManyTracesOnFigure(FileName,Timings,WaveformTypes,Channel,2,[1:10],FigureProperties,Delay,NS_GlobalConstants);
%break;
%Channels=[58:60]; %for channel 11
%Channels=[1:64];
a=find(WaveformTypes==1);
[R,signal1]=NS_AverageTraces([pwd '\data' FileName],Timings(a),Channels,TimeRange,NS_GlobalConstants);
a=find(WaveformTypes==0);
[R,signal2]=NS_AverageTraces([pwd '\data' FileName],Timings(a),Channels,TimeRange,NS_GlobalConstants);

FigureProperties=struct('TimeRange',TimeRange,'AmplitudeRange',[-400 200],'FontSize',11,'Colors',['g' 'b' 'r' 'y']);
OffsetCancellation=2;
OffsetSamples=[1:-TimeRange(1)-1];
WaveformTypes=ones(1,length(Timings));
WaveformTypes(1,Channel)=2;
WaveformTypes([9 25 57 41 64 28 37])=0; %noisy electrodes
%y=NS_PlotTracesOnArrayLayout(signal1',Channels,WaveformTypes,OffsetCancellation,OffsetSamples,1,32,FigureProperties,NS_GlobalConstants);

figure(66);
clf;
Limit2=[45:94];
exclude=73;
Limit2=[45:exclude-1 exclude+1:95];
%Limit2=[1:100];
ArtifactCancellation=1;

Timings=Timings11(1,Limit2);
WaveformTypes=ones(1,Limit);
WaveformTypes(a11)=0;
WaveformTypes=WaveformTypes(1,Limit2);
ArrayID=1;
%TimeRange=[-10 60];
FigureProperties=struct('FigureNumber',FigureNumber,'TimeRange',TimeRange,'AmplitudeRange',[-200 100],'FontSize',20,'Colors',['k' 'r' 'b' 'y']);
h=NS_PlotManyTracesOnArrayLayoutOld(FileName,Channels,Timings,WaveformTypes,ArtifactCancellation,ArrayID,66,FigureProperties,NS_GlobalConstants,3);
%break;
Timings=Timings16(1,Limit2);
WaveformTypes=ones(1,Limit)*2;
WaveformTypes(a16)=0;
WaveformTypes=WaveformTypes(1,Limit2);
h=NS_PlotManyTracesOnArrayLayoutOld(FileName,Channels,Timings,WaveformTypes,ArtifactCancellation,ArrayID,66,FigureProperties,NS_GlobalConstants,2);

break
name='C:\home\pawel\2010\Stimchip_paper\figures_two_neurons\TwoNeurons_AllChannels';
hc=figure(66);
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[10 9]);
set(hc,'PaperPosition',[0 0 10 9]);
print(hc, '-dtiff', '-r120', name);