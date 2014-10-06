ChipAddresses=[31 30];
NumberOfChannelsPerChip=32;

CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

% 1. Common variables
Channel=40;
MovieNumber=31;
Channels=Channel;
Channels=40;
TimeStart=10;
NumberOfSamples=30;
t=[1:NumberOfSamples];
Offsets=ones(length(Channels))*(-370);

% 2. Artifact calculation
%[CI,MI]=NS_Report(FileName,1,NS_GlobalConstants);
cd D:\2008-03-18saline;
ArtifactFileName='000';
%[CI,MI]=NS_Report(ArtifactFileName,1,NS_GlobalConstants);
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(ArtifactFileName,Channel,MovieNumber,NS_GlobalConstants);
Traces=NS_ReadManyTracesFromRaw(ArtifactFileName,Channels,Timings,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);
Artifact=reshape(mean(Traces),1,NumberOfSamples);
figure(2);
plot(Artifact);

[Pulse,Status]=NS_FindPulseShapeForMovie(ArtifactFileName,Channel,MovieNumber,NS_GlobalConstants);
figure(1);
clf;
subplot(2,1,1);
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,Channel,0,1,'b-',NS_GlobalConstants);

% 3. Data reading
cd D:\2008-03-21-0;
DataFileName='003';
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(DataFileName,Channel,MovieNumber,NS_GlobalConstants);
Traces=NS_ReadManyTracesFromRaw(DataFileName,Channels,Timings,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);
Data=reshape(Traces,98,NumberOfSamples);
figure(3);
plot(t,Data');

[Pulse1,Status1]=NS_FindPulseShapeForMovie(DataFileName,Channel,MovieNumber,NS_GlobalConstants);
figure(1);
subplot(2,1,2);
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse1,Status1,Channel,0,1,'b-',NS_GlobalConstants);

% 4. Artifact subtraction
Data1=Data;
for i=1:98
    Data1(i,:)=Data(i,:)-Artifact;
end
figure(4);
plot(t,Data1');