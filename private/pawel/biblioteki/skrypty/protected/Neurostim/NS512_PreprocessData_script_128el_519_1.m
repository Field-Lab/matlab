ChipAddresses=[24:31];
ChipAddresses=[31:-1:24];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

ArrayID=1502;

cd /Volumes/2012Monkey2/2012-09-27-4;
FileName='003';
%WritePath='C:\home\data\2012-09-16-4-test';
WritePath = '/Volumes/Stream-phoenix/Analysis/stim512/2012-09-27-4/stim_scan/';
WritePathFigs=WritePath;

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);

tic
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],[1:512],[1:512],0);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],[1:512],[1:512],-1);
toc