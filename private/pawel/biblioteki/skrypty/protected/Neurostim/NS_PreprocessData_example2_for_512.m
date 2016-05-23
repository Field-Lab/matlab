ChipAddresses=[24:31];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

cd D:\Home\Data\testy\2016-04-29-0;
FileName='000';
WritePath='D:\Home\Data\testy\2016-04-29-0-preproc';

WritePathFigs=WritePath;
Channels=[1:512];

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
ArrayID=1502;
%[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3_519(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[13],[1:512],[1:512],0);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main519(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:24],[1:512],[1:512],0);