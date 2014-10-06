ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
%NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
%ArrayID=500;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArrayID=1;

cd H:\2010-09-21-0;
FileName='017';
WritePath='E:\pawel\analysis\retina\2010-09-21-0\data017_2';
WritePathFigs=WritePath;

Channels=[1:64];
Movies=[1:26];
Movies=2+([1:10]-1)*5

AmplitudeRange=[-550 -250];
GoodChannels=NS_RemoveBadChannels(Channels,[4 9 25 31 57]);
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',AmplitudeRange,'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3try(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,1,Movies,[1:64],[1:64]);

FileName='020';
WritePath='E:\pawel\analysis\retina\2010-09-21-0\data020_2';
WritePathFigs=WritePath;
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3try(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,1,Movies,[1:64],[1:64]);