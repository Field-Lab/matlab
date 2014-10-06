ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

ChipAddresses=[24:31];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

cd I:\2010-09-05-test;
FileName='015';
WritePath='I:\analysis\2010-09-05-test\015';
WritePathFigs='I:\analysis\2010-08-20-0\012data\figures';

WritePathFigs=WritePath;
Channels=[1:64];

AmplitudeRange=[-100 100];
GoodChannels=NS_RemoveBadChannels(Channels,[4 9 25 31 57]);
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
%[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,1,[1:38],[5 10 17 19 27],[1:64]);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,[1:2],[1:512],[1:512]);
break;

PatternNumber=22;
Movies=[7:1:25];
Channels=[1:64];
GoodChannels=NS_RemoveBadChannels(Channels,[4 9 25 31 57]);
AdditionalBadChannels=[];
GoodChannels=NS_RemoveBadChannels(GoodChannels,AdditionalBadChannels);
ClusterFileName='C:\home\pawel\nauka\analysis\2008-12-11-0\ClusterFile_000';
Responses=[];
tic
for PatternNumber=GoodChannels
    PatternNumber;
    a=NS_ApplyLinearArtifactModel(WritePath,PatternNumber,Movies,GoodChannels,0,0,ClusterFileName);
    Responses=[Responses' a']';
end
toc