ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
%NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
%ArrayID=500;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArrayID=1;

cd J:\2010-09-21-0;
FileName='003';
WritePath='D:\Home\Pawel\analysis\2010-09-21-0\data003';
WritePathFigs=WritePath;

WritePathFigs=WritePath;
Channels=[1:64];
Movies=[1:26];

AmplitudeRange=[-550 -250];
GoodChannels=NS_RemoveBadChannels(Channels,[4 9 25 31 57]);
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',AmplitudeRange,'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v4(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,Movies,[1:64],[1:64]);

FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
DataPath=WritePath; %sciezka do danych po preprocessingu

%usun szumiace kanaly:
BadChannels=[4 9 25 57];
Channels=[1:64];
Channels=NS_RemoveBadChannels(Channels,BadChannels);

Patterns=Channels; %ktore patterny (zestawy elektrod) chcemy przeanalizowac, domyslnie [1:64]

WritePathFigsGood='D:\Home\Pawel\analysis\2010-09-21-0\proba\figures_good';
WritePathFigsBad='D:\Home\Pawel\analysis\2010-09-21-0\proba\figures_bad';

% !!!!!!!!

MinDelay=6; %valid for 50 us pulses, for 100 us should be 6 

% !!!!!!!!!!!

for i=1:length(Patterns)
    PatternNumber=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,1,[1:64],NS_GlobalConstants);
    [AmplitudesVsChannels,RMSVsChannels]=NS512_FindThresholdForLocalStimulation4PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,StimChannels,[MinDelay:40]); %7 to 76    
end