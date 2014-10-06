ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-1200 450],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
Patterns=[1];
PatternsIndexes=[2];

DataPath='D:\analysis\2008-08-26-0\data006_proba';
ArtifactDataPath='D:\analysis\2008-08-26-0\data011_proba';

PatternNumber=11;
MovieNumber=55;

figure(FigureProperties.FigureNumber);
clf;
figure(FigureProperties.FigureNumber+10);
clf;
[Channels]=NS_ShowPreprocessedData(DataPath,ArtifactDataPath,PatternNumber,MovieNumber,Patterns,PatternsIndexes,FigureProperties,NS_GlobalConstants);