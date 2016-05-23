ChipAddresses=[24:31];
ChipAddresses=[31:-1:24];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

cd I:\analysis\slices\2013-12-12-3-PH;
FileName='001';
WritePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';
WritePathFigs=WritePath;

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);

MoviesToAnalyze=[2:451];
%MoviesToAnalyze=[382:451];
PatternsToAnalyze=NS512_AllPatternsInExperiment('I:\analysis\slices\2013-12-12-3-PH\movie001',NS_GlobalConstants);

[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,MoviesToAnalyze,PatternsToAnalyze',[1:512],0);