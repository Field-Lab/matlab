ChipAddresses=[24:31];
ChipAddresses=[31:-1:24];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

%cd I:\analysis\slices\2010-09-14-0\TTX_sub;
cd G:\analysis\slices\2013-12-12-0-PH\;
FileName='001';
WritePath='J:\analysis\2013-12-12-0\data001preproc';
WritePathFigs=WritePath;
AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
MoviesToAnalyze=[2:451];
PatternsToAnalyze=NS512_AllPatternsInExperiment('J:\data\2013-12-12-0-PH\movie001',NS_GlobalConstants)
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,MoviesToAnalyze,PatternsToAnalyze',[1:512],0);

FileName='004';
WritePath='J:\analysis\2013-12-12-0\data004preproc';
WritePathFigs=WritePath;
AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
MoviesToAnalyze=[2:451];
PatternsToAnalyze=NS512_AllPatternsInExperiment('J:\data\2013-12-12-0-PH\movie004',NS_GlobalConstants)
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,MoviesToAnalyze,PatternsToAnalyze',[1:512],0);

cd G:\analysis\slices\2013-12-15-0;
FileName='004';
WritePath='J:\analysis\2013-12-15-0\data004';
WritePathFigs=WritePath;
AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
MoviesToAnalyze=[2:451];
PatternsToAnalyze=NS512_AllPatternsInExperiment('J:\data\2013-12-15-0\movie004',NS_GlobalConstants)
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,MoviesToAnalyze,PatternsToAnalyze',[1:512],0);

FileName='008';
WritePath='J:\analysis\2013-12-15-0\data008';
WritePathFigs=WritePath;
AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
MoviesToAnalyze=[2:451];
PatternsToAnalyze=NS512_AllPatternsInExperiment('J:\data\2013-12-15-0\movie008',NS_GlobalConstants)
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,MoviesToAnalyze,PatternsToAnalyze',[1:512],0);

clear
Stim512Paper_DetectAllSpikes2.m;

break
FileName='004';
WritePath='I:\analysis\slices\2013-12-12-0-PH\data004preproc';
WritePathFigs=WritePath;
AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
MoviesToAnalyze=[2:451];
PatternsToAnalyze=NS512_AllPatternsInExperiment('I:\analysis\slices\2013-12-12-0-PH\movie004',NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,MoviesToAnalyze,PatternsToAnalyze',[1:512],0);

break

MovieFile='I:\analysis\slices\2013-12-12-0-PH\movie001';
DataPath='I:\analysis\slices\2013-12-12-0-PH\data001preproc';
ArtifactDataPath=DataPath;
FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-0-PH\SpikesAnalysis\data001';
SzybkiTestDetekcjiSpikow_v5.m;

MovieFile='I:\analysis\slices\2013-12-12-0-PH\movie004';
DataPath='I:\analysis\slices\2013-12-12-0-PH\data004preproc';
ArtifactDataPath=DataPath;
FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-0-PH\SpikesAnalysis\data004';
SzybkiTestDetekcjiSpikow_v5.m;