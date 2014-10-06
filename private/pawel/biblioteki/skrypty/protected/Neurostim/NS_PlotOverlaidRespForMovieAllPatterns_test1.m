ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

cd D:\2008-08-27-4;

FileName='005';
WritePath='D:\analysis\2008-08-27-4\data005';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-1200 200],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,FigureProperties,NS_GlobalConstants);

FileName='014';
WritePath='D:\analysis\2008-08-27-4\data014';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-1200 200],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,FigureProperties,NS_GlobalConstants);

FileName='002';
WritePath='D:\analysis\2008-08-27-4\data002';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-1200 200],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,FigureProperties,NS_GlobalConstants);

FileName='013';
WritePath='D:\analysis\2008-08-27-4\data013';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-1200 200],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,FigureProperties,NS_GlobalConstants);
