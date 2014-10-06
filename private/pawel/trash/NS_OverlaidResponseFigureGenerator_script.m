ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

cd /Volumes/Bore/Data/Lauren/2008-11-12-0/;

FileName='001';
WritePath='/Volumes/Bore/Analysis/Lauren/2008-11-12-0/data001';
WritePathFigs='/Volumes/Bore/Analysis/Lauren/2008-11-12-0/data001Figs';
%WritePath='/Volumes/Lee/Analysis/Lauren/testPlotter';

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants);